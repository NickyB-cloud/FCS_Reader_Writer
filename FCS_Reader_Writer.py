#!/usr/bin/env python3
"""
FCS_Reader_Writer.py — Standalone pure-Python FCS reader/writer (numpy + pandas only).

Author: Nicholas B. Franks, PhD

This is the FCS code built into my spectral flow cytometry pipeline, pulled
out into one importable module so you can use it in your own work without
installing an FCS library.

Reads  : FCS 3.0 / 3.1 / 3.2 list-mode files, DATATYPE I (8/16/32/64-bit),
         F (float32) or D (float64), either byte order. Tested on Cytek
         Aurora-style spectral files. Values are returned LINEAR as stored
         (no $PnE log decoding, no compensation, no scaling).
Writes : FCS 3.0, float32, little-endian, Cytek-style keyword order.

Use from Python:
    from FCS_Reader_Writer import load_fcs_as_dataframe, read_fcs, write_fcs_3_0

    df = load_fcs_as_dataframe("sample.fcs")          # pandas DataFrame
    events, channels, meta = read_fcs("sample.fcs")   # numpy array, names, TEXT dict
    write_fcs_3_0("sample_copy.fcs", df)              # DataFrame -> .fcs

Use from the command line:
    python FCS_Reader_Writer.py sample.fcs                   # summary: events, channels, keywords
    python FCS_Reader_Writer.py sample.fcs --csv out.csv     # also export events to CSV
    python FCS_Reader_Writer.py sample.fcs --keywords        # print every TEXT keyword
"""

import argparse
import os
import warnings

import numpy as np
import pandas as pd


# ─────────────────────────────────────────────────────────────
# READER
# ─────────────────────────────────────────────────────────────

def parse_text_segment(text_block):
    """
    Parse the TEXT segment of an FCS file into a dict.

    Handles arbitrary delimiters, doubled (escaped) delimiters, and
    Cytek-style keywords with or without the leading "$". Keys are
    lower-cased; values are stripped strings.
    """
    if not text_block:
        return {}

    # First character is the delimiter
    delim = text_block[0]

    # FCS 3.1 spec: a doubled delimiter is a literal delimiter inside a value
    placeholder = "\x00"
    text_block = text_block.replace(delim * 2, placeholder)
    parts = text_block.strip(delim).split(delim)
    parts = [p.replace(placeholder, delim) for p in parts]

    # Odd number of tokens -> last value would be silently dropped
    if len(parts) % 2 != 0:
        warnings.warn(
            f"TEXT segment has odd number of tokens ({len(parts)}); "
            "last token will be ignored — file may be malformed.",
            UserWarning,
        )

    # key1, val1, key2, val2, ...
    text_dict = {}
    for i in range(0, len(parts) - 1, 2):
        key = parts[i].strip().lower()
        val = parts[i + 1].strip()
        text_dict[key] = val

    return text_dict


def read_fcs(path):
    """
    Read an FCS file.

    Returns
    -------
    events        : numpy array, shape (n_events, n_channels), linear values
    channel_names : list of str ($PnN, falling back to $PnS, then "Ch<i>")
    meta          : dict of every TEXT keyword (lower-case keys)
    """
    with open(path, "rb") as f:
        raw = f.read()

    # ------------------------------
    # HEADER
    # ------------------------------
    header = raw[:58].decode("ascii", errors="ignore")
    fcs_type = header[:6].strip()
    if not fcs_type.startswith("FCS"):
        raise ValueError(f"Not a valid FCS file (got header: {fcs_type!r})")

    text_start = int(header[10:18].strip() or 0)
    text_end   = int(header[18:26].strip() or 0)
    data_start = int(header[26:34].strip() or 0)
    data_end   = int(header[34:42].strip() or 0)

    # ------------------------------
    # TEXT SEGMENT
    # ------------------------------
    text_block = raw[text_start:text_end + 1].decode("ascii", errors="ignore")
    t = parse_text_segment(text_block)

    # Read both bare and $-prefixed keys
    def meta(key, default=None):
        return t.get(key, t.get(f"${key}", default))

    # TEXT overrides header offsets (header holds 0 for files > 99,999,999 bytes)
    data_start = int(meta("begindata", data_start))
    data_end   = int(meta("enddata",   data_end))

    # Fallback: scan for first binary byte
    if data_start == 0:
        for i in range(text_end + 1, len(raw)):
            if raw[i] not in b"0123456789.-+eE \t\r\n":
                data_start = i
                break

    par = int(meta("par", 0))
    tot = int(meta("tot", 0))

    data = raw[data_start:data_end + 1]

    # ------------------------------
    # NUMPY DTYPE FROM $DATATYPE / $BYTEORD / $PnB
    # ------------------------------
    dtype = str(meta("datatype", "I")).upper()
    byteord = str(meta("byteord", "1,2,3,4"))
    little_endian = byteord.startswith("1")
    endian = "<" if little_endian else ">"

    if dtype == "I":
        bits_per_channel = [int(meta(f"p{i}b", 16)) for i in range(1, max(par, 1) + 1)]
        if par > 0 and len(set(bits_per_channel)) > 1:
            raise ValueError(f"Mixed per-channel bit depths: {bits_per_channel}")

        bits = bits_per_channel[0] if par > 0 else 32
        if bits == 8:
            dt = np.dtype(f"{endian}u1")
        elif bits == 16:
            dt = np.dtype(f"{endian}u2")
        elif bits == 32:
            dt = np.dtype(f"{endian}u4")
        elif bits == 64:
            dt = np.dtype(f"{endian}u8")
        else:
            raise ValueError(f"Unsupported integer bit depth: {bits}")

    elif dtype == "F":
        dt = np.dtype(f"{endian}f4")

    elif dtype == "D":
        dt = np.dtype(f"{endian}f8")

    else:
        raise ValueError(f"Unsupported DATATYPE: {dtype}")

    arr = np.frombuffer(data, dtype=dt)

    # ------------------------------
    # INFER $PAR IF MISSING
    # ------------------------------
    if par == 0:
        print(f"⚠ {os.path.basename(path)} — $PAR missing, inferring from DATA segment")

        common_par_values = [40, 48, 50, 64, 67, 71, 73, 80]
        inferred = None
        for guess in common_par_values:
            if len(arr) % guess == 0:
                inferred = guess
                break

        if inferred is None:
            raise ValueError(f"Cannot infer PAR for {path} — file may be malformed.")

        par = inferred
        print(f"  → Inferred PAR = {par}")

    # ------------------------------
    # CHANNEL NAMES (Aurora-safe)
    # ------------------------------
    channel_names = []
    for i in range(1, par + 1):
        candidates = [
            meta(f"p{i}n"),
            meta(f"p{i}s"),
            meta(f"p{i}display"),
            meta(f"p{i}type"),
            meta(f"p{i}longname"),
            meta(f"p{i}shortname"),
        ]
        name = next((c for c in candidates if c and c.strip()), None)
        if not name:
            name = f"Ch{i}"
        channel_names.append(name)

    # ------------------------------
    # RESHAPE EVENTS
    # ------------------------------
    if tot > 0 and len(arr) == tot * par:
        events = arr.reshape((tot, par))
    else:
        n = len(arr) // par
        if n == 0:
            raise ValueError("DATA segment too small to contain any events.")
        events = arr[: n * par].reshape((n, par))

    return events, channel_names, t


def load_fcs_as_dataframe(path):
    """Read an FCS file into a pandas DataFrame (one column per channel, linear values)."""
    events, channel_names, _ = read_fcs(path)
    return pd.DataFrame(events, columns=channel_names)


# ─────────────────────────────────────────────────────────────
# WRITER
# ─────────────────────────────────────────────────────────────

def write_fcs_3_0(filename, df):
    """
    Write a DataFrame to an FCS 3.0 file (float32, little-endian).

    Column names become $PnN ("|" is replaced by "/" since "|" is the
    TEXT delimiter). $PnR is each column's max + 1. Files whose DATA
    segment ends past byte 99,999,999 get offsets in $BEGINDATA/$ENDDATA
    with zeros in the header, as the FCS 3.0 spec requires.
    """
    data = df.values.astype("float32")
    n_events, n_channels = data.shape
    channels = df.columns.tolist()

    def sanitize(name):
        return str(name).replace("|", "/")

    # ---- TEXT keywords in Cytek Aurora ordering ----
    text_items = []

    # Global keys
    text_items.append(("beginanalysis", "0"))
    text_items.append(("endanalysis", "0"))
    text_items.append(("byteord", "1,2,3,4"))
    text_items.append(("datatype", "F"))
    text_items.append(("mode", "L"))
    text_items.append(("nextdata", "0"))
    text_items.append(("tot", str(n_events)))
    text_items.append(("par", str(n_channels)))

    # Per-parameter keys
    for i, ch in enumerate(channels, start=1):
        safe = sanitize(ch)
        ch_max = df[ch].max()
        rng = str(int(ch_max) + 1) if np.isfinite(ch_max) else "4194304"

        text_items.append((f"p{i}b", "32"))
        text_items.append((f"p{i}e", "0,0"))
        text_items.append((f"p{i}n", safe))
        text_items.append((f"p{i}r", rng))
        text_items.append((f"p{i}display", "LOG"))
        text_items.append((f"p{i}type", "Raw_Fluorescence"))
        text_items.append((f"p{i}v", "0"))

    def build_text(extra_items=None):
        items = text_items.copy()
        if extra_items:
            items.extend(extra_items)
        txt = "|" + "|".join(f"{k}|{v}" for k, v in items) + "|"
        return txt.encode("ascii")

    # ---- DATA block (same bytes as struct.pack("<{n}f", ...), much faster) ----
    data_bytes = data.astype("<f4").tobytes()

    # ---- Two-pass offset resolution ----
    text_start = 256
    text_bytes = build_text()

    for _ in range(2):
        text_end = text_start + len(text_bytes) - 1
        data_start = text_end + 1
        data_end = data_start + len(data_bytes) - 1

        if data_end > 99_999_999:
            text_bytes = build_text([
                ("$BEGINDATA", str(data_start)),
                ("$ENDDATA", str(data_end)),
            ])
            hdr_data_start = "       0"
            hdr_data_end = "       0"
        else:
            hdr_data_start = f"{data_start:>8}"
            hdr_data_end = f"{data_end:>8}"

    # ---- HEADER ----
    header = (
        b"FCS3.0    "
        + f"{text_start:>8}".encode("ascii")
        + f"{text_end:>8}".encode("ascii")
        + hdr_data_start.encode("ascii")
        + hdr_data_end.encode("ascii")
        + f"{'0':>8}".encode("ascii")
        + f"{'0':>8}".encode("ascii")
        + b" " * (256 - 58)
    )

    with open(filename, "wb") as fh:
        fh.write(header)
        fh.write(text_bytes)
        fh.write(data_bytes)

    print(f"  ✓ Wrote FCS 3.0 file : {filename}")
    print(f"    Events : {n_events:,}  |  Channels : {n_channels}")


# ─────────────────────────────────────────────────────────────
# COMMAND LINE
# ─────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Inspect an FCS file (and optionally export it to CSV) "
                    "with the pure-Python reader.")
    ap.add_argument("fcs", help="path to an .fcs file")
    ap.add_argument("--csv", metavar="OUT", help="also write all events to this CSV")
    ap.add_argument("--keywords", action="store_true",
                    help="print every TEXT keyword and value")
    args = ap.parse_args()

    events, channels, meta = read_fcs(args.fcs)
    with open(args.fcs, "rb") as f:
        version = f.read(6).decode("ascii", errors="ignore")

    print(f"File      : {os.path.basename(args.fcs)}")
    print(f"Version   : {version}")
    print(f"Events    : {events.shape[0]:,}")
    print(f"Channels  : {events.shape[1]}")
    print(f"Data type : {events.dtype}")
    print()
    print(f"  {'#':>3}  {'name':<24} {'min':>14} {'median':>14} {'max':>14}")
    for i, ch in enumerate(channels):
        col = events[:, i]
        print(f"  {i + 1:>3}  {ch:<24} {col.min():>14.4g} {np.median(col):>14.4g} {col.max():>14.4g}")

    if args.keywords:
        print()
        for k, v in meta.items():
            print(f"  {k} = {v}")

    if args.csv:
        pd.DataFrame(events, columns=channels).to_csv(args.csv, index=False)
        print(f"\n  ✓ Wrote {args.csv}")


if __name__ == "__main__":
    main()
