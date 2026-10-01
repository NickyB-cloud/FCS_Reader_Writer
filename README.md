# FCS_Reader_Writer

*by Nicholas B. Franks, PhD*

**Have you ever spent half a month in tears because FlowIO and FlowKit keep breaking your scripts?**

Have you ever upgraded Python and watched `pip install flowkit` die while trying to compile FlowUtils? Pinned three packages to three versions that refuse to agree? Opened a Cytek file only to be told its TEXT segment is "malformed" by software that was fine yesterday?

Me too. So here is one file that reads and writes `.fcs` files and needs nothing but **numpy** and **pandas**. No FlowIO, no FlowKit, no FlowUtils, no compiler, no tears.

It's the reader/writer I built into my own spectral flow cytometry pipeline (Cytek Aurora data), pulled out into one file you can drop next to your script.

---

## Install

There's nothing to install except numpy and pandas:

```bash
pip install numpy pandas
```

Then download `FCS_Reader_Writer.py` and put it in the same folder as your script. That's it.

(Tested on Python 3.14 with numpy 2.4.4 and pandas 2.3.3. See `requirements.txt`.)

---

## Use it in Python

```python
from FCS_Reader_Writer import load_fcs_as_dataframe, read_fcs, write_fcs_3_0

# Read: one column per channel, linear values
df = load_fcs_as_dataframe("sample.fcs")
print(df.columns.tolist())      # ['FSC-A', 'SSC-A', 'V1-A', ...]

# Or read the raw pieces
events, channels, meta = read_fcs("sample.fcs")
# events   -> numpy array, shape (n_events, n_channels)
# channels -> list of channel names
# meta     -> dict of every TEXT keyword (keys are lower-case, e.g. meta["$tot"])

# Do something with it
gated = df[(df["FSC-A"] > 50_000) & (df["SSC-A"] < 200_000)]

# Write it back out as a normal .fcs file
write_fcs_3_0("sample_gated.fcs", gated)
```

## Use it from the command line

```bash
python FCS_Reader_Writer.py sample.fcs                 # events, channels, min/median/max per channel
python FCS_Reader_Writer.py sample.fcs --keywords      # ...plus every TEXT keyword
python FCS_Reader_Writer.py sample.fcs --csv out.csv   # ...plus export all events to CSV
```

---

## What it handles

**Reading**
- FCS 3.0, 3.1 and 3.2 list-mode files
- `$DATATYPE` I (8/16/32/64-bit), F (float32) and D (float64)
- Both byte orders (`$BYTEORD` 1,2,3,4 and 4,3,2,1)
- Keywords with or without the leading `$` (Cytek writes some of each)
- Any TEXT delimiter, including escaped (doubled) delimiters
- Large files (>100 MB) whose data offsets live in `$BEGINDATA` / `$ENDDATA`
- Channel names from `$PnN`, falling back to `$PnS`, then `Ch1`, `Ch2`, ...

**Writing**
- FCS 3.0, float32, little-endian, from any pandas DataFrame
- Column names become `$PnN` (a `|` in a name becomes `/`, since `|` is the delimiter)
- Handles files over 100 MB

## What it does *not* do

This reader is deliberately simple. It returns values **exactly as stored** (linear) and does **not** apply:
- `$PnE` log decoding or `$PnG` gain
- compensation (`$SPILLOVER`)
- arcsinh, logicle or any other transform
- gating, FlowJo workspaces, or anything else FlowKit does

Files with a different bit depth on each channel are refused with an error rather than read wrong.

If you need those features, FlowKit is the right tool. If you just need your data in a DataFrame and back out again, this is the easier route.

---

## Does it actually work?

Before release I checked that:
- the writer's output is **byte-for-byte identical** to the writer in my analysis pipeline, on both a small file and a 104 MB file
- the reader gives **identical** arrays, channel names and keywords to my pipeline's reader on real Cytek Aurora raw and unmixed files
- the reader correctly reads test files with 8/16/32-bit integer, float32 and float64 data in both byte orders

### Proof it works

This reader/writer ran the entire spectral flow cytometry pipeline behind my PhD dissertation: raw Cytek Aurora files through gating, unmixing and final figures. You can find the dissertation through the **Tufts University Tisch Library** by searching:

> Nicholas Franks, Impact of Microbial Modification of Lipid A on Dendritic Cell Functional Responses

---

## License and citation

MIT License. Free to use, modify and share; just keep the copyright notice. See `LICENSE`.

If this saved you some tears, please cite it. GitHub's **"Cite this repository"** button (from `CITATION.cff`) gives a ready-made citation:

> Franks, N. B. *FCS_Reader_Writer: a standalone pure-Python FCS file reader and writer* (v1.0.0).
