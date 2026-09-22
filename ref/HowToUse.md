# Input file naming and required columns (Tercen)

Quick reference for naming the files you upload from Tercen for the autoreport app.


## QC

**In a report with only either TRs or BRs, name the QC export(s) as:**

```
QC_<PTK|STK>_[<number>]_<TR|BR>.csv
```

* `number`: use different number if there are multiple QC exports from different normalizations (e.g. log & VSN). This is rarely used.

**In a report with both TRs and BRs:**

- first export the TRs and add the factor defining biological replicate (e.g. Biol_Rep or Sample name). 

- Then TRs are averaged to BRs. Export mean BR values.

  - Barcode + Row --&gt;Biol_Rep or Sample name / sample id (a column containing "sample")

  - exported value: ".value" 


Name the QC export as:

```
QC_<PTK|STK>_[<number>]_<TR|BR>[_<Log|LogCmb|VSN|VSNCmb>].csv
```
- `Log` / `LogCmb` / `VSN` / `VSNCmb` : case-insensitive. Since averaging results in ".value" column, need to define what were the original values.


**Required columns**:

| Type | Column |
|---|---|
| Sample |`Barcode`+`Row`, or a column name containing "sample" (case-insensitive), or a column named exactly `Biol_Rep` |
| ID | `ID` |
| value column(s) | any of: `.logTransformed`, `.identity`, `.CmbCor`, `.value` |
| Conditions | any of: `Test Condition`, `Supergroup`, Sample name, Biol_Rep or other column(s) |

## Limma

```
Limma_<PTK|STK>_<number>_<dataset name>.csv
```
- `Dataset name`: Any name. Must not contain underscores. Must be identical between the PTK and STK files. e.g. TvsC, or WP1. 
- Supergroup names and comparisons come from the file content (contrast and optionally Supergroup column), not the filename.

## UKA

```
UKA_<PTK|STK>_<number>_<Dataset name>.csv
```

- `Dataset name`: Any name. Must not contain underscores. Must be identical between the PTK and STK files. e.g. TvsC, or WP1

**Required columns**:

- `Sgroup_contrast` - this is needed to extract comparisons!
