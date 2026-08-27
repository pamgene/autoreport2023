# Input file naming and required columns (Tercen)

Quick reference for naming the files you upload from Tercen for the autoreport app.


## QC

In a report with only either TRs or BRs

```
QC_<PTK|STK>_[<number>]_<TR|BR>.csv
```

In a report with both TRs and BRs, where the TRs are averaged:
```
QC_<PTK|STK>_[<number>]_<TR|BR>[_<Log|LogCmb|VSN|VSNCmb>].csv
```
- `Log` / `LogCmb` / `VSN` / `VSNCmb` : case-insensitive

### Required columns:

| Column | Type |
|---|---|
| `Barcode`, `Row` or: column name containing "sample" | Definition of sample |
| `ID` | - |
| any of: `.logTransformed`, `.identity`, `.CmbCor`, `.value` | value columns |
| any of: `Test Condition`, `Supergroup`, or other column(s) | Conditions |

## Limma

```
Limma_<PTK|STK>_<number>_<Group>.csv
```
- Any name. Must not contain underscores. Must be identical between the PTK and STK files. 
- Supergroup names and comparisons come from the file content, not the filename.

## UKA

```
UKA_<PTK|STK>_<number>_<anything>.csv
```
