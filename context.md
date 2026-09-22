# Domain context: autoreport2023

Terminology used throughout this codebase and its docs (`ref/HowToUse.md`,
`R/00_GeneralFunctions.R`, `R/01_BasicProcessing.R`, `R/02_PhosphositeAnalysis.R`). Get
these wrong and QC/Limma/TT/MTvC/UKA logic reads as more confusing than it is.

## Condition ≠ Comparison

- **Condition**: a single experimental group/treatment/state. Examples: `T1`, `T2`,
  `Control`, `Parental`. In QC files this is (part of) the `Test Condition` column value;
  it's what a sample/replicate *is*.
- **Comparison**: a pairing of two conditions being tested against each other. Written as
  `X vs Y` (or `XvsY`), e.g. `T1 vs Control`, `T vs C`. This is what stats apps (Limma, TT,
  MTvC, UKA) report *on* - it's not itself a condition, and grouping data "by comparison"
  is a different operation than grouping "by condition".

`"Test vs Control"` is a comparison of the `Test` and `Control` conditions - not itself a
condition. Don't use `X vs Y`-shaped strings as condition examples.

## Supergroup

A grouping/stratification factor that comparisons are nested within - e.g. distinguishes
different cell lines, timepoints, or experimental batches when a study contains more than
one. Comparisons are often reported prefixed by their Supergroup, e.g. `Sgroup1 - T vs C`
(UKA's `Sgroup_contrast` column - Supergroup + contrast concatenated - exists specifically
to carry both at once). A study with only one Supergroup may omit it or use a placeholder
name; it's still a distinct concept from both Condition and Comparison.

## Where "condition" means something slightly wider (QC variability calc)

For the Data Variability Indicator's per-condition SD calculation
(`classify_qc_columns()`/`compute_condition_sd()` in `R/01_BasicProcessing.R`), "condition"
means *whatever combination of columns groups replicate samples into distinct experimental
units* - normally `Supergroup` + `Test Condition`, but in a technical-replicate (TR) QC
file that also carries a biological-replicate factor (e.g. `Biological REP`, values like
`BR1`/`BR2`/`CRL1`), that factor is itself folded into the condition grouping (each
biological replicate becomes its own distinct condition for this specific calculation) -
see `ref/median_sd_variability_notes.md`. This is a QC-calculation-specific usage, not a
redefinition of "condition" elsewhere in the app.

## Assay type

`PTK` (protein tyrosine kinase) / `STK` (serine-threonine kinase) - the two kinase activity
assay types a PamGene study can include, independently of each other. Most tables/flags are
computed per assay type.

## Replicate type (TR / BR)

- **TR** (technical replicate): repeated measurement of the *same* biological sample (e.g.
  multiple PamChip arrays/spots). Sample axis is `Barcode`+`Row`.
- **BR** (biological replicate): distinct biological samples (donors, independent
  experiments). When a BR file is produced by averaging up TR data, the per-array sample
  columns are gone - replaced by a `Biol_Rep`/"sample name"-ish column identifying which
  biological replicate each averaged row came from.

See `ref/HowToUse.md` for the exact QC filename/column conventions, and
`C:\Users\dschuller\.claude\plans\new-feature-the-qc-swirling-thompson.md` for how TR vs BR
drives QC flagging and the Data Variability Indicator.
