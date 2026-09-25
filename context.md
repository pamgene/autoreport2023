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

## QC set

The QC files sharing one number token in their filename (`QC_<PTK|STK>_<number>_<TR|BR>...`,
e.g. all `QC_*_01_*` files vs all `QC_*_02_*` files). Numbered exports exist only when a
study was exported under more than one normalization approach (e.g. 01 = Log, 02 = VSN - see
`ref/HowToUse.md`), so one QC set = one normalization approach, covering PTK and STK and TR
and BR alike. PTK vs STK never makes two QC sets on its own. Unnumbered files (e.g.
`QC_PTK.csv`) and BioNavigator QC files form a single QC set.

## Executive Summary deck

The editable PowerPoint output (`01_REPORTS/03_ExecutiveSummary_PamDx_<yymmdd>.pptx`), built
with `officer` from the PamDx slide template next to the Main Report and Supplement. Holds
only summary content (QC tables, Data Variability Indicator, phosphosite table); dot plots,
Coral trees and conclusions are pasted/written by hand into space the template leaves empty.

## Marker shape

A shape in the slide template, named in PowerPoint's Selection Pane and showing `<...>` text,
that the code replaces with generated content (table, image, date). Its box defines the
content's position and maximum size - resizing the marker in the template changes the output
without code changes. Shapes that aren't marker shapes (static labels, "Replace this text.")
are never touched.

## Continuation slide

An extra "Results (continued)" slide holding the rows of a phosphosite table that don't fit
the `table_psite_analysis` marker shape; contains only the title and the table chunk (header
rows repeated).
