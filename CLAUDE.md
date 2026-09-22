# CLAUDE.md

## Local R environment for dev/testing

Production is pinned to R 4.1.0 (`Dockerfile`'s `rocker/r-ver:4.1.0`, matching `renv.lock`).
On Windows dev machines, R 4.1.0 itself is typically not installed — only later versions
(e.g. 4.2/4.3/4.4), each with its own `renv/library/R-<version>/` folder. These folders are
**not equally valid substitutes**:

- Check which one's packages were actually restored from `renv.lock` (e.g. does it have
  `ggplot2` at the locked version?) — that one approximates production; treat it as the
  project's real local environment.
- A folder that only has a few packages installed ad hoc (e.g. `testthat`, which isn't even
  in `renv.lock`) is not renv-managed and is not a substitute for the real one, even if it
  happens to be more convenient (e.g. it's the one with the test runner installed).
- Missing dev-only tooling (like `testthat` and its dependencies) in the real environment
  should be installed there, rather than routing tests through whichever environment
  already has it installed.
- This never changes the code itself — plain R files/`.Rmd`s aren't version-specific. It
  only affects which local `Rscript.exe` you invoke to run tests/scripts during
  development. Production continues to run on R 4.1.0 via the unchanged Docker image.
- Known drift (as of 2026-09): R-4.3's renv library has `flextable` 0.9.2 installed, not
  the `renv.lock`-pinned 0.6.10. This breaks `render_ref_qc_table()`'s `footnote()` call
  with a "subscript out of bounds" error in `tests/testthat.R` — unrelated to whatever
  feature you're testing if that's the only failure you see; production isn't affected
  (its Docker image installs the pinned 0.6.10). Fix by reinstalling the pinned version if
  it starts blocking real verification work.

## Domain modeling - building glossary into CONTEXT.md

Use precise terms of the domain / project instead of vague ones. Use /mattpocock-domain-modeling skill to build the terms of the domain / project into the CONTEXT.md glossary. When you understand new terms from conversation, suggest a definition and ask the user whether you can update the CONTEXT.md. 

## PLAN FORMAT
When writing or updating a plan file (plan mode), structure it as:

1. Goal — 1-few sentences, the actual underlying question/objective.

2. Key decisions: a concise, "technical manager level" overview with bullet points.

3. Requirements:

- `R1: When <trigger>, the <system> shall <response>.`
- `R2: While <state>, the <system> shall <response>.`

4. Implementation:

- `S1: <edit(s) that satisfy R1/R2>.`
- `S2: USER checkpoint: review/commit chunk 1.`

5. Verification:

- `V1 (R1,R2): npm test`

6. Rejected approaches, IF relevant. List the working approaches that have been reverted. Only include approaches that have a substantial effect on methodology or design — not every minor UI or presentation iteration. Each entry: what it was, why it was rejected, in one or two lines. No code references needed here — keep it conceptual.

### REQUIREMENTS FORMAT (EARS)

Write compact, testable requirements about the system/component under change (not the agent). Name the system explicitly (e.g. `InkNewWizard`, `bep new`, `WizardState`).

- Ubiquitous (always true): `The <system> shall <response>.`
- State-driven: `While <precondition(s)>, the <system> shall <response>.`
- Event-driven: `When <trigger>, the <system> shall <response>.`
- Optional feature/scope: `Where <feature/scope applies>, the <system> shall <response>.`
- Unwanted behavior: `If <unwanted condition>, then the <system> shall <mitigation>.`
- Complex: `While <precondition(s)>, when <trigger>, the <system> shall <response>.`

Practical rules:

- Use requirement IDs (`R1`, `R2`, ...) so implementation and verification can reference them.
- Prefer observable behavior and invariants; avoid file/function names unless they are part of the external contract.
- Include: key outputs/deliverables (with exact path).

### IMPLEMENTATION PLAN FORMAT

Describe *how* you'll satisfy the requirements as concrete steps (agent actions), chunked into small git-committable units when appropriate.

- Size the steps to the change: use as few steps as needed for small fixes, and break larger changes into multiple git-committable chunks.
- Keep one concrete outcome per step (code change, test addition, verification, or user checkpoint).
- Include a USER checkpoint step for major or risky changes, consistent with the workflow above.

### VERIFICATION FORMAT

Include explicit checks that map back to the requirements.

- Each verification item should reference one or more requirement IDs (`R#`) and name the check (`npm test`, `npm run build`, or targeted manual validation).

### General rules:

- the language of the plan should use the terms from CONTEXT.md
- No filler text. Cut sentences that restate the obvious. Every sentence should carry information the reader doesn't already have.
- When revising an existing plan file, treat it as a **cumulative document** — update/extend it, don't replace the whole file with only the newest increment. The plan should always read as the full current picture, not a diff.
- If a description (e.g. of a rejected approach) gets too long, you should offer it to replace it into a different relevant documentation.
- Keep it concise and scannable: prefer short bullets and small tables over prose paragraphs. Use subtitles to make it skimmable.
- Write as if the reader has no memory of the conversation that produced it — spell out context instead of saying "as discussed" or assuming prior state is known.
- Don't restate implementation history as commentary in the Overview (e.g. "this column is carried over from an older version") — describe what's true and relevant now.
- Before a substantial rewrite/restructuring of an existing plan (not a small edit), ask whether the user wants a timestamped snapshot of the current version saved first (`yy-mm-dd-hh-mm_<filename>.md`, same folder) — don't snapshot automatically every time.