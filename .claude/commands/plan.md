You are helping plan edits to this project before making any changes.

The user's request is: $ARGUMENTS

## Step 1 — Identify the target file(s)

Based on the request, identify which file(s) are most relevant:
- LaTeX text edits → `paper/PowerSim/PowerSim_Paper2.tex` or `paper/PowerSim/supplementary.tex`
- R code edits → files under `code/` or `examples/`
- Documentation → files under `doc/`

Read the relevant file(s) or sections. Do NOT make any edits yet.

## Step 2 — Route to planner based on file type

**If the target is a `.tex`, `.md`, or `.bib` file (prose/text):**

Run Codex in headless mode to generate the plan:

```bash
codex exec \
  --ephemeral \
  --sandbox read-only \
  -o /tmp/dcpower_plan.txt \
  "You are a scientific writing assistant for a biostatistics paper about circadian rhythm power analysis (DCPower). 
   
   The user wants to: $ARGUMENTS
   
   Read the relevant section of the file provided below and produce a structured edit plan with:
   1. WHAT to change (quote the exact current text)
   2. WHERE it appears (section name + approximate line range)  
   3. WHY the change is needed
   4. PROPOSED NEW TEXT (write the exact replacement)
   
   Be concise. Do not make any edits — only plan.
   
   File content:
   $(cat <RELEVANT_FILE>)"
```

Then read `/tmp/dcpower_plan.txt` and present the plan clearly.

**If the target is an `.R`, `.sh`, or other code file:**

Generate the plan yourself (do not call Codex):
1. Read the file and identify the relevant function/section
2. State exactly what needs to change and why
3. Show the proposed diff (old → new)
4. Note any downstream files that may be affected

## Step 3 — Present the plan

Format the output as:

---
**FILE:** `path/to/file.tex` (lines X–Y)  
**CHANGE:** one-sentence summary  
**CURRENT TEXT:**  
> exact quoted text

**PROPOSED TEXT:**  
> exact replacement

**REASON:** why this change is correct/needed
---

If multiple changes are needed, number them (Change 1, Change 2, ...).

## Step 4 — Wait for approval

After presenting the plan, ask:
> "Proceed with these changes? (yes / modify / cancel)"

Do NOT make any edits until the user confirms.
