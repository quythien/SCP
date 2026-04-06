You are helping draft or revise a figure caption for the DCPower paper.

The user's request is: $ARGUMENTS

## Step 1 — Identify the figure

Based on the request, identify which figure is being captioned:
- `bootstrap_summary.pdf` → Fig 1 (two-stage vs bootstrap)
- `multi_dataset_dr_power.pdf` → Fig 3 (DR power across datasets)
- `bm_tradeoff.pdf` → Fig 4 (B vs m tradeoff)
- `dr_power.pdf` → DR power stratified
- `dp_power.pdf` → DP power stratified
- `fourier_robustness.pdf` → Supplementary: cosinor assumption robustness

Read the current caption from `paper/PowerSim/PowerSim_Paper2.tex` (or `supplementary.tex`).

## Step 2 — Gather figure context

Read the relevant Results section text surrounding the figure to understand:
- What the figure shows
- Key findings being illustrated
- Any specific numbers mentioned in the text

## Step 3 — Route to Codex to draft the caption

```bash
codex exec \
  --ephemeral \
  --sandbox read-only \
  -o /tmp/dcpower_caption.txt \
  "You are a scientific writing assistant for a biostatistics paper on circadian rhythm power analysis.

   Draft a figure caption following these conventions:
   - Start with a bold one-sentence overview of what the figure shows
   - Then describe each panel (A, B, C, ...) in order
   - Include key quantitative findings where relevant
   - End with interpretation (what the reader should conclude)
   - Tone: concise, precise, passive voice where appropriate
   - Length: 3-6 sentences
   - No LaTeX formatting — plain text only
   
   Figure: $ARGUMENTS
   
   Current caption (if any):
   $(grep -A 10 'caption{' paper/PowerSim/PowerSim_Paper2.tex | head -20)
   
   Surrounding results text:
   $(grep -B 5 -A 20 'fig:' paper/PowerSim/PowerSim_Paper2.tex | head -40)"
```

Read `/tmp/dcpower_caption.txt` and display the drafted caption.

## Step 4 — Review

Show the drafted caption and ask:
> "Use this caption? (yes / modify / cancel)"

If modify, ask what to change and pass back to Codex for revision.

## Step 5 — Insert into paper

Once approved:
1. Locate the `\caption{...}` block for the figure in `PowerSim_Paper2.tex`
2. Replace it with the approved text (properly formatted in LaTeX)
3. Show the diff before saving
