Example upload pilot for the SCP Shiny app ("Upload my own pilot").

  example_expression.csv  2500 genes x 60 samples; genes in rows, samples in
                          columns (first column = gene IDs, header = sample IDs).
  example_tod.csv         one collection time of day (hours) per sample, one per
                          line, NO header (matches the upload reader).

These are SIMULATED to resemble a GTEx Adrenal pilot (generated via
SCP::generatePilotData); they are NOT real GTEx expression or time-of-death,
so they are freely shareable. Swap in your own CSVs of the same shape.
