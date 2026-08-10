# Building the cssig example data

```bash
Rscript data-raw/run-all.R
```

Run from the package root. Roughly 20 minutes on a modern laptop, dominated by
the forward Wright-Fisher stage. Stage 1 caches one `.rds` per chromosome in
`cache/`, so an interrupted run resumes; delete `cache/` to force a rebuild.

Needs `scrm` and `rehh`, both CRAN packages and both `Suggests` only — they are
used to *generate* the data and are never loaded by the installed package.

## Files

| File | What it does |
|---|---|
| `00-parameters.R` | Every constant. Nothing downstream hard-codes a value. |
| `01-simulate-genome.R` | Stage A: neutral coalescent background with `scrm`. Stage B: forward Wright-Fisher breed formation with sweeps. |
| `02-constituent-stats.R` | Chip-style SNP ascertainment, then F_ST / ΔDAF / ΔSAF via the package's own functions and XP-EHH via `rehh`. |
| `03-package-data.R` | Rounds, compresses and writes `data/`. |
| `run-all.R` | Driver; writes `provenance.txt`. |

## Design notes

**Why simulate genotypes rather than statistics.** F_ST, XP-EHH and ΔDAF are
correlated because they read the same underlying genealogy, and CSS exists to
exploit that. Simulating the statistics directly would require inventing that
correlation, making any demonstration built on them circular.

**Uniform treatment of all breeds.** Every breed runs forward for the same
number of generations at the same effective size, whether or not it carries a
sweep. If only swept breeds ran forward, they would accumulate extra drift and
F_ST would be raised across the whole region rather than at the sweep — a
dataset that looks convincing and is worthless.

**Conditioning on non-loss.** A beneficial mutation entering at a single copy is
lost with probability roughly 1 − s. Runs that lose it are retried; the retry
count per chromosome is recorded in `provenance.txt`, because conditioning is a
modelling choice and should be visible.

**Chunked coalescent.** Chromosomes are simulated in 10 Mb chunks to bound
memory. Chunks are independent, so LD does not carry across a chunk boundary.
With cattle LD decaying within 1–2 Mb this affects only windows straddling a
boundary, roughly 5% of them.

**Causal variants are excluded** from the final panel, as they would be on a
real SNP chip. The signal comes from hitchhiking neighbours, which is what the
source papers detect.

**Selection coefficients are large** (s = 0.10–0.50) because the window is only
40 generations. Under additive fitness an allele needs
`t = (2/s) * log[(p1/(1-p1)) / (p0/(1-p0))]` generations to go from p0 to p1, so
a complete sweep from one copy in 40 generations needs s near 0.5. These are
artificial-selection values, not natural-selection ones.

**Pre-thinning does not bias the frequency spectrum.** Sites carried into the
forward stage are chosen at random within spatial bins, not by MAF. Chip-like
frequency bias is introduced deliberately and separately, by the
discovery-panel ascertainment in `02-constituent-stats.R`.
