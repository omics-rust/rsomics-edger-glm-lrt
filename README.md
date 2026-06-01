# rsomics-edger-glm-lrt

Negative-binomial GLM fit and likelihood-ratio test of a coefficient or contrast
— a Rust port of edgeR's `glmFit()` + `glmLRT()` + `topTags(sort.by="none")`.

For each gene it fits the full design by NB IRLS (log link, per-sample offset
`log(lib.size × norm.factors)`), refits the reduced model (the tested coefficient
dropped, or the design rotated so a contrast becomes one coefficient), and
reports:

| column | meaning |
|---|---|
| `logFC` | log2 fold change of the tested coefficient/contrast (prior.count 0.125 shrinkage, edgeR's `predFC`) |
| `logCPM` | average log2-CPM (`aveLogCPM`, prior.count 2, dispersion 0.05) |
| `LR` | likelihood-ratio statistic = drop in NB deviance from the reduced model |
| `PValue` | chi-square upper tail of `LR` on the tested degrees of freedom |
| `FDR` | Benjamini-Hochberg adjusted p-value (with `--fdr`) |

## Usage

```
rsomics-edger-glm-lrt counts.tsv --design design.tsv [--coef N | --contrast c.tsv] \
    [--dispersion D | --dispersion-file f.tsv] [--norm-factors f.tsv] [--fdr] [-o de.tsv]
```

- `counts.tsv` — header `gene<TAB>sample…`, one integer-count row per gene.
- `design.tsv` — header of coefficient names, then one numeric row per sample.
- `--coef N` — 1-based design column to test (default: the last). `--contrast`
  takes a per-coefficient weight vector instead.
- `--dispersion` — common NB dispersion (default 0.05); `--dispersion-file`
  gives per-gene values (gene order).

```
rsomics-edger-glm-lrt counts.tsv --design design.tsv --coef 2 --dispersion 0.1 --fdr -o de.tsv
```

## Origin

This crate is an independent Rust reimplementation of edgeR's `glmFit`/`glmLRT`
based on:

- The published method: McCarthy DJ, Chen Y, Smyth GK, "Differential expression
  analysis of multifactor RNA-Seq experiments with respect to biological
  variation", *Nucleic Acids Research* 40(10):4288-4297, 2012.
  DOI: 10.1093/nar/gks042. Robinson MD, Smyth GK (2008) for the NB framework.
- The public edgeR R-level interface (function signatures and the documented
  `predFC` prior-count shrinkage / `aveLogCPM` definitions), observed
  black-box.
- Black-box behaviour testing against the edgeR binary (`glmFit` + `glmLRT` +
  `topTags`).

No source code from the GPL edgeR/limma C internals was read or used during
implementation; the IRLS NB-GLM fit and the LRT statistic are implemented from
the published method. The special functions (pchisq via the incomplete gamma)
are ports of the public-domain numerical algorithms R itself uses, not edgeR
code. Test fixtures are independently generated count matrices.

LR/PValue/FDR/logCPM match edgeR to ~1e-6; logFC matches to the slack of
edgeR's own `tol=1e-6` Levenberg stopping (ours converges to the augmented MLE).

License: MIT OR Apache-2.0.
Upstream credit: edgeR <https://bioconductor.org/packages/edgeR/> (GPL ≥ 2).
