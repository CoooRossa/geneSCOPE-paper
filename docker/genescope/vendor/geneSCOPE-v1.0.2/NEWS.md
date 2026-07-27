# geneSCOPE 1.0.2

- Lee's L uses the canonical S2 normalization throughout the R and native-code paths.
- The default permutation mode is the all-grid joint row shuffle used in version 1.0.0. Block-constrained joint shuffling remains available with `use_blocks = TRUE`.
- `computeL()` and `getTopLvsR()` now share the same all-grid permutation default.
- `getTopLvsR()` excludes non-finite pair statistics before defining the eligible pair universe and records the selected and eligible family sizes in its provenance.
- L-vs-r results expose `Delta` and `delta_fdr`; the existing `fdr` column is retained as a compatible alias.
