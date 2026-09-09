# Reference output

Numeric logs from the runs the paper's tables were read from, kept so that a re-run can be compared line
by line without access to the authors' machine. `analysis/output/` itself is gitignored and regenerated.

| File | Script | Paper |
|---|---|---|
| `13_hierarchical_conformal.txt` | `analysis/13_hierarchical_conformal.R` | Table 9 (design 1, `cov_assess`/`width_assess` columns) and Table 10 (design 2) |
| `13_hierarchical_conformal_summary.csv` | same | the same numbers as a table |

The runs are seed-fixed (`set.seed(1)`, `seed = 1` in every `ranger()` call), so a re-run on the same
package versions reproduces them exactly; across `ranger` versions small differences in the third
decimal are possible.
