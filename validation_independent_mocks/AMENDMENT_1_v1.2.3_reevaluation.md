# Amendment 1: post hoc re-evaluation of sets A, D and E with MetagenApp v1.2.3

Date: 2026-11-12. Committed before this re-evaluation is run. It does not replace or
modify the pre-specified plan (`PREREGISTRATION.md`, commit e353b80) or its results.

## 1. Why this amendment exists (disclosure)

The pre-specified analysis of MetagenApp v1.1 was run and its results were seen:
v1.1 discarded almost all reads of sets A, D and E because primer trimming was anchored
to the first base of the merged read, and these libraries carry spacers or staggered
(phased) primers before the forward primer. That result stands as the primary result
and is reported as obtained.

Sets A, D and E were then used to diagnose the failure and to test a fix (unanchored
primer trimming). An intermediate analysis (v1.1 with only `trim_primers.py` replaced)
was run on these sets and on the original mock communities; its exact code, parameters
and outputs are archived separately. MetagenApp v1.2.3 includes that fix together with
other changes made since v1.1 (per-sample de novo chimera removal by default, and
genus-level reporting of unresolved assignments).

Therefore this re-evaluation is **not blind**: the libraries have already been used to
develop the fix. It evaluates the currently released version as a whole; differences
from v1.1 cannot be attributed to the primer fix alone.

## 2. What is run

- Sets A, D and E only (the same 10 libraries and inputs as the pre-specified analysis).
  The main benchmark and the original mock communities are not re-run; the manuscript
  keeps citing the version that produced them (v1.1 / v1.1.1).
- **MetagenApp v1.2.3** as released (tag v1.2.3, commit b97a9c6;
  doi:10.5281/zenodo.24170815), with the same reference files
  (doi:10.5281/zenodo.23326733), default settings of v1.2.3 (including
  `--chimera-mode per-sample`), and the same command-line parameters as the
  pre-specified analysis:
  `--mode ref --marker 16S --classifier wang --primer-f CCTACGGGNGGCWGCAG
  --primer-r GACTACHVGGGTATCTAATCC --min-length 400 --max-length 440 --max-ambigs 0
  --max-poly 8 -t 16`
- mothur and QIIME 2 are not re-run; their pre-specified results are used for comparison.

## 3. Metrics and reporting

Same metrics, same formula and same expected composition as Section 6 of the
pre-specified plan, per set and pooled. Results are reported whatever their direction,
alongside (never in place of) the pre-specified v1.1 results, labelled as a post hoc,
non-blind re-evaluation of v1.2.3. No parameter is changed after seeing these results.
