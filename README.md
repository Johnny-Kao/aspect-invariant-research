# ASPECT invariant research

Independent research repository (not a fork). Upstream: https://github.com/geodynamics/aspect at `540a201c460759e160968bfc03323db8e7b76367`. Initial research repository contains test harnesses, not a complete upstream mirror. Harnesses retrieve pinned upstream files, including the Pyrolite table, and compile selected arithmetic using minimal stand-ins for ASPECT/deal.II data structures. ASPECT source is GPL-2.0-or-later; upstream attribution and license apply.

| Candidate | Invariant | Check |
|---|---|---|
| AS-04 | Depth and translation independence | Extraction reaction arithmetic, original vs candidate correction |
| AS-05 | Phase transition index consistency | Density interpolation, original vs candidate correction |
| AS-06 | Local water mass-fraction conservation | Zero solubility reaction, original vs candidate correction |
| AS-07 | Lookup domain and derivative consistency | Actual upstream table, baseline counterexamples only |

Run `python3 research/check_as04_as07.py` or specify a candidate (e.g. AS-04). CI success means the expected counterexamples and controls were observed, not that upstream passes the corrected invariant. No upstream submission or full-model test is implied. AS-05 overlaps open upstream PR #7331; AS-06 is zero-solubility only; AS-07 derivative policy and endpoint correction require additional design. Full native regression testing remains separate. 
