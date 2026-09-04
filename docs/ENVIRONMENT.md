# Environment and data sources

Recorded 2026-09-04 so that every step can be reproduced.

## Compute
- 4 CPUs, 15 GB RAM, Linux container. Python 3.11.
- Installed: sympy, numpy, scipy, networkx, z3-solver, python-sat, ortools, gmpy2, mpmath, numba.
- No Lean, Sage, GAP, or PARI. Workflow agent concurrency is capped at 2 (CPUs - 2).

## Network (egress proxy)
- **Works:** WebSearch (search snippets frequently quote erdosproblems.com, arXiv, OEIS pages), `github.com`, `raw.githubusercontent.com`, `git clone` of public GitHub repos, PyPI.
- **Blocked:** erdosproblems.com, en.wikipedia.org, oeis.org, arxiv.org, export.arxiv.org, mathoverflow.net, quantamagazine.org, terrytao.wordpress.com, tagteam.harvard.edu, api.semanticscholar.org, math.columbia.edu, emergentmind.com, *.github.io.

## Local copies used for problem statements and status
Cloned into the session scratchpad (not committed here; all public):
- `teorth/erdosproblems` — `data/problems.yaml`, ground truth for the status of all 1217 Erdős problems.
  Status counts on 2026-09-04: open 592, proved 334, disproved 139, solved 99, decidable 9, falsifiable 25, verifiable 7, plus independence-related categories.
- `google-deepmind/formal-conjectures` — Lean statements with docstrings (624 Erdős problems + OEIS, Green, MathOverflow, Wikipedia, and other collections).
- `mehmetmars7/Erdosproblems-llm-hunter` — archive of prior frontier-LLM attempts (TeX) on Erdős problems, useful for seeing what has been tried.

## Status vocabulary (from the erdosproblems database CONTRIBUTING.md)
- **decidable**: open, but reduced to a finite computation nobody has completed.
- **falsifiable**: open; a finite counterexample would disprove it.
- **verifiable**: open; a finite example would prove it.
