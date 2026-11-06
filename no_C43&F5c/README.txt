Supplementary files for “Exact Turán densities in triple systems”

Contents
--------
- no_K4F5c.ipynb: SageMath notebook used for the flag-algebra computations.
- Prop3-3_verify.ipynb: supplementary verification notebook for Proposition 3.3.
- prop_3_1_k4f5c.pickle, prop_3_2_k4f5c.pickle, prop_3_3_k4f5c.pickle:
  flag-algebra certificates for Propositions 3.1, 3.2, and 3.3.
- validator.py: independently reconstructs the objective and positivity
  assumptions for the certificates.
- verifier.py: checks the certificate identities, positive semidefiniteness,
  nonnegative slacks, and nonnegativity of every assumption multiplier.
- verify_exact_3_1_3_2.py: an additional exact-rational checker for the
  certificates of Propositions 3.1 and 3.2 (Python standard library only).
- repair_prop_3_2.py: exact reconstruction of the six-vertex objective and
  seven positivity rows, required by verify_exact_3_1_3_2.py.

The Proposition 3.2 certificate uses exact rational coefficients for the
minimum-degree threshold 4641/10000. The statement and numerical constants
of Proposition 3.2 are unchanged.

Suggested commands (run from this directory)
---------------------------------------------
python validator.py prop_3_1_k4f5c.pickle
python validator.py prop_3_2_k4f5c.pickle
python verifier.py prop_3_2_k4f5c.pickle --range 122 122
python verify_exact_3_1_3_2.py prop_3_1_k4f5c.pickle
python verify_exact_3_1_3_2.py prop_3_2_k4f5c.pickle

The commands using --range check only the selected identities, not all
constraints. To verify every constraint for Propositions 3.1 or 3.2, run
verify_exact_3_1_3_2.py without --range.

The displayed outputs saved in Prop3-3_verify.ipynb are historical;
they should not be confused with a fresh run of the current scripts.
For Proposition 3.3, the full verifier is computationally expensive.

The three .pickle certificate files are bundled in No_K4F5c.zip.
Extract the archive before running the verification scripts.
