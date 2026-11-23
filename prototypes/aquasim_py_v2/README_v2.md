# aquasim_py v2: SOD_All port (SEDCOMP + formula evaluator + processes)
New in v2 (v1 is unchanged):
- aquasim_py/formula.py: CONSTVAR, FORMVAR (^, if/then/else/endif -> np.where, log10 ...), REALLISTVAR/VARLISTVAR LINEAR, PROGVAR (t,z,theta,calcnum,zone); evaluates whole grid columns at once.
- aquasim_py/sod_model.py: SEDCOMP InitCond/Delta port (highres=FALSE, Qinterf=0), AQCOMP CalcRates, banded FD Jacobian (make_jac).
  The DAE algebraic rows (Q, porosity, top/bottom BC) are solved explicitly, so this is an ODE in inner cells (248 x 30 = 7440 unknowns).
- run_sod_v2.py <aqu> <out.npz> [rtol]   (scipy BDF, rtol=atol=1e-6)
- compare_v2.py <out.npz> ref_profiles_v1.zip
- massbal_v2.py <aqu> <out.npz>
Result (box, 8 cores): 198 s (AQUASIM reference binary: 515 s). Max |diff|/max(profile) per variable is at most 1.8e-3 at all 8 times; the REF lists have 4 significant digits.
