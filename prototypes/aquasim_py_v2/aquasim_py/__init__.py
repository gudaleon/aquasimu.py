"""aquasim_py v1 - Python port skeleton of the AQUASIM 2.1g engine (aqcore).

Module map (C++ aqcore -> Python):
  fileio/readwords/all       -> aqufile.py   (.aqu brace-format reader/writer)
  aqvar/formvar/othervar/varsys/symbol -> variables.py (state, constant, formula, real-list, probe vars; formulas via sympy/numexpr)
  aqproc/procsys             -> processes.py (dynamic + equilibrium processes, stoichiometry)
  aqcomp + mixcomp/advcomp/filmcomp/soilcomp/rivcomp/lakecomp/sedcomp -> compartments/<name>.py
  aqlink/linksys             -> links.py
  aquasys/calcsys/statesys/compsys -> system.py (assemble global residual G(t,y,y')=0)
  dassl/das/linpack/eispack/numpar -> solver.py (DAE: scikits.odes IDA / assimulo IDA; ODE fallback scipy BDF)
  fit/fitsys/stat/zero/gamma -> fitting.py, sensitivity.py (scipy.optimize / lmfit)
  plotsys/drawpts/plotpts    -> plotting.py (matplotlib)
"""
__version__ = "0.1.0-v1"
