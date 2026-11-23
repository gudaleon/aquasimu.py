# NOTICE — biowenp_py_AN_v1

This package is a Python port, with additions, of **BioWEN-P**, the phosphorus reactive-transport model for the
Elm Drive bioretention cell (Mississauga, Ontario) written in MATLAB by Bowen Zhou (Copyright 2020-2023 Bowen Zhou).

Please cite the original work:

* Zhou, B., Shafii, M., Parsons, C.T., Passeport, E., Rezanezhad, F., Lisogorsky, A., Van Cappellen, P. (2023).
  Modeling multi-year phosphorus dynamics in a bioretention cell: Phosphorus partitioning, accumulation, and export.
  *Science of the Total Environment* 876, 162749. https://doi.org/10.1016/j.scitotenv.2023.162749
* Zhou, B., Shafii, M., Parsons, C., Passeport, E., Rezanezhad, F., Lisogorsky, A., Van Cappellen, P. (2023).
  Modeling multi-year phosphorus dynamics in a bioretention cell: phosphorus partitioning, accumulation, and export
  [model scripts and modelled results]. Federated Research Data Repository. https://doi.org/10.20383/103.0640

Licence: the original `.m` files state the GNU Lesser General Public License v3 or (at your option) any later
version; the FRDR dataset README states CC BY 4.0 for the data. This port is a derivative of the LGPL code and is
therefore distributed under **LGPL-3.0-or-later** (see `COPYING.LESSER` and `COPYING`). It comes WITHOUT ANY
WARRANTY. Modified by Alex Neumann / EnviroManagement, November 2026: translated to Python; tridiagonal solver
instead of dense `A\B`; modular layout; mass-balance ledgers; augmentation modules (`flowthrough_AN_v1.py`,
`priors_AN_v1.py`) for saturated eggshell + woodchip flow-through reactors, which are NOT part of BioWEN-P.

The original MATLAB files are NOT redistributed here. They are read (unmodified) from the local FRDR snapshot only
for the Octave equivalence check.

Augmentation priors: Lee et al. (2022) Chemosphere 287:132267 (calcined eggshell); Addy et al. (2016) J. Environ.
Qual. 45:873 and Schipper et al. (2010) Ecol. Eng. 36:1532 (woodchip denitrification).
