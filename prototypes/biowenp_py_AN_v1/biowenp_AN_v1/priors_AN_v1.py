# Part of biowenp_py_AN_v1 (LGPL-3.0-or-later). Augmentation for eggshell + woodchip reactors (not in BioWEN-P).
"""Literature priors for the augmented flow-through model. These are PRIORS for calibration (PEST++/DREAM),
not fixed values and not field performance.

Eggshell P sorption -- Lee et al. (2022) Chemosphere 287:132267, doi:10.1016/j.chemosphere.2021.132267
  (calcined eggshell, CES-800; values from the abstract/highlights):
  * Langmuir maximum adsorption capacity (batch): 108.2 mg P/g  = 108.2 g P/kg calcined eggshell
  * pH 3 -> 11: adsorbed amount 99.8 -> 62.3 mg P/g
  * dynamic column with real lake water: >90 % removal up to 158 h; 11.5 mg P/g adsorbed
  * eggshell calcined at <= 700 C: < 1 mg P/g (negligible) -> calcination state matters
  * mechanism: Ca elution and hydroxylapatite precipitation (so "sorption" includes precipitation)
  K_L and the pseudo-first-order rate constant are NOT available from the abstract; wide priors are used and
  must be replaced by Maria's own batch isotherm/kinetic tests (or the full text of Lee et al.).
  Zhou et al. (2023) Qmax = 683.5 mg/kg is for SANDY filter media and is deliberately NOT used here.

Woodchip nitrate removal -- Addy et al. (2016) J. Environ. Qual. 45:873, doi:10.2134/jeq2015.07.0399:
  bed removal rate R = 1.79 exp(0.0766 T) g N per m3 BED per day (Q10 = 2.15) -> R20 = 8.3 g N/m3/d,
  theta = exp(0.0766) = 1.080; ~2 g N/m3/d below 6 C.  Schipper et al. (2010) Ecol. Eng. 36:1532: 2-22 g N/m3/d.
  Converted to the pore-water basis of Maria's eq. (10)-(11) by dividing by the woodchip porosity.
"""
import numpy as np

LEE2022 = dict(qmax_batch_g_per_kg=108.2, qmax_column_g_per_kg=11.5, qmax_pH3=99.8, qmax_pH11=62.3,
               uncalcined_lt=1.0)
ADDY2016 = dict(a=1.79, b=0.0766, q10=2.15)


def eggshell_priors(eggshell_mass_fraction=None):
    """Priors (name: (median, low, high, transform, source)). q_max in g P per kg EGGSHELL. If a mass fraction of
    eggshell in a mixed bed is given, a derived per-kg-bed capacity is also returned."""
    pri = {
        "q_max": (np.sqrt(11.5 * 108.2), 11.5, 108.2, "log", "Lee2022: column 11.5 to batch 108.2 g/kg"),
        "K_L": (0.5, 0.01, 10.0, "log", "no value in Lee2022 abstract; wide, replace with batch isotherm"),
        "k_ad": (2.0, 0.05, 50.0, "log", "1/d; no value in Lee2022 abstract; wide"),
        "k_mi": (0.0, 0.0, 0.5, "lin", "g P/kg/d zero-order sorbed->mineral (Zhou-type); 0 by default"),
        "theta_P": (1.02, 1.0, 1.06, "lin", "weak temperature effect on sorption; assumption"),
    }
    if eggshell_mass_fraction is not None:
        f = eggshell_mass_fraction
        pri["q_max_per_kg_bed"] = (pri["q_max"][0] * f, 11.5 * f, 108.2 * f, "log", f"derived, eggshell fraction {f}")
    return pri


def woodchip_priors(porosity=0.65):
    r20_bed = ADDY2016["a"] * np.exp(ADDY2016["b"] * 20)
    return {
        "vmax_N20": (r20_bed / porosity, 2.0 / porosity, 22.0 / porosity, "log",
                     "g N/m3 pore/d; Addy2016 R20=8.3 bed basis; Schipper2010 2-22"),
        "k_N20": (1.0, 0.1, 5.0, "log", "1/d first-order alternative; wide"),
        "K_N": (2.0, 0.2, 10.0, "log", "g N/m3 Michaelis-Menten half-saturation; wide"),
        "theta_N": (float(np.exp(ADDY2016["b"])), 1.04, 1.12, "lin", "Addy2016 Q10=2.15 -> theta=1.080"),
    }


def prior_table():
    rows = []
    for grp, d in (("eggshell", eggshell_priors(0.04)), ("woodchip", woodchip_priors())):
        for k, (med, lo, hi, tr, src) in d.items():
            rows.append(dict(group=grp, parameter=k, median=med, low=lo, high=hi, transform=tr, source=src))
    return rows
