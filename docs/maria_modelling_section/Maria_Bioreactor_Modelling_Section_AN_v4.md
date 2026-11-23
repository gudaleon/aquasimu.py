# Modelling: mass-balance tool for eggshell–woodchip bioreactors at wetland outflows

*Canada Water Agency Great Lakes Freshwater Ecosystem Initiative — restoring and protecting coastal areas.*

**Purpose.** This component will estimate how much nitrogen (N) and phosphorus (P) modular bioreactors can trap, and later recycle, at each installed bioreactor and at the wetland outflows they treat in the Severn Sound and Georgian Bay coastal area, in partnership with Georgian Bay Forever and Severn Sound. It supports eligible program activities: numerical modelling with bio-physical data collection; pre- and post-implementation assessment, monitoring, and reporting; and transferable science on wetland and nearshore water quality, including nutrients. The reactors use calcined eggshell and woodchips, building on pilots in the Bay of Quinte and the Lake Erie catchment near Windsor. A batch result of up to 80% phosphate removal with 3–5 total weight% eggshell may inform priors and design, but will not be treated as field performance.

**Mass-balance framework.** The analysis opens with an explicit water, N, and P mass balance at each installed bioreactor and at the wetland outflow it treats. Untreated bypass during high flow is part of the balance. A calibration will not be accepted until these balances close to a reported residual, within a specified relative tolerance of $10^{-6}$. Pore volume is treated as constant, so the water storage change is zero. Over a simulation period the closures are

$$
V_{\mathrm{in}} = V_{\mathrm{t}} + V_{\mathrm{b}} + \varepsilon_Q \qquad (1)
$$

where $V_{\mathrm{in}}$, $V_{\mathrm{t}}$ and $V_{\mathrm{b}}$ are the inflow, treated-outflow and bypass volumes (m³) and $\varepsilon_Q$ is the water residual (m³);

$$
L_{N,\mathrm{in}} = L_{N,\mathrm{t}} + L_{N,\mathrm{b}} + L_{N,\mathrm{rem}} + \Delta S_N + \varepsilon_N \qquad (2)
$$

where $L_{N,\mathrm{in}}$, $L_{N,\mathrm{t}}$, $L_{N,\mathrm{b}}$ and $L_{N,\mathrm{rem}}$ are the inflow, treated-outflow, bypass and removed nitrate-N loads (g), $\Delta S_N$ is the change in dissolved N storage (g) and $\varepsilon_N$ is the N residual (g); and

$$
L_{P,\mathrm{in}} = L_{P,\mathrm{t}} + L_{P,\mathrm{b}} + L_{P,\mathrm{rem}} + \Delta S_P + \varepsilon_P \qquad (3)
$$

where the P terms (g) are defined analogously, and $L_{P,\mathrm{rem}}$ includes phosphorus transferred to the eggshell solid. Each residual is a **closure check** and will be reported, not absorbed into a removal term. Intermediate checks and the full input tables will be published alongside the outputs. The model will be versioned, inputs archived, and results recomputed when data are updated. Comparisons with field monitoring are consistency checks on model components, not “validation,” unless the spatial and temporal scale of the observations matches the scale of the prediction.

**Process model.** The process representation follows AQUASIM, Eawag’s program for identifying and simulating aquatic systems (Reichert 1994, 1998). Each bioreactor is a mixed-reactor or advective-diffusive (longitudinal flow-through) compartment, with an optional biofilm compartment if depth gradients in the media matter. Concentrations are expressed in g/m³ (equivalent to mg/L), flows in m³/d, and time in days.

Inflow above the design flow bypasses the unit and is remixed with the treated outlet:

$$
Q_{\mathrm{t}} = \min(Q, Q_{\mathrm{design}}),\qquad Q_{\mathrm{b}} = Q - Q_{\mathrm{t}} \qquad (4)
$$

where $Q$ is the inflow, $Q_{\mathrm{t}}$ the treated flow and $Q_{\mathrm{b}}$ the bypass (all m³/d). The flow-weighted outlet concentration is

$$
C_{\mathrm{out}} = (Q_{\mathrm{t}} C_{\mathrm{t}} + Q_{\mathrm{b}} C_{\mathrm{in}}) / Q \qquad (5)
$$

where $C_{\mathrm{t}}$ is the treated-outlet concentration and $C_{\mathrm{in}}$ the inflow concentration (g/m³).

For a completely mixed reactor (one cell, pore volume $V$, m³) the solute balance is

$$
V\,\frac{\mathrm{d}C}{\mathrm{d}t} = Q_{\mathrm{t}}(C_{\mathrm{in}} - C) - V r \qquad (6)
$$

where $C$ is the pore-water concentration (g/m³) and $r$ is the volumetric removal rate (g/m³/d). A woodchip or eggshell bed with longitudinal flow is represented as $n$ mixed cells of volume $V_c = V/n$ (tanks in series with optional physical dispersion), which is a finite-volume form of the advection–dispersion equation

$$
V_c\,\frac{\mathrm{d}C_i}{\mathrm{d}t} = Q_{\mathrm{t}}(C_{i-1} - C_i) + E(C_{i+1} - 2C_i + C_{i-1}) - V_c r_i \qquad (7)
$$

$$
E = (\alpha_L / L)\, n\, Q_{\mathrm{t}} \qquad (8)
$$

where $C_i$ is the concentration in cell $i$ (g/m³), $E$ is the inter-cell dispersive exchange (m³/d), $\alpha_L$ is the longitudinal dispersivity (m) and $L$ is the bed length (m). The inlet cell uses a Danckwerts (total-flux) condition with $Q_{\mathrm{t}} C_{\mathrm{in}}$. When $n = 1$ and $E = 0$, equation (7) reduces to equation (6).

Nitrate removal and phosphorus sorption are temperature-corrected with an Arrhenius $\theta$-factor at a reference temperature $T_{\mathrm{ref}} = 20\ ^\circ\mathrm{C}$:

$$
f(T) = \theta^{T - T_{\mathrm{ref}}} \qquad (9)
$$

where $T$ is water temperature (°C) and $\theta$ is dimensionless. First-order and Michaelis–Menten nitrate-N removal rates (g N/m³/d) are

$$
r_N = k_{N,20}\, f_N(T)\, C_N \qquad (10)
$$

$$
r_N = v_{\max,N}\, f_N(T)\, \frac{C_N}{K_N + C_N} \qquad (11)
$$

where $k_{N,20}$ is the first-order rate at $T_{\mathrm{ref}}$ (1/d), $v_{\max,N}$ is the maximum Michaelis–Menten rate (g N/m³/d), $K_N$ is the half-saturation constant (g/m³) and $C_N$ is the nitrate-N concentration (g/m³). These rates are informed by the woodchip literature cited below, but are calibrated rather than fixed.

Phosphorus retention on calcined eggshell is represented by kinetic Langmuir sorption with a finite solid capacity that depletes as the media loads (Lee et al. 2022):

$$
\frac{\mathrm{d}q}{\mathrm{d}t} = f_P(T)\bigl[k_{\mathrm{ads}} C_P (q_{\max} - q) - k_{\mathrm{des}} q\bigr] \qquad (12)
$$

$$
r_P = (\rho_b / \theta_w)\,\frac{\mathrm{d}q}{\mathrm{d}t} \qquad (13)
$$

where $q$ is sorbed P (g P/kg media), $q_{\max}$ is the sorption capacity (g P/kg), $C_P$ is the dissolved phosphate-P concentration (g/m³), $k_{\mathrm{ads}}$ is the adsorption rate (m³/g/d), $k_{\mathrm{des}}$ is the desorption rate (1/d), $\rho_b / \theta_w$ is the media mass per unit pore volume (kg/m³) and $r_P$ is the corresponding aqueous sink (g P/m³/d). An equivalent relaxation form, $\mathrm{d}q/\mathrm{d}t = f_P(T)\, k_{\mathrm{rel}}\bigl[q_{\max} C_P/(K_P + C_P) - q\bigr]$, is available when a Langmuir isotherm is approached at a finite rate $k_{\mathrm{rel}}$ (1/d). Sensitivity analysis will be used to show which rate constants, residence-time terms, and temperature coefficients are constrained by the monitoring data.

The modelling framework is implemented as a version-controlled Python code base, aquasimu.py, that ports the AQUASIM framework and is maintained in a GitHub repository (https://github.com/gudaleon/aquasimu.py), ensuring reproducibility of model structure, inputs, and results.

**Calibration: PEST++ and DREAM.** Both PEST++ and DREAM will be used to calibrate the same process model against the monitoring data, and their estimates will be compared. PEST++ will provide gradient-based parameter estimation, identifiability analysis, and linear uncertainty analysis (Doherty 2018; White et al. 2020). DREAM (DiffeRential Evolution Adaptive Metropolis), including the DREAM(ZS) variant, will provide Bayesian Markov chain Monte Carlo sampling, posterior parameter distributions, and credible intervals on annual load reductions (Vrugt et al. 2009; Vrugt 2016). In both cases the Python process model is run as the inner simulator: the calibrator proposes parameters, the model returns simulated concentrations and loads, and mass-balance closure is checked before a result is accepted. Priors and starting values will be informed by the literature and by the Quinte and Windsor work.

The PEST++ objective is the weighted sum of squares

$$
\Phi(\boldsymbol{\theta}) = \sum_i w_i^2 \bigl(y_i - \hat{y}_i(\boldsymbol{\theta})\bigr)^2 \qquad (14)
$$

where $y_i$ is an observation, $\hat{y}_i$ is the corresponding model prediction, $\boldsymbol{\theta}$ is the parameter vector and $w_i = 1/\sigma_i$ is the weight based on the observation standard deviation $\sigma_i$ (same units as $y_i$). The corresponding Gaussian log-likelihood used in DREAM (up to an additive constant when $\sigma_i$ is fixed) is

$$
\ell(\boldsymbol{\theta}) = -\frac{1}{2}\sum_i \left(\frac{y_i - \hat{y}_i(\boldsymbol{\theta})}{\sigma_i}\right)^2 \qquad (15)
$$

A heteroscedastic form is also used, in which the error standard deviation depends on the simulated value,

$$
\sigma_i = s_0 + s_1\,\hat{y}_i,\qquad
\ell(\boldsymbol{\theta}, s_0, s_1) = -\sum_i \log\sigma_i - \frac{1}{2}\sum_i \left(\frac{y_i - \hat{y}_i}{\sigma_i}\right)^2 \qquad (16)
$$

where $s_0$ has the same units as $y_i$ and $s_1$ is dimensionless. Separate $(s_0, s_1)$ pairs may be estimated for N and for P.

Annual retained loads are obtained by scaling the removed mass over a simulation of duration $\tau$ days:

$$
M_X = \frac{365}{\tau}\cdot\frac{1}{1000}\,L_{X,\mathrm{rem}} \qquad (17)
$$

where $X$ denotes N or P, $L_{X,\mathrm{rem}}$ is the removed load (g) and $M_X$ is the retained mass (kg/yr). The main results will be $M_N$ and $M_P$ at each installed bioreactor and at the wetland outflow, with uncertainty reported from both PEST++ and DREAM. Where supported by the data, the same relationships will be used to scale estimates to the contributing subwatershed and the Severn Sound and Georgian Bay coastal context. Cost-effectiveness is

$$
c_X = E / M_X \qquad (18)
$$

where $E$ is project expenditure over the project period and $c_X$ is the cost per kilogram of N or of P removed.

**Deliverables.** (i) Kilograms of N and P retained per year at each installed bioreactor and wetland outflow, with uncertainty from PEST++ and from DREAM and a reported closure residual. (ii) Cost-effectiveness, expressed as cost per kilogram of N and of P removed, calculated from project expenditures over the project period. (iii) Recycled-nutrient potential of spent eggshell media, with the same uncertainty treatment. (iv) A transferable, version-controlled modelling tool and an archive of inputs, intermediate checks, and outputs.

**Monitoring.** Monitoring will be designed to close the water, N, and P balances at each installed bioreactor. It will include flow at the reactor inlet, treated outlet, and any bypass; water temperature and dissolved oxygen; and nitrate (NO₃⁻), ammonium (NH₄⁺), total phosphorus (TP), and soluble reactive phosphorus (SRP). Sampling frequency will be set to capture event flows and the seasonal temperature range. Organic N and N₂O will be included where needed to close the organic and gaseous terms of the N balance. Spent-media sampling at the end of operating seasons will estimate P accumulated on eggshell and remaining woodchip N-removal capacity.

**Timeline.** Work will proceed over the project period, for up to four years as provided by the program. Year 1: compile existing pilot and batch data; implement the process model in the version-controlled code base; design the monitoring program with project partners; begin baseline outflow measurements. Year 2: support installation; intensive paired inflow–outflow–bypass monitoring through the first operating season; first joint PEST++ and DREAM calibration and closure report. Years 3 and 4: additional seasons, including high-flow events; updated estimates; subwatershed and coastal scaling; cost-effectiveness and recycled-nutrient estimates; release of the transferable tool and archive.

**References**

Addy, K., A.J. Gold, L.E. Christianson, M.B. David, L.A. Schipper, and N.A. Ratigan. 2016. Denitrifying bioreactors for nitrate removal: A meta-analysis. *Journal of Environmental Quality* 45: 873–881. https://doi.org/10.2134/jeq2015.07.0399

Doherty, J. 2018. *PEST: Model-Independent Parameter Estimation — User Manual*, 7th ed. Watermark Numerical Computing, Brisbane, Australia.

Lee, J.-I., J.-M. Kim, S.-C. Yoo, E.H. Jho, C.-G. Lee, and S.-J. Park. 2022. Restoring phosphorus from water to soil: Using calcined eggshells for P adsorption and subsequent application of the adsorbent as a P fertilizer. *Chemosphere* 287: 132267. https://doi.org/10.1016/j.chemosphere.2021.132267

Reichert, P. 1994. AQUASIM – A tool for simulation and data analysis of aquatic systems. *Water Science and Technology* 30(2): 21–30. https://doi.org/10.2166/wst.1994.0025

Reichert, P. 1998. *AQUASIM 2.0 – User Manual. Computer Program for the Identification and Simulation of Aquatic Systems.* Swiss Federal Institute for Environmental Science and Technology (EAWAG), Dübendorf, Switzerland. ISBN 3-906484-16-5.

Schipper, L.A., W.D. Robertson, A.J. Gold, D.B. Jaynes, and S.C. Cameron. 2010. Denitrifying bioreactors—An approach for reducing nitrate loads to receiving waters. *Ecological Engineering* 36: 1532–1543. https://doi.org/10.1016/j.ecoleng.2010.04.008

Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, B.A. Robinson, J.M. Hyman, and D. Higdon. 2009. Accelerating Markov chain Monte Carlo simulation by differential evolution with self-adaptive randomized subspace sampling. *International Journal of Nonlinear Sciences and Numerical Simulation* 10(3): 273–290. https://doi.org/10.1515/IJNSNS.2009.10.3.273

Vrugt, J.A. 2016. Markov chain Monte Carlo simulation using the DREAM software package: Theory, concepts, and MATLAB implementation. *Environmental Modelling & Software* 75: 273–316. https://doi.org/10.1016/j.envsoft.2015.08.013

White, J.T., R.J. Hunt, M.N. Fienen, and J.E. Doherty. 2020. Approaches to highly parameterized inversion: PEST++ Version 5, a software suite for parameter estimation, uncertainty analysis, management optimization and sensitivity analysis. *U.S. Geological Survey Techniques and Methods* 7-C26, 52 p. https://doi.org/10.3133/tm7C26
