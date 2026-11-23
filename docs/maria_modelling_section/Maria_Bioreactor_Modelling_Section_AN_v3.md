# Modelling: mass-balance tool for eggshell–woodchip bioreactors at wetland outflows

*Canada Water Agency Great Lakes Freshwater Ecosystem Initiative — restoring and protecting coastal areas.*

**Purpose.** This component will estimate how much nitrogen (N) and phosphorus (P) modular bioreactors can trap, and later recycle, at each installed bioreactor and at the wetland outflows they treat in the Severn Sound and Georgian Bay coastal area, in partnership with Georgian Bay Forever and Severn Sound. It supports eligible program activities: numerical modelling with bio-physical data collection; pre- and post-implementation assessment, monitoring, and reporting; and transferable science on wetland and nearshore water quality, including nutrients. The reactors use calcined eggshell and woodchips, building on pilots in the Bay of Quinte and the Lake Erie catchment near Windsor. A batch result of up to 80% phosphate removal with 3–5 total weight% eggshell may inform priors and design, but will not be treated as field performance.

**Mass-balance framework.** The analysis opens with an explicit water, N, and P mass balance at each installed bioreactor and at the wetland outflow it treats. Untreated bypass during high flow is part of the balance. A calibration will not be accepted until these balances close to a reported residual.

- **Water:** inflow = treated outflow + untreated bypass + change in storage + residual.
- **Nitrogen:** inflow load = outflow load + bypass load + denitrification to N₂/N₂O + plant or biomass uptake + immobilization in the woodchip matrix + change in storage + residual.
- **Phosphorus:** inflow load = outflow load + bypass load + sorption and precipitation on calcined eggshell (including hydroxyapatite formation) + plant or biomass uptake + change in storage + residual.

Each residual is a **closure check** and will be reported, not absorbed into a removal term. Intermediate checks and the full input tables will be published alongside the outputs. The model will be versioned, inputs archived, and results recomputed when data are updated. Comparisons with field monitoring are consistency checks on model components, not “validation,” unless the spatial and temporal scale of the observations matches the scale of the prediction.

**Process model.** The process representation follows AQUASIM, Eawag’s program for identifying and simulating aquatic systems (Reichert 1994, 1998). Each bioreactor will be represented as a mixed-reactor or advective-diffusive (longitudinal flow-through) compartment, with an optional biofilm compartment if depth gradients in the media matter. High-flow bypass will be a linked untreated pathway, so it remains in the mass balance. N and P removal will be user-defined first-order or Michaelis–Menten kinetics with temperature dependence, informed by the woodchip and eggshell literature cited below, but calibrated rather than fixed. Sensitivity analysis will be used to show which rate constants, residence-time terms, and temperature coefficients are constrained by the monitoring data.

The modelling framework is implemented as a version-controlled Python code base, aquasimu.py, that ports the AQUASIM framework and is maintained in a GitHub repository (https://github.com/gudaleon/aquasimu.py), ensuring reproducibility of model structure, inputs, and results.

**Calibration: PEST++ and DREAM.** Both PEST++ and DREAM will be used to calibrate the same process model against the monitoring data, and their estimates will be compared. PEST++ will provide gradient-based parameter estimation, identifiability analysis, and linear uncertainty analysis (Doherty 2018; White et al. 2020). DREAM (DiffeRential Evolution Adaptive Metropolis), including the DREAM(ZS) variant, will provide Bayesian Markov chain Monte Carlo sampling, posterior parameter distributions, and credible intervals on annual load reductions (Vrugt et al. 2009; Vrugt 2016). In both cases the Python process model is run as the inner simulator: the calibrator proposes parameters, the model returns simulated concentrations and loads, and mass-balance closure is checked before a result is accepted. Priors and starting values will be informed by the literature and by the Quinte and Windsor work. The main results will be estimated kilograms of N and P retained per year at each installed bioreactor and at the wetland outflow, with uncertainty reported from both methods. Where supported by the data, the same relationships will be used to scale estimates to the contributing subwatershed and the Severn Sound and Georgian Bay coastal context.

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
