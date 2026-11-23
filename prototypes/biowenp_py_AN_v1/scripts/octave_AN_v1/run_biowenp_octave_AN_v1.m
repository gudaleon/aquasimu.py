% Octave driver for the MATLAB-vs-Python consistency check (biowenp_py_AN_v1).
% Runs the ORIGINAL, UNMODIFIED BioWEN-P code from the read-only FRDR snapshot on a SYNTHETIC forcing:
%   - the hydrologic loop is evaluated from the text of ModelMainloop.m (PET block + main loop + truncation),
%   - the P model is the original function BioC_P_RTM.m (called directly), with Discharge.m and InflowEstimator.m.
% Usage: octave --no-gui -q run_biowenp_octave_AN_v1.m <snapshot_dir> <forcing.csv> <params.json> <out.mat>
args = argv();
snap = args{1}; fcsv = args{2}; fjson = args{3}; fout = args{4};
addpath(snap);
addpath(fileparts(mfilename('fullpath')));   % month.m helper
F = dlmread(fcsv, ',', 1, 0);    % columns: datenum, precip, temp
prm = jsondecode(fileread(fjson));
fn = fieldnames(prm);
global A_BC dfm_fil dfm n_fm rho_fm b_Sr b_bs n_bs alpha_POP alpha_PIP alpha_TDP beta_Microbial beta_Labile beta_Stable beta_SorbedP beta_MineralP beta_RecalcitrantP ...
    beta_SRP C_overflow lambda alpha_L Kf Pu_Depth I_Pu_max C_SRP_Pu_min Km_Pu k_LF k_ad Qmax_ad kL_ad k_mi Um Km_Im ...
    f_debris k_Death k_Hy lambda_s d_s
for k = 1:numel(fn)
  eval([fn{k} ' = prm.(fn{k});']);
end
dfm_fil = dfm; A_BMP = A_BC + A_PP; b_bs = hfm_u - 1.6; d_s = (hfm_u - b_Sr - dfm)*100;
Sr_crit = u - b_Sr; Qd_Sr_Net = []; ElmDrive = 1; dt = 1; Timestep = 1;
Model_Output = struct('Time', F(:,1), 'Precip', F(:,2), 'Temp', F(:,3));
Model_Output.Precip(isnan(Model_Output.Precip)) = 0;
Model_Output.Temp(isnan(Model_Output.Temp)) = 3;
stime = length(Model_Output.Precip);
src = fileread(fullfile(snap, 'ModelMainloop.m'));
i1 = strfind(src, 'if PETMethod == 0 || PETMethod == 1'); i2 = strfind(src, '% Calculate daily max and min');
eval(src(i1(1):i2(1)-1));
j1 = strfind(src, '% % Initiate vectors for simulated variables'); j2 = strfind(src, '% Remove index for unit hydrograph');
eval(src(j1(1):j2(1)-1));
Model_Output = BioC_P_RTM(Model_Output.s, Model_Output.C_TP_in, 0, Model_Output.I, dt, stime, dx, Macro_percent, SorptionEquilibrium, Model_Output);
O = struct();
keep = {'s','I','ET','L','Sp','Qo','Qd','Sr','Q_Ex','Q_in','Q_in_d','C_TP_in','Ss','C_macro','C_macro_POP','C_macro_PIP', ...
        'qe_ad','qe_mineral','qe_Recalcitrant','M_POP_Labile','M_POP_Stable','Plant_P','TPflux','TDPflux','Cout_TP', ...
        'Cout_TDP','C_Sr_TP','Ce_TDP_FM','Ce_PP_FM','TPExflux','Precip_flux','Adsorbed_flux'};
for k = 1:numel(keep)
  O.(keep{k}) = Model_Output.(keep{k});
end
save('-v7', fout, '-struct', 'O');
printf('octave run done: %d steps\n', stime);
