% This script sets the values of input parameters of the system
% you can change their values here, then save to parameters2.mat
% ------------------------------------------------------------------

params.noSearchAgents = 30;
params.noAnten = 4;

params.logNormalMean = 0;
params.logNormalDeviation = 8.0;

params.noRealizations = 13; %200; %300;

params.beta_t = 0.5;
params.beta_e = 1 - params.beta_t;
params.beta = [params.beta_t params.beta_e];


params.n0 = db2lin(-114 - 30);
params.B_k = 1e6;

params.p_min = 1e-8;
params.p_max = 0.25;
params.P_SBS_max = 39.81; % 46 dBm
params.P_SBS_min = 0.25*10^(-3);

params.f0 = 1e9* 8;
params.D_n = 1*420e3;
params.C_n = 1000e6;
params.kappa = 5e-27;
params.zeta = 1;

params.lamda = 1e14;
params.nu = 1e14;
params.P_tol = 1.001;

params.maxIter_woa = 300;
params.maxIter = 1500; 		%1500
% params.noSubcs = 5;
% params.noBSs   = 5;

params.Adet = 1;

params.f_local = 1e9*[0.5 0.8 1];
params.f_user = zeros(1000, 1);
for i = 1:1000
    params.f_user(i) = params.f_local(randi(length(params.f_local), 1));
end

params.gam_dl_th = 1;
params.R_th      = 2.6e7; % threshold to compute DL utility
