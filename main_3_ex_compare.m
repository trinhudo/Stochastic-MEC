%------------------------------
% compare proposed algorithm with exhaustive search in term of convergence behavior and runtime
% -----------------------------
clear all
close all

tic
load('parameter_settings.mat')

rng('default')

noSearchAgents = 30;
params.maxIter = 13; %1500;
params.maxIter_woa = 13; %100;

%NoUsers = 2:7; % values of N
NoUsers = 13; % 6 % to check convergence curve, we can plot the convergence curve based on the shown results on Command Window

M_ul = 2;
M_dl = 2;
noBSs   = M_ul + M_dl;
noSubcs = 3;
params.noSubcs = noSubcs;
noAnten = 4;

noRealizations = 1;

doTol = 0;

% po: percentage offloading
% su: system utility
dbstop if error

% MECNOMA21 means : we started this work in 2021 with NOMA idea, but later
% improved with the stochatic idea.

po_MECNOMA21   = zeros(length(NoUsers), noRealizations);   % 5 x noReal matrix
su_MECNOMA21   = zeros(length(NoUsers), noRealizations);   % 5 x noReal matrix
time_MECNOMA21 = zeros(length(NoUsers), noRealizations);   % 5 x noReal matrix

po_EX   = zeros(length(NoUsers), noRealizations);   % 5 x noReal matrix
su_EX   = zeros(length(NoUsers), noRealizations);   % 5 x noReal matrix
time_EX = zeros(length(NoUsers), noRealizations);   % 5 x noReal matrix

for iN = 1:length(NoUsers)

    users_no = NoUsers(iN);
    %     name = sprintf('../Conver_behave/position_data/pos_BS_UEs_%dUE.mat', users_no);
    %     load(name);

    for iReal = 1:noRealizations
        fprintf('iReal:%i/%i    iN:%i/%i',iReal,noRealizations,NoUsers(iN),NoUsers(length(NoUsers)));

        UEs.total = [2 6 10];
        while UEs.total(2) ~= floor(UEs.total(3)/2) % force N_ul = N_dl trick to get average quicker
            [UE_BS, UEs, BS] = location_voronoi(users_no, M_ul, M_dl, 0);
            % UE_BS_   == N_active x M matrix % matrix of relation of UEs and SBSs
            % UEs == 1x1 struct
            %       UEs.active   == N_active_ue x 2 matrix == (N_ul + N_dl) x 2 matrix
            %       UEs.inactive == N_inactive x 2 matrix
            %       UEs.inBS     == 1 x N_active_ue  : SBS that covers the active UEs
            %       UEs.total    == 1 x 2 matrix == [N_ul N_dl N]
            %       UEs.d        == N_active x N_active : distances between active UEs
            % BS  == 1x1 struct
            %       BS.positions == N_sbs x 2 matrix
            %       BS.SBS       == N_sbs x 1 cell : save the positions of UEs that the SBS covers
            %       BS.total = [M_ul M_dl M]
            %       BS.d     == M x M == distances between SBSs and DL SBSs
        end
        N_ul = UEs.total(1);
        sys_voronoi{iN}.UE_BS = UE_BS;
        sys_voronoi{iN}.UEs = UEs;
        sys_voronoi{iN}.BS = BS;

        [ChannelGain, ~] = channelMod(UEs, BS, noAnten, noSubcs, logNormalMean, logNormalDeviation);
        % ChannelGain == struct with
        %       hArray 	  == N x M x K cell,
        %                  each cell is a L x 1 vector  == vector of channel gain
        %                   (each SBS has L antennas)
        %       h2h       == N x N x M x K matrix
        %                   ex: h2h(1,1,m,k) = ||h_{1m}^k||^2
        %                       h2h(1,2,m,k) = |h_{1m}^k'*h_{2m}^k|
        %       h_UE      == N_ul x N_dl x K matrix
        %       G_SBS     == M_ul x M_dl x K cell
        %                       each cell == L (ul) x L (dl) matrix
        % ~     	  == N x M matrix 	== distance from UEs to SBSs


        t = randi(800, 1);
        var.f_l = params.f_user(t: t+N_ul-1);
        T_l = params.C_n ./ var.f_l;
        E_l = params.kappa .* params.C_n .*(var.f_l) .^2;

        var.eta     = params.beta_t .* params.D_n ./ (T_l);
        var.theta   = params.beta_e .* params.D_n ./ (params.zeta .* E_l);

        var.Adet = 1;

        [var.lb_woa, var.ub_woa, var.P_SBS_min, var.P_SBS_max, fobj_woa, fobj_woa_dl, fobj_bwoa] = getFunctionDetails2('SIC_MEC', UEs, BS, UE_BS, noSubcs, ChannelGain, params, var);
        %   function in ..\

        % Exhaustive search
        fprintf("\n Exhautive search \n")
        [leader_score_bwoa, leader_pos_bwoa, time] = exhaustive2(UEs, BS, UE_BS, fobj_bwoa, fobj_woa, fobj_woa_dl, ChannelGain.h2h, params, var);
        % function in ..\

        po_EX(iN, iReal) = sum(sum(leader_pos_bwoa(1:N_ul,:) ))/N_ul;
        su_EX(iN, iReal) = leader_score_bwoa;
        time_EX(iN, iReal) = time;
        leader_score_bwoa
        var.ex_lead = leader_score_bwoa;

        fprintf("BWOA \n")
        [BWOA_result, WOA_result, time] = BWOA4('WOA_SIC_MEC', doTol, UEs, BS, UE_BS, fobj_bwoa, fobj_woa, fobj_woa_dl, ChannelGain.h2h, params, var);
        % function in ..\WOA_voronoi
        po_MECNOMA21(iN, iReal) = sum(sum(BWOA_result.leader_pos(1:N_ul,:)))/users_no;
        su_MECNOMA21(iN, iReal) = BWOA_result.leader_score;
        time_MECNOMA21(iN, iReal) = time;

        BWOA.curve{iN} = BWOA_result.conver_curve;
        BWOA_result.leader_score


    end
end
BWOA.po = mean(po_MECNOMA21, 2);
BWOA.su = mean(su_MECNOMA21, 2);
BWOA.time = mean(time_MECNOMA21, 2);

EX.po = mean(po_EX, 2);
EX.su = mean(su_EX, 2);
EX.time = mean(time_EX, 2);

% save('results\Script_compare.mat', 'BWOA', 'EX', 'NoUsers','sys_voronoi', 'noBSs', 'noSubcs');

