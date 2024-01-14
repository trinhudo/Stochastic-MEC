%==========================================================================
% This function returns the channel gain matrix in linear and UEs-BSs distances
%=================================================================================

% Output:
% ChannelGain == struct with
% hArray 	  == N x M x K cell,
%                  each cell is a L x 1 vector  == vector of channel gain
%                   (each SBS has L antennas)
% h2h         == N x N x M x K matrix
%                   ex: h2h(1,1,m,k) = ||h_{1m}^k||^2
% h2h(1,2,m,k)== |h_{1m}^k'*h_{2m}^k|
% h_UE        == N_ul x N_dl x K matrix
% G_SBS       == M_ul x M_dl x K cell
% each cell   == L (ul) x L (dl) matrix
% r_nm 		  == N x M matrix  == distance from UEs to SBSs

function [ChannelGain, r_nm] = channelMod(UEs, BS, noAnten, noSubcs, logNormalMean, logNormalDeviation)
% UEs == 1x1 struct
% (first N_ul elements of UEs are UL UEs, last N_dl elements of UEs are DL UEs)
%       UEs.active   == N_active x 2 matrix
%                               (1st col == x-coordinate
%                                2nd col == y-coordinate)
%       UEs.inactive == N_inactive x 2 matrix
%                               (1st col == x-coordinate
%                                2nd col == y-coordinate)
%       UEs.inBS     == 1 x N_active  : SBS that covers the active UEs
%                              example: UEs.inBS(2) = 4 means...
%                                         ...UE 2 in coverage of SBS 4
%       UEs.d        == N_active x N_active : distances between active UEs
%                      == (N_ul + N_dl) x (N_ul + N_dl)
% BS  == 1x1 struct
% (first M_ul elements of SBSs are UL SBSs, last M_dl elements of SBSs are DL SBSs)
%       BS.positions == N_sbs x 2 matrix
%                               (1st col == x-coordinate
%                                2nd col == y-coordinate)
%       BS.SBS       == N_sbs x 1 cell : save the positions of UEs that the SBS covers
%                       example: BS.SBS{1} == [150 100;
%                                              120 200;
%                                             -125 100]
%                                     --> SBS1 covers the UEs at (150,100),
%                                                     (120,200), (-125,100)
% noSubcs  == double == number of subchannels
% logNormalMean = 0; logNormalDeviation = 8;

% only active UEs are considered
N = size(UEs.active, 1); % number of UEs
M = size(BS.positions, 1); % number of SBSs
K = noSubcs;
L = noAnten;

M_ul = BS.total(1);
M_dl = BS.total(2);

r_nm = zeros(N, M); 	% N x M matrix == distances from UEs to SBSs
% the row is the index of UEs and the col is the index of SBSs

% Calculate the value of r_nm
for i = 1: M
    r_nm(:,i) = sqrt(sum((BS.positions(i,1:2)-UEs.active).^2,2) + 15^2);  % [m] - SBS's height =15m
end

%% return the real distance between UEs and BSs
dArray = zeros(N,M,K); % transform 2D matrix into 3D matrix to calculate hArray
for k = 1:K
    dArray(:,:,k) = r_nm;   % [m]
end
% Caculate the channel gain between the eNodeB and mobile users
% Channel gain hArray represents pathloss and log-normal shadowing
hArray    = cell(N,M,K);    % each cell is a L x 1 vector
for n = 1:N
    for m = 1:M
        for k = 1:K
            hArray{n,m,k} = zeros(L,1);
            C = normrnd(logNormalMean, logNormalDeviation, 25, L);
            shadowing = mean(C)'; % L x 1

            % The path loss model for small cells is h = 140.7 + 36.7*log10(d) (small cell)
            %hArray{n,m,k} = -140.7 - 36.7.*log10(dArray(n,m,k)/1000) + shadowing; % L x 1
            % uplink

            %hArray{n,m,k} = db2lin(hArray{n,m,k} - 5);

            % Return the channel gain (normalized by the noise which is 95 dBm - 30 = 65 dB) in linear
            % hArray = db2lin(hArray - 5); % UL
            % hArray = db2lin(hArray);

            hArray{n,m,k} = db2pow(-22.7 - 26*log10(2.4) - 36.7*log10(dArray(n,m,k)/1000)) + shadowing; %

        end
    end
end

h2h = zeros(N,N,M,K); % ex: h2h(1,1,m,k) = ||h_{1m}^k||^2
% h2h(1,2,m,k) = |h_{1m}^k'*h_{2m}^k|

for k = 1:K
    h_k = hArray(:,:,k); % N x M cell
    for m = 1:M
        h_mk = h_k';     % M x N cell
        h_mk = h_mk(m,:); % 1 x N cell of Lx1 matrix
        h_mk = cell2mat(h_mk); % L x N matrix
        h2h(:,:,m,k) = h_mk'*h_mk;
    end
end
h2h = sqrt(h2h);


N_ul = UEs.total(1);
N_dl = UEs.total(2);
h_UE = zeros(N_ul, N_dl, K); % channel between UEs
for n_ul = 1 : N_ul
    for n_dl = 1:N_dl
        for k = 1: K
            C_ = normrnd(logNormalMean, logNormalDeviation, 25, 1);
            shadowing = mean(C_);
            h_UE(n_ul,n_dl,k) = db2pow(-22.7 - 26*log10(2.4) - 36.7*log10(UEs.d(n_ul,N_ul + n_dl)/1000)) + shadowing;
        end
    end
end

G_SBS = cell(M_ul, M_dl, K); % each cell is a L x L matrix
for m_ul = 1:M_ul
    for m_dl = 1:M_dl
        for k = 1:K
            G_SBS{m_ul,m_dl,k} = zeros(L,L);
            C = normrnd(logNormalMean, logNormalDeviation, L,L,25);
            shadowing = mean(C,3); % L x L

            % The path loss model for small cells is h = 140.7 + 36.7*log10(d) (small cell)
            %hArray{n,m,k} = -140.7 - 36.7.*log10(dArray(n,m,k)/1000) + shadowing; % L x 1
            % uplink

            %hArray{n,m,k} = db2lin(hArray{n,m,k} - 5);

            % Return the channel gain (normalized by the noise which is 95 dBm - 30 = 65 dB) in linear
            % hArray = db2lin(hArray - 5); % UL
            % hArray = db2lin(hArray);

            G_SBS{m_ul,m_dl,k} = db2pow(-22.7 - 26*log10(2.4) - 36.7*log10(BS.d(m_ul,m_dl+M_ul)/1000)) + shadowing; % LxL

        end
    end
end

ChannelGain.h2h    = h2h;
ChannelGain.hArray = hArray;
ChannelGain.h_UE   = h_UE;
ChannelGain.G_SBS  = G_SBS;
ChannelGain.r_nm   = r_nm;

% %--------------Plotting the locations----------------%
%  LabledPlotting(xx_n, yy_n, teta,tera, UE_position, eNB_position, cellRadiusMin, cellRadiusMax)
end
