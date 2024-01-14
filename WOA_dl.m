% ---------------------------------------------------------
% WOA for power transmission allocation P == M_dl x N_dl matrix
%----------------------------------------------------------

% Output:
% leaderScore: value of obj function after this code == double
% leaderPos == N_dl x M_dl matrix
% convergenceCurve == 1 x maxIter  matrix = value of obj function after each iteration

function [leaderScore, leaderPos, convergenceCurve] = WOA_dl(noSearchAgents, noUsers, noBSs, UE_BS_, maxIter, var, fobj, posi_p_ul, X_)
% noSearchAgents: number of whales
% noUsers = N_dl  = number of UEs
% noBSs   = M_dl  = number of SBSs
% UE_BS_  = N x M
% P_SBS_max    = 1 x M = maximum transmit power of SBSs
% posi_p_ul   == N_ul x 1 == power allocation for UL
% X_      == (N_ul + M_dl) x K matrix == association matrix

% start_idx_M = max(1, size(X_,1) - noBSs + 1);
% X = X_(start_idx_M:end, :); % == M_dl x K

start_idx_m = max(1, size(UE_BS_,2) - noBSs + 1);
start_idx_n = max(1, size(UE_BS_,1) - noUsers + 1);
UE_BS = UE_BS_(start_idx_n:end, start_idx_m:end); % N_dl x M_dl

% maxIter     = 300; % 150
leaderPos   = zeros(noUsers, noBSs); % N_dl x M_dl
leaderScore = -inf;

leader_score_pre = leaderScore;

convergenceCurve = zeros(1, maxIter);

% ======================== Initialization =================

% new channel model
lb = var.P_SBS_min.*UE_BS;         % N_dl x M_dl == lower bound of power transmission of M SBSs to N UEs
ub = var.P_SBS_max.*UE_BS;         % N_dl x M_dl == upper bound

% % old channel model
% lb = P_SBS_min.*UE_BS / 10^32;         % N x M == lower bound of power transmission of M SBSs to N UEs
% ub = P_SBS_max.*UE_BS / 10^20;         % N x M == upper bound

posi_p = zeros(noUsers, noBSs, noSearchAgents);	% N_dl x M_dl x noSA matrix
% posi_p == N_dl x M_dl x noSA
for nSA = 1:noSearchAgents
    posi_p(:,:,nSA) = UE_BS .* (1/noUsers *rand(size(posi_p(:,:,nSA))).*(ub-lb) + lb);
end

% ======================== Loop ===========================
% Loop counter
t = 0;
todoTol = 1; % =0 to run all iteration
delta = 1e3;
flag = 0;

% Main loop
while t < maxIter && flag < 10

    % Return back the search agents that go beyond the boundaries of the search space
    tmp = posi_p;
    flag4lb = tmp < lb;
    flag4ub = tmp > ub;
    posi_p = tmp.*(~(flag4lb + flag4ub)) + lb.*flag4lb + ub.*flag4ub;

    % Calculate objective function for each search agent
    for i = 1:noSearchAgents
        fitness = fobj(posi_p_ul, posi_p(:,:, i),  X_);
        % posi_p_ul == N_ul x 1
        % posi_p(:,:,i) == N_dl x M_dl
        % X_ == (N_ul + M_dl) x K

        % Update the leader
        if fitness > leaderScore
            leaderScore = fitness;
            leaderPos   = posi_p(:,:, i);  % N_dl x M_dl
        end
    end


    % a decreases linearly from 2 to 0
    a = 2 - t*(2/maxIter);
    % a2 linearly decreases from -1 to -2 to calculate t
    a2 = -1 + t*(-1/maxIter);

    % Update the position of each search agents
    for i = 1:noSearchAgents
        r1 = rand();
        r2 = rand();

        A = 2*a*r1 - a;
        C = 2*r2;

        % parameters for spiral updating position
        b = 1;
        l = (a2 - 1)*rand + 1;

        p = rand();

        for n = 1:noUsers   % 1:N_dl
            for m = 1:noBSs   % 1:M_dl
                if UE_BS(n,m) ~= 1
                    continue  % only consider transmit power of available UE-BS association
                end
                % follow the shrinking encircling mechanism or prey search
                if p < 0.5
                    % search for prey (exploration phase)
                    if abs(A) >= 1
                        randLeaderIndex = floor(noSearchAgents*rand + 1);
                        X_rand = posi_p(:,:, randLeaderIndex); 		% -> X_rand == N x M matrix
                        D_X_rand = abs(C*X_rand(n,m) - posi_p(n,m,i)); 	% double
                        posi_p(n,m, i) = X_rand(n,m) - A*D_X_rand;
                    elseif abs(A) < 1
                        D_Leader = abs(C*leaderPos(n,m) - posi_p(n,m, i));   	% D_Leader==double %% leaderPos == N x 1
                        posi_p(n,m,i) = leaderPos(n,m) - A*D_Leader;
                    end
                elseif p >= 0.5
                    distance2Leader = abs(leaderPos(n,m) - posi_p(n,m,i));
                    posi_p(n,m,i) = distance2Leader*exp(b.*l).*cos(l.*2*pi) + leaderPos(n,m);
                end
            end
        end
    end

    % increase the iteration index by 1
    t = t + 1;
    convergenceCurve(1,t) = leaderScore;

    if todoTol == 1 && (t>150) && abs(leaderScore - leader_score_pre) < delta
        flag = flag + 1;
        convergenceCurve = convergenceCurve(1, 1:t);
    else
        flag = 0;
    end
    %         fprintf('WOA iter:%i, leaderScore:%i, flag:%i\n', t, leaderScore, flag)
    leader_score_pre = leaderScore;

end
%     plot(1:size(convergenceCurve,2),convergenceCurve);
% plot conver curve of WOA in one iter,
% uncomment and set breakpoint to se the figure

