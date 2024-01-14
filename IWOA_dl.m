% ---------------------------------------------------------
% WOA for power transmission allocation P == M x N matrix
%----------------------------------------------------------

% Output:
% leaderScore: value of obj function after this code == double
% leaderPos == M x N matrix
% convergenceCurve == 1 x maxIter  matrix = value of obj function after each iteration

function [leaderScore, leaderPos, convergenceCurve] = IWOA_dl(noSearchAgents, noUsers, noBSs, UE_BS_, maxIter, var, fobj, posi_p_ul, X)
% noSearchAgents: number of whales
% noUsers = N_dl  = number of UEs
% noBSs   = M_dl  = number of SBSs
% UE_BS_  = N x M
% P_SBS_max    = 1 x M = maximum transmit power of SBSs
% posi_p_ul   == N_ul x 1 == power allocation for UL
% X      == (N_ul + M_dl) x K matrix == association matrix

start_idx_m = max(1, size(UE_BS_,2) - noBSs + 1);
start_idx_n = max(1, size(UE_BS_,1) - noUsers + 1);
UE_BS = UE_BS_(start_idx_n:end, start_idx_m:end); % N_dl x M_dl

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
delta = 1e-4;
flag = 0;
PS_flag1 = 0; % population of search agents flags
PS_flag2 = 0;
PS_max = 40; %40;
PS_min = 10;
alpha = 0.25;
gamma_non = 20;

% Main loop
while t < maxIter && flag < 10
    if t > 3
        if PS_flag1 == 0 || PS_flag2 == 0
            n_inc  = 0;
        end
        if PS_flag2 == 2 % "increase"
            % add n_inc search agents
            n_inc = round(noSearchAgents * (PS_max - noSearchAgents)^2 / PS_max^2);
            if n_inc >0
                % compute distance
                distant = zeros(1, noSearchAgents); % save the distance of SAs to the best SA
                for jj = 1: noSearchAgents
                    distant(jj) = norm(posi_p(:,:,jj) - leaderPos);
                end
                [~, II] = sort(distant, 'ascend'); % II: index of SAs from nearest to furthest -> to the best SA

                SA_bests = zeros(1,n_inc); % 1 x n_inc vector containting indexes of n_inc best SAs from n_in groups

                if n_inc == 1
                    n_inc_ = n_inc +1; % if n_inc =1, we still need to get 2 SAs for generating the new SA
                else n_inc_ = n_inc;
                end
                % divide II into n_inc groups
                temp = [1:n_inc_, randi(n_inc_, 1, noSearchAgents -n_inc_)]; % random the index of groups, make sure each group appears at least once
                temp = temp(randperm(length(temp))); % shuffle elements in temp

                for ii = 1: n_inc_
                    group_ii = II(temp == ii);  % vector containing indexes of SAs that are in group ii
                    SA_bests(ii) = group_ii(1);        % add the best SA at the group to S
                end

                for ii = 1:n_inc
                    % pick 2 random SAs in SA_bests
                    SA_picked = randperm(n_inc_, 2); % pick 2 indexes of SA

                    % generate position of new SA
                    posi_SAnew = alpha^0.5* posi_p(:,:,SA_bests(SA_picked(1))) + (1-alpha)^0.5 *posi_p(:,SA_bests(SA_picked(2)));
                    posi_p(:,:,noSearchAgents +ii) = posi_SAnew;
                end
                noSearchAgents = noSearchAgents + n_inc;
            end
        end

        if PS_flag1 == 1 || PS_flag2==1 % "decrease"  % delete furthest SAs
            n_dec = round(noSearchAgents * (PS_max - noSearchAgents)^2 / PS_max^2);

            % compute distance
            distant = zeros(1, noSearchAgents);
            for jj = 1: noSearchAgents
                distant(jj) = norm(posi_p(:,:,jj) - leaderPos);
            end
            [~, II] = sort(distant, 'descend'); % II: index of SAs from furthest to nearest -> to the best SA
            posi_p(:,:, II(1:n_dec)) = [];        % delete the furthest SAs

            noSearchAgents = noSearchAgents - n_dec;
        end
    end
    % Return back the search agents that go beyond the boundaries of the search space
    tmp = posi_p;
    flag4lb = tmp < lb;
    flag4ub = tmp > ub;
    posi_p = tmp.*(~(flag4lb + flag4ub)) + lb.*flag4lb + ub.*flag4ub;

    % Calculate objective function for each search agent
    for i = 1:noSearchAgents
        fitness = fobj(posi_p_ul, posi_p(:,:, i), X);
        % posi_p_ul == N_ul x 1
        % posi_p(:,:,i) == N_dl x M_dl
        % X == (N_ul + M_dl) x K

        % Update the leader
        if fitness > leaderScore
            leaderScore = fitness;
            leaderPos = posi_p(:,:, i);  % N_dl x M_dl
        end
    end


    % a decreases linearly from 2 to 0
    % a = (1- t/(gamma_non * maxIter)) * (1+ 1/(1- gamma_non* t /maxIter));
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

        for n = 1:noUsers
            for m = 1:noBSs
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

    leader_pos_SA{t} = leaderPos;   % cell of N x 1
    %       == position of best SA in generation t
    if t >2
        if ((sum(sum(leader_pos_SA{t} ~= leader_pos_SA{t-1})) >0 ) && ...
                (sum(sum(leader_pos_SA{t-1} ~= leader_pos_SA{t-2})) >0) && ... % update 2 gen consecutively
                noSearchAgents > PS_min)
            PS_flag1 = 1; %"decrease";
        else
            PS_flag1 = 0;
        end
        if ((sum(sum((leader_pos_SA{t} ~= leader_pos_SA{t-1}) == 0 ))) && ... % not update 1 gen
                noSearchAgents == PS_max)
            PS_flag2 = 1; % "decrease";
        elseif ((sum(sum(leader_pos_SA{t} ~= leader_pos_SA{t-1}) ==0 )) && ...
                noSearchAgents < PS_max)
            PS_flag2 = 2;% "increase";
        else
            PS_flag2 = 0;
        end
        if noSearchAgents > PS_max
            PS_flag2 =1;
        end
    end

    if todoTol == 1 && leaderScore<10 && abs(leaderScore - leader_score_pre) < delta     && (t>150)
        flag = flag + 1;
        convergenceCurve = convergenceCurve(1, 1:t);
    else
        flag = 0;
    end
    %        fprintf('WOA iter:%i, leaderScore:%i, flag:%i\n', t, leaderScore, flag)
    leader_score_pre = leaderScore;

end
% plot(1:size(convergenceCurve,2),convergenceCurve);
% hold on;
% plot conver curve of WOA in one iter,
% uncomment and set breakpoint to se the figure

