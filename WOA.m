% ---------------------------------------------------------
% WOA for power transmision P == N x 1 matrix
%----------------------------------------------------------

% Output:
% leaderScore: value of obj function after this code == double
% leaderPos == N x 1 matrix
% convergenceCurve == 1 x maxIter  matrix = value of obj function after each iteration

function [leaderScore, leaderPos, convergenceCurve, cci_SBS] = WOA(noSearchAgents, N_ul, maxIter, var, fobj, Posi_P, X_)
% Input:
% noSearchAgents: number of whales
% N_ul  = number of UL UEs
% var == struct with
%       ub_woa == N_ul x 1 matrix == upper bound transmit power p_i^{max}
%       lb_woa == N_ul x 1 matrix == lower bound transmit power p_i^{min}
% Posi_P == N_dl x M_dl matrix == power allocation for DL
% X_ == (N_ul + M_dl) x K matrix == association matrix

% X = X_(1:N_ul, :);  % N_ul x K matrix

leaderPos   = zeros(N_ul, 1);
leaderScore = inf;

leader_score_pre = leaderScore;

convergenceCurve = zeros(1, maxIter);
noSubcs = size(X_,2);
cci_SBS = cell(1, noSubcs);

% ======================== Initialization =================
%     ub = ub/10^2;
%     lb = lb/10^4;
% If each variable has a different lb and ub
posi_p = 0.5*ones(N_ul, noSearchAgents).*(var.ub_woa - var.lb_woa)/10^2 + var.lb_woa; %rand(noUsers, noSearchAgents).*(ub - lb) + lb;
%         posi_p = ones(noUsers, noSearchAgents).*(ub); %rand(noUsers, noSearchAgents).*(ub - lb) + lb;

% ======================== Loop ===========================
% Loop counter
t = 0;
todoTol = 1; % =0 to run all iteration
delta = 1e-4;
flag = 0;

% Main loop
while t < maxIter && flag < 10

    % Return back the search agents that go beyond the boundaries of the search space
    tmp = posi_p;
    flag4lb = tmp < var.lb_woa;
    flag4ub = tmp > var.ub_woa;
    posi_p = tmp.*(~(flag4lb + flag4ub)) + var.lb_woa .*flag4lb + var.ub_woa .*flag4ub;

    % Calculate objective function for each search agent
    for i = 1:noSearchAgents
        [fitness, cci_SBS_tmp] = fobj(posi_p(:, i), Posi_P, X_);
        % posi_p(:, i) == N_ul x 1 matrix _ matrix of ofloading power
        % P_SBS        == N_dl x M_dl matrix _ matrix of broadcasting power
        % X_           == (N_ul + M_dl) x K binary matrix _ matrix of association

        % Update the leader
        if fitness < leaderScore
            leaderScore = fitness;
            leaderPos = posi_p(:, i);
            cci_SBS = cci_SBS_tmp;
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

        for n = 1:N_ul
            % follow the shrinking encircling mechanism or prey search
            if p < 0.5
                % search for prey (exploration phase)
                if abs(A) >= 1
                    randLeaderIndex = floor(noSearchAgents*rand + 1);
                    X_rand = posi_p(:, randLeaderIndex); 		% -> X_rand == N_ul x 1 matrix
                    D_X_rand = abs(C*X_rand(n) - posi_p(n,i)); 	% double
                    posi_p(n, i) = X_rand(n) - A*D_X_rand;
                elseif abs(A) < 1
                    D_Leader = abs(C*leaderPos(n) - posi_p(n, i));   	% D_Leader==double %% leaderPos == N_ul x 1
                    posi_p(n,i) = leaderPos(n) - A*D_Leader;
                end
            elseif p >= 0.5
                distance2Leader = abs(leaderPos(n) - posi_p(n,i));
                posi_p(n,i) = distance2Leader*exp(b.*l).*cos(l.*2*pi) + leaderPos(n);
            end
        end
    end

    % increase the iteration index by 1
    t = t + 1;
    convergenceCurve(1,t) = leaderScore;

    if todoTol == 1 && leaderScore<10 && abs(leaderScore - leader_score_pre) < delta    && (t>150)
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
% uncomment and set breakpoint to see the figure

