%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA102
% Project Title: Implementation of Particle Swarm Optimization in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
%
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
%
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

% each earch agent is a matrix of N x M
function [leaderScore, leaderPos, convergenceCurve] = PSO_dl(noSearchAgents, N_dl, M_dl, UE_BS_, MaxIt, var, fobj, posi_p_ul, X)

start_idx_m = max(1, size(UE_BS_,2) - M_dl + 1);
start_idx_n = max(1, size(UE_BS_,1) - N_dl + 1);
UE_BS = UE_BS_(start_idx_n:end, start_idx_m:end); % N_dl x M_dl


VarMin = var.P_SBS_min.*UE_BS;         % Lower Bound of Variables
% P_SBS_min.*UE_BS
VarMax = var.P_SBS_max.*UE_BS;         % Upper Bound of Variables
% P_SBS_max.*UE_BS;         % N_dl x M_dl == upper bound

%% PSO Parameters
% MaxIt = 300;      % Maximum Number of Iterations
% noSearchAgents = 30;   % nPop=100;        % Population Size (Swarm Size)
% PSO Parameters
w     = 1;            % Inertia Weight
wdamp = 0.99;     % Inertia Weight Damping Ratio
c1    = 1.5;         % Personal Learning Coefficient
c2    = 2.0;         % Global Learning Coefficient
% If you would like to use Constriction Coefficients for PSO,
% uncomment the following block and comment the above set of parameters.
% % Constriction Coefficients
% phi1=2.05;
% phi2=2.05;
% phi=phi1+phi2;
% chi=2/(phi-2+sqrt(phi^2-4*phi));
% w=chi;          % Inertia Weight
% wdamp=1;        % Inertia Weight Damping Ratio
% c1=chi*phi1;    % Personal Learning Coefficient
% c2=chi*phi2;    % Global Learning Coefficient
% Velocity Limits
VelMax = 0.1*(VarMax - VarMin); % N x M
VelMin = - VelMax;

%% Initialization
empty_particle.Position      = []; % N x M matrix
empty_particle.Cost          = [];
empty_particle.Velocity      = []; % N x M matrix
empty_particle.Best.Position = [];
empty_particle.Best.Cost     = [];

particle = repmat(empty_particle, noSearchAgents, 1); % nSA x 1 matrix
% each element is a class

GlobalBest.Cost = -inf;
leaderScore = -inf;
leaderPos   = zeros(N_dl, M_dl);
leader_score_pre = leaderScore;
convergenceCurve = zeros(1, MaxIt);

todoTol = 1; % =0 to run all iteration
it      = 0;
delta   = 1e3;
flag    = 0;

for i = 1 : noSearchAgents

    % Initialize Position
    particle(i).Position = 1/sum(sum(UE_BS))* UE_BS.* unifrnd(VarMin, VarMax, size(UE_BS));
    % N_dl x M_dl matrix

    % Initialize Velocity
    particle(i).Velocity = zeros(size(UE_BS));

    % Evaluation
    particle(i).Cost    = fobj(posi_p_ul, particle(i).Position, X);
    % posi_p_ul == N_ul x 1
    % particle(i).Position == N_dl x M_dl
    % X == (N_ul + M_dl) x K

    % ??      % Update Personal Best ??
    particle(i).Best.Position = particle(i).Position;   % N x M
    particle(i).Best.Cost     = particle(i).Cost;       % double

    % Update Global Best
    if particle(i).Best.Cost > GlobalBest.Cost

        GlobalBest  = particle(i).Best; %  struct with .Cost == double
        %     and .Position  == N x M
        leaderScore = GlobalBest.Cost;
        leaderPos   = particle(i).Best.Position;
    end

    % Evaluation
    particle(i).Cost     = fobj(posi_p_ul, particle(i).Position, X);
    % posi_p_ul == N_ul x 1
    % particle(i).Position == N_dl x M_dl
    % X_ == (N_ul + M_dl) x K
end

%% PSO Main Loop
while (it < MaxIt && flag<10)
    it = it+1;

    for i = 1 : noSearchAgents

        % Update Velocity
        particle(i).Velocity = UE_BS.* (w*particle(i).Velocity ...
            + c1*rand(size(UE_BS)).*(particle(i).Best.Position - particle(i).Position) ...
            + c2*rand(size(UE_BS)).*(GlobalBest.Position - particle(i).Position));

        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity, VelMin); % N x M
        particle(i).Velocity = min(particle(i).Velocity, VelMax);

        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;

        % Velocity Mirror Effect
        IsOutside = (particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside) = - particle(i).Velocity(IsOutside);

        % Apply Position Limits
        particle(i).Position = max(particle(i).Position, VarMin); % N x M
        particle(i).Position = min(particle(i).Position, VarMax);

        % Evaluation
        particle(i).Cost = fobj(posi_p_ul, particle(i).Position, X);
        % posi_p_ul == N_ul x 1
        % particle(i).Position == N_dl x M_dl
        % X_ == (N_ul + M_dl) x K

        % Update Personal Best
        if particle(i).Cost > particle(i).Best.Cost

            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost     = particle(i).Cost;

            % Update Global Best
            if particle(i).Best.Cost > GlobalBest.Cost

                GlobalBest = particle(i).Best; % struct with .Position == N x m
                %             .Cost     == double
                leaderScore         = GlobalBest.Cost;
                leaderPos           = GlobalBest.Position;

            end

        end

    end

    convergenceCurve(it) = GlobalBest.Cost;

    %         disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(convergenceCurve(it))]);

    if todoTol == 1 && (it>150) && abs(leaderScore - leader_score_pre) < delta
        flag = flag + 1;
        convergenceCurve = convergenceCurve(1, 1:it);
    else
        flag = 0;
    end

    leader_score_pre = leaderScore;

    w = w * wdamp;

end
% plot(1:length(convergenceCurve), convergenceCurve);
end