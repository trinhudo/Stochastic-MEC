% Reference for original PSO code:
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

function [leaderScore, leaderPos, convergenceCurve] = PSO(noSearchAgents, noUsers, MaxIt, var, fobj, Posi_P, X)
% Input:
% noSearchAgents: number of whales
% noUsers = N_ul  = number of UL UEs
% var == struct with
%       ub_woa == N_ul x 1 matrix == upper bound transmit power p_i^{max}
%       lb_woa == N_ul x 1 matrix == lower bound transmit power p_i^{min}
% Posi_P == N_dl x M_dl matrix == power allocation for DL
% X == (N_ul + M_dl) x K matrix == association matrix


%% PSO Parameters
% MaxIt = 300;      % Maximum Number of Iterations
%     noSearchAgents = 100;   % nPop=100;        % Population Size (Swarm Size)
% PSO Parameters
w     = 1;           % Inertia Weight
wdamp = 0.99;        % Inertia Weight Damping Ratio
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
VelMax = 0.1*(var.ub_woa - var.lb_woa); % N_ul x 1
VelMin = - VelMax;      % N_ul x 1


% ======================== Initialization =================

empty_particle.Position      = []; % N_ul x 1 matrix
empty_particle.Cost          = [];
empty_particle.Velocity      = []; % N_ul x 1 matrix
empty_particle.Best.Position = [];
empty_particle.Best.Cost     = [];

particle = repmat(empty_particle, noSearchAgents, 1); % nSA x 1 matrix
% each element is a class

GlobalBest.Cost = inf;
leaderPos   = zeros(noUsers, 1);
leaderScore = inf;

leader_score_pre = leaderScore;
convergenceCurve = zeros(1, MaxIt);

todoTol = 1; % =0 to run all iteration
it      = 0;
delta   = 1e-6;
flag    = 0;

for i = 1 : noSearchAgents

    % Initialize Position
    particle(i).Position = 0.5*rand(noUsers, 1).*(var.ub_woa - var.lb_woa)/100 + var.lb_woa;
    % N_ul x 1 matrix

    % Initialize Velocity
    particle(i).Velocity = zeros(noUsers,1); % N_ul x 1

    % Evaluation
    particle(i).Cost    = fobj(particle(i).Position, Posi_P, X);

    % ??      % Update Personal Best ??
    particle(i).Best.Position = particle(i).Position;   % N_ul x 1
    particle(i).Best.Cost     = particle(i).Cost;       % double

    % Update Global Best
    if particle(i).Best.Cost < GlobalBest.Cost

        GlobalBest  = particle(i).Best; %  struct with .Cost == double
        %     and .Position  == N x 1
        leaderScore = GlobalBest.Cost;
        leaderPos   = particle(i).Best.Position;
    end

    % Evaluation
    particle(i).Cost     = fobj(particle(i).Position, Posi_P, X);
end

%% PSO Main Loop
while (it < MaxIt && flag<10)
    it = it+1;

    for i = 1 : noSearchAgents

        % Update Velocity
        particle(i).Velocity = (w*particle(i).Velocity ...
            + c1*rand(noUsers,1).*(particle(i).Best.Position - particle(i).Position) ...
            + c2*rand(noUsers,1).*(GlobalBest.Position - particle(i).Position));

        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity, VelMin); % N_ul x 1
        particle(i).Velocity = min(particle(i).Velocity, VelMax);

        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;

        % Velocity Mirror Effect
        IsOutside = (particle(i).Position<var.lb_woa | particle(i).Position> var.ub_woa);
        particle(i).Velocity(IsOutside) = - particle(i).Velocity(IsOutside);

        % Apply Position Limits
        particle(i).Position = max(particle(i).Position, var.lb_woa); % N_ul x 1
        particle(i).Position = min(particle(i).Position, var.ub_woa);

        % Evaluation
        particle(i).Cost = fobj(particle(i).Position, Posi_P, X);
        % particle(i).Position == N_ul x 1 matrix _ matrix of ofloading power
        % P_SBS        == N_dl x M_dl matrix _ matrix of broadcasting power
        % X_           == (N_ul + M_dl) x K binary matrix _ matrix of association

        % Update Personal Best
        if particle(i).Cost < particle(i).Best.Cost

            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost     = particle(i).Cost;

            % Update Global Best
            if particle(i).Best.Cost < GlobalBest.Cost

                GlobalBest = particle(i).Best; % struct with .Position == N x 1
                %             .Cost     == double
                leaderScore         = GlobalBest.Cost;
                leaderPos           = GlobalBest.Position;

            end

        end

    end

    convergenceCurve(it) = GlobalBest.Cost;

    %         disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(convergenceCurve(it))]);

    if todoTol == 1 && leaderScore<10 && abs(leaderScore - leader_score_pre) < delta    && (it>150)
        flag = flag + 1;
        convergenceCurve = convergenceCurve(1, 1:it);
    else
        flag = 0;
    end

    leader_score_pre = leaderScore;

    w = w * wdamp;

end
%     plot(1:length(convergenceCurve), convergenceCurve);
%     hold on
end