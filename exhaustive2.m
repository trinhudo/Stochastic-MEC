%--------------------------------------------------------------------------------
% Algorithm to solve ODSTCA pb by
% using EXHAUSTIVE SEARCH to solve SA pb, instead of using BWOA. (still using WOA to solve TPC problem)
%--------------------------------------------------------------------------------

% Output:
% leader_score_bwoa	== double		== obtained value of maximum U
% leader_pos_bwoa	== (N_ul + M_dl) x K
% time				== double

function [leader_score_bwoa, leader_pos_bwoa, time] = exhaustive2(UEs, BS, UE_BS, fobj_bwoa, fobj_woa, fobj_woa_dl, h2h, params, var)
% Input:

tic
sa 	= zeros(UEs.total(1)+BS.total(2), params.noSubcs); % (N_ul + M_dl) x K matrix
n 	= 1;
cnt = 0;
noSubcs = params.noSubcs;

N_ul = UEs.total(1);
N_dl = UEs.total(2);
M_dl = BS.total(2);

WOA_rs    = 0;
WOA_rs_dl = 0;

leader_pos_bwoa 	= zeros(size(sa));
leader_score_bwoa 	= -inf;
bwoa                = -inf;

% leader_score_woa
pos_woa = zeros(N_ul, 1);
leader_pos_woa = zeros(N_ul, 1);
pos_woa_dl = zeros(N_dl, M_dl);
leader_pos_woa_dl = zeros(N_dl, M_dl);
leader_pos_woa_dl(:,:) = params.P_SBS_min;


phi = @(y,a,x,eta) y*log2(1 + a*x) - (a/log(2))*(eta + y*x)/(1 + a*x);
fmin = @(y,a,x,eta) (eta+y*x)/(params.B_k*log2(1 + a*x));
varepsilon = 1e-5 ;

TRY(n);

    function [] = solution(sa1)
        flag = 0;
        % COMPUTING RESOURCE ALLOCATION
        fitbwoa = fobj_bwoa(sa1);

        %         % condition 1
        % 		if fitbwoa > 1e2 || fitbwoa <= leader_score_bwoa
        % 				flag = 1;
        % 		end
        %
        % 		% condition 2
        %         if (flag == 0)
        %             WOA_tmp = 0;
        %                 p_tmp = zeros(N_ul,1);
        %                 % p_tmp = zeros(noUsers,1);
        %                 % bisection
        %                 for n1 = 1: (N_ul)
        %                     m = find(UE_BS(n1,:)>0); % SBS that covers UE n
        %                     for k = 1:noSubcs
        %                         if sa(n1, k) == 0
        %                             continue;
        %                         end
        %                         if phi(var.theta(n1), h2h(n1, n1, m, k)/(params.n0), params.p_max, var.eta(n1)) <= 0
        %                             p_tmp(n1) = params.p_max;
        %                         else
        %                             p_s = params.p_min;
        %                             p_t = params.p_max;
        %                             while (abs(p_t - p_s) > varepsilon)
        %                                 p_l = (p_t + p_s)/2;
        %                                 if phi(var.theta(n1), h2h(n1, n1, m, k)/(params.n0), p_l, var.eta(n1)) <= 0
        %                                     p_s = p_l;
        %                                 else
        %                                     p_t = p_l;
        %                                 end
        %                             end
        %                             p_tmp(n1) = (p_s + p_t)/2;
        %                         end
        %                         WOA_tmp = WOA_tmp + fmin(var.theta(n1), h2h(n1, n1, m, k)/(params.n0), p_tmp(n1), var.eta(n1));
        %                     end
        %                 end
        %             if (fitbwoa - WOA_tmp) <= leader_score_bwoa
        %                 flag = 1;
        %             end

        if flag == 0
            % TRANSMIT POWER ALLOCATION
            if N_ul > 0
                [WOA_rs, pos_woa, ~] = WOA(params.noSearchAgents, ...
                    N_ul, params.maxIter_woa, var, fobj_woa, leader_pos_woa_dl, sa1);
            end

            if N_dl>0
                if (cnt==0)
                    leader_pos_woa = pos_woa;
                end
                [WOA_rs_dl, pos_woa_dl, ~] = WOA_dl(params.noSearchAgents, ...
                    N_dl, M_dl, UE_BS, params.maxIter_woa, var, fobj_woa_dl, leader_pos_woa, sa1);
            end
        end


        bwoa = fitbwoa - WOA_rs + WOA_rs_dl; % double

        cnt = cnt + 1;
        if bwoa > leader_score_bwoa
            leader_score_bwoa = bwoa;
            % make sure DL association to be reasonable
            fl_ = sum(UE_BS,1)>0;         % 1 x M
            fl_ = fl_(1,BS.total(1)+1 :end)';  % M_dl x 1
            sa1((N_ul+1):end, :) = sa((N_ul+1):end, :) .* fl_;
            leader_pos_bwoa   = sa1;
            leader_pos_woa 	  = pos_woa;
            leader_pos_woa_dl    = pos_woa_dl;
            fprintf('exhaustive leader: %i\n', bwoa);
        end
    end

    function [] = TRY(n)
        for j = 0:params.noSubcs
            if j ~= 0
                sa(n, j) = 1;
            end
            if n == N_ul
                if M_dl == 0
                    sa1 = sa;
                    solution(sa1);
                elseif M_dl>0
                    TRY2(1);
                end
            else
                TRY(n + 1);
            end
            if j ~= 0
                sa(n, j) = 0;
            end
        end
    end

    function [] = TRY2(m)
        for k = 1:params.noSubcs
            sa(N_ul+m, k) = 1;
            if m == M_dl
                sa1 = sa;
                % make sure DL association to be reasonable
                fl_ = sum(UE_BS,1)>0;         % 1 x M
                fl_ = fl_(1,BS.total(1)+1 :end)';  % M_dl x 1
                sa1((N_ul+1):end, :) = sa1((N_ul+1):end, :) .* fl_;
                solution(sa1);
            else
                TRY2(m+1);
            end
            sa(N_ul+m, k) = 0;
        end
    end


toc;
time = toc;
end