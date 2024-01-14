%------------------------------------------------------------------
% Obtain the objective function and search domain for WOA and BWOA.
% Version C4 in fraction form
% matrix in the dimension of N x M+1 x K
%------------------------------------------------------------------

% Output:
% lb_woa = N x 1 matrix = lower bound of power transmission of UEs
% ub_woa = N x 1 matrix = upper bound of power transmission of UEs
% fobj_woa = 'string' = @function name of woa
% fobj_bwoa = 'string' = @function name of bwoa

function [lb_woa, ub_woa, P_SBS_min, P_SBS_max, fobj_woa, fobj_woa_dl, fobj_bwoa] = getFunctionDetails2(F, UEs, BS, UE_BS, noSubcs, ChannelGain, params, var)
% h2h in ChannelGain

% F function name: F: depends on each simulation schemes {ARJOA, ODSTCA, IOJOA, OFDMA}
% network size: noUsers x noBSs x noSubcs == N x M x K
% UE_BS   == N x M matrix   == binary matrix of relation of UEs and BSs
% ChannelGain = struct with
% h2h     == N x N x M x K matrix
%               ex: h2h(1,1,m,k) = ||h_{1m}^k||^2
%                   h2h(1,2,m,k) = |h_{1m}^k'*h_{2m}^k|
% params == struct with
% users' power budget: [p_min, p_max]
% tolerant power for NOMA: P_tol == double
% thermal noise: n0
% B_k: bandwidth of subchannel
% f0 == 1 x M matrix:  total computing resource of BSs
% penalty factor for constraint dealing: <double> nu (here, in this code; but in pdf, they are nu and lambda)
% perference in time and energy: beta = [beta_t beta_e] == N x 2 matrix
% matrix defines offloading decision of IOJOA: Adet 	== N x 1 matrix
% var == struct with
% eta, theta, f_l: parameters of local computing (in (26))  == N x 1 matrix

%     noUsers = size(UEs.active, 1);   % N
%     noBSs   = size(BS.positions, 1); % M

noUsers = UEs.total(3);
N_ul  = UEs.total(1);
N_dl  = UEs.total(2);

M_ul  = BS.total(1);
M_dl  = BS.total(2);

h2h = ChannelGain.h2h;
G_SBS = ChannelGain.G_SBS;
hArray = ChannelGain.hArray;
h_UE = ChannelGain.h_UE;

noSubcs = params.noSubcs;

lb_woa = zeros(N_ul, 1);
ub_woa = zeros(N_ul, 1);
lb_woa(:) = params.p_min;
ub_woa(:) = params.p_max;

P_SBS_min = zeros(1, M_dl);
P_SBS_max = zeros(1, M_dl);
P_SBS_min(:) = params.P_SBS_min;  % of SBS
P_SBS_max(:) = params.P_SBS_max;

PMin = zeros(N_ul, 1); % of UE
PMax = zeros(N_ul, 1);

PMin(:) = params.p_min;
PMax(:) = params.p_max;

switch F
    case 'ALCA'
        fobj_woa = @FWOA_ul;
        fobj_woa_dl = @FWOA_dl;
        fobj_bwoa = @FBWOA_ALCA;
    case 'SIC_MEC'
        fobj_woa = @FWOA_ul;
        fobj_woa_dl = @FWOA_dl;
        fobj_bwoa = @FBWOA_SICMEC;
    case 'ARJOA'
        fobj_woa = @FWOA_ul;
        fobj_bwoa = @FBWOA_ARJOA;
        fobj_woa_dl = @FWOA_dl;
    case 'IOJOA'
        fobj_woa = @FWOA_ul;
        fobj_woa_dl = @FWOA_dl;
        fobj_bwoa = @FBWOA_IOJOA;
    case 'OFDMA'
        fobj_woa = @FWOA_ul;
        fobj_woa_dl = @FWOA_dl;
        fobj_bwoa = @FBWOA_OFDMA;
end

% Beamforming
%% Calculate beamforming matrix w_mn == Lx1 vector
W_beam = cell(1, M_dl); % 1 x M_dl cell, each cell is a L x N_m x K matrix
% N_m is the number of UEs in cell m
% hArray == NxMxK cell - each cell is 1 Lx1 vector
% Need to extract H_m == L x N_m  x K from hArray (N_m == number of UEs in cell m)

% H_incell = cell(1, noBSs); % 1 x M cell, each cell is L x N_m x K matrix
%  == channel gains of UEs in 1 cell

for m = 1 : M_dl % consider DL SBSs
    % get indexes of UEs in cell m
    UEs_m = UE_BS(:,m+ M_ul)>0;
    if sum(UEs_m > 0)
        % get channel gains of UEs in cell m
        temp     = hArray(UEs_m', m+ M_ul, :);     % N_m x 1 x K cell of L x 1 vector
        temp     = permute(temp, [2 1 3]);   % 1 x N_m x K cell of L x 1 vector
        H_incell = cell2mat(temp); %         % L x N_m x K matrix
        %     H_incell{m} = cell2mat(temp);  % L x N_m x K matrix

        for kk = 1:noSubcs
            W_beam{m} (:,:,kk) = H_incell(:,:, kk) / ((H_incell(:,:, kk)' * H_incell(:,:, kk)));
            % L x N_m matrix
        end
    end
end
UE_BS_dl = UE_BS((N_ul+1):end, (M_ul+1):end); % N_dl x M_dl
BS_broad = find(sum(UE_BS_dl, 1)>0); % 1 x ?? vector of broadcasting SBSs (SBS cells that have UEs)
% range 1:M_dl

    function o = FBWOA_SICMEC(X)
    	% X == (N_ul + M_dl) x K matrix = asociation matrix
        X_ul = X(1:N_ul, :);
        % == N_ul x K
        % channel allocation matrix for UL
        X_dl = X((N_ul+1):end, :);
        % M_dl x K
        % channel allocation matrix for DL

    	% constraint dealing
        % Constraint of UL
        A_n  = sum(X_ul,2); % == N_ul x 1 matrix
        fc2  = A_n -1;		% == constraint 2
        flag_fc2 = fc2>0;

        pnal_fc2 = sum(params.nu .* flag_fc2.*(fc2.^2));

        % Constraint of DL
        A_m   = sum(X_dl,2); % == M_dl x 1 matrix
        tol_UE_dl = sum(UE_BS_dl,1)'; % M_dl x 1 == (each element) number of UEs in that cell
        have_UE = tol_UE_dl>0; % flag to check that there is UE(s) in the cell
        fc_dl = abs(A_m -1); % == M_dl x 1 matrix
        flag_fc_dl = fc_dl>0;

        pnal_fc_dl = sum(params.nu .*flag_fc_dl .*(fc_dl.^2) .*have_UE);

     	% objective function
        % get from GRA (33)
        tmp = sqrt(params.beta(:, 1) .* var.f_l); % numerator in (22)
        % == N_ul x 1 matrix
        % UE_BS association in UL
        A_ul = UE_BS(1:N_ul, 1:M_ul);  % == N_ul x M_ul
        % associations between UL UEs and UL SBSs

        A_nm  = repmat(A_ul, 1, 1, noSubcs);    % == N_ul x M_ul x K
        A_nk_temp = reshape(X_ul, N_ul, 1, noSubcs); % == N_ul x 1    x K
        A_nk  = repmat(A_nk_temp, 1, M_ul); % == == N_ul x M_ul x K
        A_nmk = A_nm.* A_nk; % == N_ul x M_ul x K
        % UE-SBS-subc associations

        A  	= sum(A_nmk,3).*tmp;     % N_ul x M_ul matrix
        A1 	= (sum(A,1)).^2;  		 % 1 x M_ul matrix
        V_AF= sum(A1./ params.f0);
        A2  = sum(sum(A_nmk,3),2);

    	beta_sum = sum(sum(params.beta, 2).*A2); % (beta_t + beta_e)

    	o = beta_sum - V_AF - pnal_fc_dl - pnal_fc2; % obj function has other subtrahends (values of WOA and U_dl obj functions)
        % but we perform it in line 104 code BWOA.m

    end


    function [o, cci_SBS_] = FWOA_ul(P, P_SBS, X)
        % P       == N_ul x 1 matrix _ matrix of ofloading power
        % P_SBS   == N_dl x M_dl matrix _ matrix of broadcasting power
        % X       == (N_ul + M_dl) x K binary matrix _ matrix of association

        % cci_SBS_ == 1 x K cell, each cell is 1 x M_ul vector

        % calculate P_SBSdl
        % P_SBSdl == M_dl x K cell
        %           == each cell is 1 x L matrix
        %               == each element is | P^{SBS, 1/2}_{m',k}[l] |^2

        P_SBSdl = cell(M_dl,noSubcs); % M_dl x K cell
        % each cell is 1 x L vector
        % == power of SBSs on subchannels
        cci_SBS_ = cell(1, noSubcs);

        for mm = BS_broad % range 1:M_dl
            if (~ isempty(W_beam{mm})) % if DL cell mm is not empty
                idx = find(UE_BS(:, mm+M_ul) >0);
                P_SBS_m = P_SBS(idx-N_ul, mm); % N_m x 1
                P_SBS_m_sqrt = sqrt(P_SBS_m);
                %    each cell is N_m x 1 vector

                for k1 = 1: noSubcs
                    % W_beam == 1 x M_dl cell, each cell is a L x N_m x K matrix
                    % N_m is the number of UEs in cell m
                    P_SBSdl{mm, k1} = W_beam{mm}(:,:,k1) * P_SBS_m_sqrt; % L x 1
                    % (L x N_m) * (N_m x 1)
                    P_SBSdl{mm, k1} = P_SBSdl{mm, k1}'.^2; % 1 x L
                end
            end
        end

    	o = 0;

    	rs = 0;
    	% C3:
        fc3 = P - PMax;
        flag_fc3 = fc3 >= 0; %== G(f(x)) in (31)
        pnal_fc3 = sum(params.nu .* flag_fc3.*(fc3.^2));  %% penalty function for C3
        rs = rs + pnal_fc3;
        % mu = hArray/(n0*W); %% normalized channel gain

        % C4:
        % g4 = 0;
        for k = 1: noSubcs
            [UE_off,~] = find(X(1:N_ul,k)>0 ); % == ?1 x 1 vector
            % find UEs occupying subchannel k
            if (sum(UE_off)==0)
                continue 			% no UE offload via subchannel k --> go check k+1
            end
            BS_off = UEs.inBS(UE_off')';  % == ?1 x 1 vector
            % the corresponding UL SBSs
            BS_no = unique(BS_off);  % == ?2 x 1 vector
            % UL SBSs that have UEs using subchannel k

            Xk = zeros(UEs.total(1), BS.total(1)); % == N_ul x M_ul
            % matrix of UL associations
            % between UEs and SBSs on
            % subchannel k
            % DL SBSs using subchannel k
            BS_dl_k = find(X(N_ul+1:end, k) >0); % == ?4 x 1 vector == SBSs that have DL UEs
            % range 1:M_dl

            for ii = 1:length(UE_off)
                Xk(UE_off(ii), BS_off(ii)) = 1;
            end

            I_nmk = zeros(N_ul,1);  	% intercell and intracell interference
            % to UEs using subchannel k
            % N x 1 matrix

            gamma_k = zeros(N_ul, M_ul); % N_ul x M_ul matrix == SNR matrix

            % cci_SBS == co-channel interference from DL SBSs to UL SBSs
            cci_SBS = zeros(1, M_ul); % 1 x M_ul
            % cci_SBS(m) = 0 if there is no  DL SBS using subc k

            for m_idx = 1 : length(BS_no) % consider BSs that have offloading UEs, not all BSs
                BS_idx = BS_no(m_idx);

                % Calculate the cci from DL SBSs to UL SBS BS_idx
                for m_dl_k = BS_dl_k'  % BS_dl_k ranges in 1:M_dl
                    dbstop if error
                    temp_mat = abs(G_SBS{BS_idx, m_dl_k, k}).^2 .* P_SBSdl{m_dl_k, k};
                    % == L(ul) x L(dl) matrix (double)
                    cci_SBS(BS_idx) = cci_SBS(BS_idx) + sum(sum(temp_mat));
                end


                UE_off_m = find(Xk(:,BS_idx)>0);  % indexes of the UEs that offloading tasks to BS_off(m)
                % == find N_m,k == n's where n \in N_m,k
                % == ?3 x 1 vector
                % calculate the interference that UEs offloading to BS_idx suffer
                for i = 1:length(UE_off_m)		   % UE_off_m(i) = n \in N_m,k
                    flag_less = zeros(N_ul, 1); % N x 1 matrix  		%% flag to check whether UE (n) satisfy SIC condition
                    for n_ = 1:N_ul
                        flag_less(n_, 1) = h2h(n_,n_,BS_idx,k) <= h2h(UE_off_m(i),UE_off_m(i),BS_idx,k); % consider BS number BS_idx
                        % ||h_{nm}^k||^2 <= ||h_{UE(i)m}^k||^2
                    end
                    flag_less(UE_off_m(i)) = 0;

                    %interfer_k(UE_off_m(i)) = sum(flag_less.*P.*Hk(:,BS_idx)); % double
                    %                    I_nmk(UE_off_m(i),1) = flag_less'*(P.*(h2h(UE_off_m(i),:,BS_idx,k)'.^2/(h2h(UE_off_m(i),UE_off_m(i),BS_idx,k))));
                    %  P == N x 1 matrix
                    %  h2h(i,:,BS_idx,k) == 1 x N matrix
                    %  I_nmk == N x 1 matrix
                    I_nmk(UE_off_m(i),1) = (flag_less'*(P.*diag(h2h(1:N_ul, 1:N_ul, BS_idx,k))));

                    g4		= params.P_tol *(flag_less'*(P.*diag(h2h(1:N_ul, 1:N_ul, BS_idx,k)))) - P(UE_off_m(i))*h2h(UE_off_m(i), UE_off_m(i), BS_idx,k); %double
                    flag_g4 = g4 >0;
                    pnal_g4 = sum(10^16* params.nu *flag_g4*(g4 ^2)); 	% double
                    rs = rs + pnal_g4;
                end
                % I_nmk == N x 1 matrix == noise from UL UEs


                % Calculate transmission rate
                for bs = 1: M_ul
                    gamma_k(:, bs) = (P.*diag(h2h(1:N_ul, 1:N_ul,bs,k)))./(params.noAnten * params.B_k * params.n0 + I_nmk + cci_SBS(bs)); % N_ul x 1
                end
                % gamma_k == N_ul x M_ul matrix

                Rk = params.B_k *log2(1+gamma_k);  % N_ul x M_ul matrix

                Wk = sum((var.eta + var.theta .*P).*sum(Xk./Rk,2)); % W in (31)
                % (N_ulx1 + N_ulx1.* N_ulx1)____N_ulxM_ul ./ N_ulxM_ul

                rs = rs + Wk;
            end
            cci_SBS_{k} = cci_SBS; % 1 x M_ul
            % cci_SBS_ == 1 x K cell
        end

    	o = o + rs;
    end

%% DL
    function o = FWOA_dl(P_ul, P_SBS, X)
        % P_ul    == N_ul x 1 matrix _ matrix of transmission power
        % P_SBS   == N_dl x M_dl matrix _ matrix of broadcasting power
        % X       == (N_ul + M_dl) x K binary matrix _ matrix of association

        X_dl = X((N_ul+1):end ,:); % M_dl x K

        P_SBSdl = cell(M_dl,noSubcs); % M_dl x K cell

        for mm = BS_broad % range 1:M_dl
            idx = (UE_BS((N_ul+1):end, mm+M_ul) >0);
            if (sum(idx)>0)  % if N_m>0
                P_SBS_m = P_SBS(idx, mm); % N_m x 1
                P_SBS_m_sqrt = sqrt(P_SBS_m);
                %    each cell is N_m x 1 vector

                for k1 = 1: noSubcs
                    % W_beam == 1 x M_dl cell, each cell is a L x N_m x K matrix
                    % N_m is the number of UEs in cell m
                    P_SBSdl{mm, k1} = W_beam{mm}(:,:,k1) * P_SBS_m_sqrt; % L x 1
                    % (L x N_m) * (N_m x 1)
                    P_SBSdl{mm, k1} = P_SBSdl{mm, k1}'.^2; % 1 x L
                end
            end
        end


        % Constraint
        % Constraint l
        pnal_fc8 = 0;

        P_mn = P_SBS.*UE_BS_dl ; % N_dl x M_dl  % only consider the transmit power of SBS to its associating UEs
        for m = BS_broad % range 1: M_dl
            norm_w2  = zeros(1, length(W_beam{m}(1,:,1)) );    % 1 x N_m
            W_beam_m = 0;
            for k = 1: noSubcs
                W_beam_m = W_beam_m + X_dl(m, k) * W_beam{m}(:,:,k); % L x N_m
                % with m pre-defined, X(m,k) =1 for just 1 value of k
            end
            for n = 1: length(W_beam{m}(1,:,1)) % for n = 1 to N_m
                norm_w2(n) = norm(W_beam_m(:,n)) ^2; % 1 x N_m
                %  norm ^2 of w_mn
            end
            P_m = P_mn(:,m);     % N_dl x 1
            P_m = P_m(P_m>0);    % N_m x 1
            fc8_m = 0;
            for nn = 1: length(P_m) % for nn = 1 to N_m
                fc8_mn = P_m(nn) * norm_w2(nn);
                fc8_m  = fc8_m + fc8_mn;
            end
            fc8_m = fc8_m - P_SBS_max(m);
            flag_fc8m = fc8_m >0;
            pnal_fc8 = pnal_fc8 + params.nu .* flag_fc8m.*(fc8_m.^2);
        end

        penal_fc9 = 0;
        gamma_dl  = zeros(M_dl, N_dl); % == M_dl x N_dl
        R_dl      = zeros(M_dl, N_dl, noSubcs);  % M_dl x N_dl x K

        % Calculate cci from UL UEs to DL UEs
        % P_ul X UE_BS
        % h_UE == N_ul x N_dl x K matrix
        cci_dl = zeros(N_dl, 1);
        for nn = 1: N_dl
            % find which cell is it in
            sbs_m =  find(UE_BS_dl(nn,:)>0); % range 1:M_dl

            % find subc k it is occupying
            k_m = find(X_dl(sbs_m,:)>0); % X_dl == M_dl x K

            if (length(k_m) == 1)
                % find UL UEs also occupying subc k X
                ue_k = find(X(1:N_ul, k_m)>0);
                for ii = 1:length(ue_k)
                    cci_dl(nn) = cci_dl(nn) + P_ul(ue_k(ii)) * abs(h_UE(ue_k(ii), nn ,k_m))^2;
                end
            end
        end

        for k = 1:noSubcs
            Xk = X((N_ul+1):end,k);		% Xk == M_dl x 1 vector == association vector regarding to subchannel k

            % Find the BSs broadcasting via subchannel k
            BS_no  = find(Xk>0);  	%  ?1 x 1 vector
            %  == indexes of BSs that broadcasting via subchannel k
            if (isempty(BS_no))
                continue 			% no SBS broadcast via subchannel k --> go check k+1
            elseif length(BS_no) == 1
                % find UEs in cell
                UEs_m = find(UE_BS_dl(:,BS_no) >0);  % ?2 x 1 vector
                % UE_BS_dl == N_dl x M_dl

                % Constraint m
                for n = UEs_m'
                    fc9_mn     = params.B_k * params.n0 * params.gam_dl_th - P_SBS(n,BS_no);
                    flag_fc9mn = fc9_mn >0;
                    penal_fc9_mn = sum(sum(params.nu .* flag_fc9mn.*(fc9_mn.^2)));
                    penal_fc9    = penal_fc9 + penal_fc9_mn;
                    % if penal_fc9< 1*10^(-7)
                    %     penal_fc9 = 0;
                    % end

                    gamma_dl(BS_no,n)= P_SBS(n,BS_no) / (cci_dl(n) + params.B_k * params.n0);
                    R_dl(BS_no,n,k) = params.B_k .* log2(1+gamma_dl(BS_no,n));
                end
                % interf = 0; % interference
            elseif length(BS_no) > 1 % exist inter-cell interference
                for m = BS_no'
                    % find UEs in cell m
                    UEs_m = find(UE_BS_dl(:,m) >0);  % ?2 x 1 vector
                    % UE_BS_dl == N_dl x M_dl
                    % find cell j also use subchannel k; j <> m
                    BS_not_m = BS_no';
                    BS_not_m(BS_not_m == m) = [];  % 1 x ?3 vector

                    for  n = UEs_m'
                        inn = 0;
                        for j = BS_not_m
                            P_j = P_SBS(:,j);                % N   x 1 vector
                            P_j = P_j(UE_BS_dl(:,j) >0);    % N_j x 1 vector
                            % change the dimension to N_j x 1 to be in
                            % harmony with W{j}(:,:,k) == L x N_j
                            for i = 1:length(P_j) % for i = 1 to N_j
                                inn = inn + P_j(i) * abs( hArray{n,j,k}'* W_beam{j}(:,i,k))^2;
                            end
                        end
                        inn = inn + params.B_k * params.n0 + cci_dl(n);
                        gamma_dl(m,n)= P_SBS(n,m) / inn;
                        fc9_mn       = params.gam_dl_th *inn - P_SBS(n,m);
                        flag_fc9mn   = fc9_mn >0;
                        penal_fc9_mn = sum(sum(params.nu.* flag_fc9mn.*(fc9_mn.^2)));
                        penal_fc9    = penal_fc9 + penal_fc9_mn;

                        R_dl(m,n,k) = params.B_k .* log2(1+gamma_dl(m,n));
                        % if SBS m not occupy subc k, R_dl(m,:,k) = 0
                    end
                end

            end

        end

        % % calculate threshhold function
        %     R_temp = R_dl > params.R_th; % M_dl x N_dl x K
        %     R_count = sum(sum(R_temp, 3), 2); % M_dl x 1
        % % calculate number of UEs in cells
        %     num_UEs_dl = sum(UE_BS_dl,1)';  % M_dl x 1
        %
        %     o_temp = zeros(M_dl,1);

        % for jj = 1:M_dl
        %     if num_UEs_dl(jj)>0
        %         o_temp(jj) = R_count(jj)/num_UEs_dl(jj);
        %     end
        % end
        % o = sum(o_temp) - penal_fc9 - pnal_fc8;
        o = sum(sum(sum(R_dl))) / 1e8  - penal_fc9 - pnal_fc8;
    end

    function o = FBWOA_ALCA(X)
    	% X == (N_ul + M_dl) x K matrix = asociation matrix
        X_dl = X((N_ul+1):end, :);
        % M_dl x K
        % channel allocation matrix for DL

    	% constraint dealing
        % Constraint of UL
        pnal_fc2 = 0;

        % Constraint of DL
        A_m   = sum(X_dl,2); % == M_dl x 1 matrix
        tol_UE_dl = sum(UE_BS_dl,1)'; % M_dl x 1 == (each element) number of UEs in that cell
        have_UE = tol_UE_dl>0;  % M_dl x 1
        % flag to check that there is UE(s) in the cell (flag for not-null cell)
        fc_dl = abs(A_m -1); % == M_dl x 1 matrix
        flag_fc_dl = fc_dl>0;

        pnal_fc_dl = sum(params.nu .*flag_fc_dl .*(fc_dl.^2) .*have_UE);

     	% objective function
        V_AF= 0;

    	beta_sum = 0; % (beta_t + beta_e)

    	o = beta_sum - V_AF - pnal_fc_dl - pnal_fc2; % obj function has other subtrahends (values of WOA and U_dl obj functions)
        % but we perform it in line 104 code BWOA.m

    end


%% All Remote Joint Optimization Algorithm
    function o = FBWOA_ARJOA(X)
        % X == (N_ul + M_dl) x K matrix = asociation matrix
        X_ul = X(1:N_ul, :);
        % == N_ul x K
        % channel allocation matrix for UL
        X_dl = X((N_ul+1):end, :);
        % M_dl x K
        % channel allocation matrix for DL

    	% constraint dealing
        % Constraint of UL : all UEs have to offload
        A_n  = sum(X_ul,2); % == N_ul x 1 matrix
        fc2  = A_n -1;		% == constraint 2
        flag_fc2 = fc2 ~=0;
        pnal_fc2 = sum(params.nu .* flag_fc2.*(fc2.^2));

        % Constraint of DL
        A_m   = sum(X_dl,2); % == M_dl x 1 matrix
        tol_UE_dl = sum(UE_BS_dl,1)'; % M_dl x 1 == (each element) number of UEs in that cell
        have_UE = tol_UE_dl>0; % flag to check that there is UE(s) in the cell
        fc_dl = abs(A_m -1); % == M_dl x 1 matrix
        flag_fc_dl = fc_dl>0;
        pnal_fc_dl = sum(params.nu .*flag_fc_dl .*(fc_dl.^2) .*have_UE);


     	% objective function
        % get from GRA (33)
        tmp = sqrt(params.beta(:, 1) .* var.f_l); % numerator in (22)
        % == N_ul x 1 matrix
        % UE_BS association in UL
        A_ul = UE_BS(1:N_ul, 1:M_ul);  % == N_ul x M_ul
        % associations between UL UEs and UL SBSs

        A_nm  = repmat(A_ul, 1, 1, noSubcs);    % == N_ul x M_ul x K
        A_nk_temp = reshape(X_ul, N_ul, 1, noSubcs); % == N_ul x 1    x K
        A_nk  = repmat(A_nk_temp, 1, M_ul); % == == N_ul x M_ul x K
        A_nmk = A_nm.* A_nk; % == N_ul x M_ul x K
        % UE-SBS-subc associations

        A  	= sum(A_nmk,3).*tmp;     % N_ul x M_ul matrix
        A1 	= (sum(A,1)).^2;  		 % 1 x M_ul matrix
        V_AF= sum(A1./ params.f0);
        A2  = sum(sum(A_nmk,3),2);

    	beta_sum = sum(sum(params.beta, 2).*A2); % (beta_t + beta_e)

    	o = beta_sum - V_AF - pnal_fc_dl - pnal_fc2; % obj function has other subtrahends (values of WOA and U_dl obj functions)
        % but we perform it in line 104 code BWOA.m
    end

%% Independent Offloading Joint Optimization Algorithm
    function o = FBWOA_IOJOA(X) % A was determined
        % X == (N_ul + M_dl) x K matrix = asociation matrix
        % var.Adet == N_ul x 1
        X_ul = X(1:N_ul, :);
        % == N_ul x K
        % channel allocation matrix for UL
        X_dl = X((N_ul+1):end, :);
        % M_dl x K
        % channel allocation matrix for DL

    	% constraint dealing
        % Constraint of UL: constraint with Adet
        A_n  = sum(X_ul,2); % == N_ul x 1 matrix
        fc2  = A_n - var.Adet;		% == constraint
        flag_fc2 = fc2>0;

        pnal_fc2 = sum(params.nu .* flag_fc2.*(fc2.^2));

        % Constraint of DL
        A_m   = sum(X_dl,2); % == M_dl x 1 matrix
        tol_UE_dl = sum(UE_BS_dl,1)'; % M_dl x 1 == (each element) number of UEs in that cell
        have_UE = tol_UE_dl>0; % flag to check that there is UE(s) in the cell
        fc_dl = abs(A_m -1); % == M_dl x 1 matrix
        flag_fc_dl = fc_dl>0;
        pnal_fc_dl = sum(params.nu .*flag_fc_dl .*(fc_dl.^2) .*have_UE);

     	% objective function
        % get from GRA (33)
        tmp = sqrt(params.beta(:, 1) .* var.f_l); % numerator in (22)
        % == N_ul x 1 matrix
        % UE_BS association in UL
        A_ul = UE_BS(1:N_ul, 1:M_ul);  % == N_ul x M_ul
        % associations between UL UEs and UL SBSs

        A_nm  = repmat(A_ul, 1, 1, noSubcs);    % == N_ul x M_ul x K
        A_nk_temp = reshape(X_ul, N_ul, 1, noSubcs); % == N_ul x 1    x K
        A_nk  = repmat(A_nk_temp, 1, M_ul); % == == N_ul x M_ul x K
        A_nmk = A_nm.* A_nk; % == N_ul x M_ul x K
        % UE-SBS-subc associations

        A  	= sum(A_nmk,3).*tmp;     % N_ul x M_ul matrix
        A1 	= (sum(A,1)).^2;  		 % 1 x M_ul matrix
        V_AF= sum(A1./ params.f0);
        A2  = sum(sum(A_nmk,3),2);

    	beta_sum = sum(sum(params.beta, 2).*A2); % (beta_t + beta_e)

    	o = beta_sum - V_AF - pnal_fc_dl - pnal_fc2; % obj function has other subtrahends (values of WOA and U_dl obj functions)
        % but we perform it in line 104 code BWOA.m
    end

%% OFDMA
    function o = FBWOA_OFDMA(X)
        % X == (N_ul + M_dl) x K matrix = asociation matrix
        X_ul = X(1:N_ul, :);
        % == N_ul x K
        % channel allocation matrix for UL
        X_dl = X((N_ul+1):end, :);
        % M_dl x K
        % channel allocation matrix for DL

    	% constraint dealing
        % Constraint of UL
        A_n  = sum(X_ul,2); % == N_ul x 1 matrix
        fc2  = A_n -1;		% == constraint 2
        flag_fc2 = fc2>0;

        pnal_fc2 = sum(params.nu .* flag_fc2.*(fc2.^2));

        % Constraint of DL
        A_m   = sum(X_dl,2); % == M_dl x 1 matrix
        tol_UE_dl = sum(UE_BS_dl,1)'; % M_dl x 1 == (each element) number of UEs in that cell
        have_UE = tol_UE_dl>0; % flag to check that there is UE(s) in the cell
        fc_dl = abs(A_m -1); % == M_dl x 1 matrix
        flag_fc_dl = fc_dl>0;
        pnal_fc_dl = sum(params.nu .*flag_fc_dl .*(fc_dl.^2) .*have_UE);

        % Constraint of OFDMA
        A_k 	  = sum(X,2); % == 1 x K
        fc_ofdm   = A_k-1;
        flag_ofdm = A_k>1; % == 1 x K

        pnal_ofdm = sum(params.nu .*flag_ofdm.*(fc_ofdm.^2));

     	% objective function
        % get from GRA (33)
        tmp = sqrt(params.beta(:, 1) .* var.f_l); % numerator in (22)
        % == N_ul x 1 matrix
        % UE_BS association in UL
        A_ul = UE_BS(1:N_ul, 1:M_ul);  % == N_ul x M_ul
        % associations between UL UEs and UL SBSs

        A_nm  = repmat(A_ul, 1, 1, noSubcs);    % == N_ul x M_ul x K
        A_nk_temp = reshape(X_ul, N_ul, 1, noSubcs); % == N_ul x 1    x K
        A_nk  = repmat(A_nk_temp, 1, M_ul); % == == N_ul x M_ul x K
        A_nmk = A_nm.* A_nk; % == N_ul x M_ul x K
        % UE-SBS-subc associations

        A  	= sum(A_nmk,3).*tmp;     % N_ul x M_ul matrix
        A1 	= (sum(A,1)).^2;  		 % 1 x M_ul matrix
        V_AF= sum(A1./ params.f0);
        A2  = sum(sum(A_nmk,3),2);

    	beta_sum = sum(sum(params.beta, 2).*A2); % (beta_t + beta_e)

    	o = beta_sum - V_AF - pnal_fc_dl - pnal_fc2 - pnal_ofdm; % obj function has other subtrahends (values of WOA and U_dl obj functions)
        % but we perform it in line 104 code BWOA.m
    end
end