%----------------------------------------------------------------------------------------
%%% Initialze positions of whales (subchannel assignment in the case of A --> posi_a
%%%								  UEs' transmiting power in the case of P --> posi_p )
%----------------------------------------------------------------------------------------

% Output:
% posi_a = (N_ul + M_dl) x K x noSA  matrix	== positions of whales in 3-D searching for the best associations

function [posi_a] = Initialization2(functionname, noSubcs, UEs, BS, UE_BS, noSearchAgents, Adet)
% functionname 	== 'string'
% noSubcs = K 			   == number of subchannel
% UE_BS    == N_active x M matrix % matrix of relation of UEs and SBSs
%                                 % contains just 0's and 1's
% UEs == 1x1 struct
%       UEs.total    == 1 x 2 matrix == [N_ul N_dl, N]
% BS  == 1x1 struct
%       BS.total = [M_ul M_dl M]
% Adet == N_ul x 1

noUE_ul = UEs.total(1);
noBS_ul = BS.total(1);
noBS_dl = BS.total(2);

% Adet	  = N x 1 matrix

posi_a = zeros(noUE_ul + noBS_dl, noSubcs, noSearchAgents); % (N_ul + M_dl) x K x noSA matrix

% setup posi_a_dl
posi_a_dl = zeros(noBS_dl, noSubcs, noSearchAgents);
SBS_busy  = sum(UE_BS,1)>0; % SBSs that contain active UEs
% M x 1

% for i = 1:noSearchAgents
%     k = 0;
%     for m = 1:noBS_dl
%         if (SBS_busy(1,noBS_ul + m)) % just considering DL SBSs
%             posi_a_dl(m,k+1,i) = 1;     % allocate 1 subchannel for 1 DL SBS
%             k = rem(k+1, noSubcs);      % if no. DL SBSs < no. subchannels, return to subchannel 1
%         end
%     end
% end

% setup posi_a_ul
posi_a_ul = zeros(noUE_ul, noSubcs, noSearchAgents);
switch functionname
    case 'ALCA'
        k_count = 0;
        for i = 1:noSearchAgents
            k_count = k_count+1;
            k = k_count;
            for m = 1:noBS_dl
                if (SBS_busy(1,noBS_ul + m)) % just considering DL SBSs
                    k = rem(k+1, noSubcs);      % if no. DL SBSs < no. subchannels, return to subchannel 1
                    posi_a_dl(m,k+1,i) = 1;     % allocate 1 subchannel for 1 DL SBS
                end
            end
        end
        posi_a = [posi_a_ul; posi_a_dl];
    case {'WOA_SIC_MEC', "IWOA_SIC_MEC", 'PSO_SIC_MEC'}
        k_count = 0;
        for i = 1:noSearchAgents
            k_count = k_count+1;
            k = k_count;
            for m = 1:noBS_dl
                if (SBS_busy(1,noBS_ul + m)) % just considering DL SBSs
                    k = rem(k+1, noSubcs);      % if no. DL SBSs < no. subchannels, return to subchannel 1
                    posi_a_dl(m,k+1,i) = 1;     % allocate 1 subchannel for 1 DL SBS
                end
            end
            for n = 1:noUE_ul
                if rand > 0.5  % probability of offloading is 50%
                    rand_idx = rem(floor((noSubcs-k)*rand + 1) +k, noSubcs) +1; % randomly allocate 1 subchannel to the UL UE
                    % prior to the subchannels that have not been occupied by DL SBSs
                    posi_a_ul(n, rand_idx, i) = 1;
                end
            end
        end
        posi_a = [posi_a_ul; posi_a_dl];
        %            posi_a = zeros(noUsers, noBSs, noSubcs, noSearchAgents);


        % All Remote Joint Optimization Algorithm
    case 'ARJOA'
        k_count = 0;
        for i = 1:noSearchAgents
            k_count = k_count+1;
            k = k_count;
            for m = 1:noBS_dl
                if (SBS_busy(1,noBS_ul + m)) % just considering DL SBSs
                    k = rem(k+1, noSubcs);      % if no. DL SBSs < no. subchannels, return to subchannel 1
                    posi_a_dl(m,k+1,i) = 1;     % allocate 1 subchannel for 1 DL SBS
                end
            end
            for n = 1:noUE_ul
                rand_idx = rem(floor((noSubcs-k)*rand + 1) +k, noSubcs) +1; % randomly allocate 1 subchannel to the UL UE
                % prior to the subchannels that have not been occupied by DL SBSs
                posi_a_ul(n, rand_idx, i) = 1;
            end
        end
        posi_a = [posi_a_ul; posi_a_dl];
        % Independent Offloading Joint Optimization Algorithm
    case 'IOJOA'
        k_count = 0;
        for i = 1:noSearchAgents
            k_count = k_count+1;
            k = k_count;
            for m = 1:noBS_dl
                if (SBS_busy(1,noBS_ul + m)) % just considering DL SBSs
                    k = rem(k+1, noSubcs);      % if no. DL SBSs < no. subchannels, return to subchannel 1
                    posi_a_dl(m,k+1,i) = 1;     % allocate 1 subchannel for 1 DL SBS
                end
            end
            for n = 1:noUE_ul
                if Adet(n) == 1
                    rand_idx = rem(floor((noSubcs-k)*rand + 1) +k, noSubcs) +1; % randomly allocate 1 subchannel to the UL UE
                    % prior to the subchannels that have not been occupied by DL SBSs
                    posi_a_ul(n, rand_idx, i) = 1;
                end
            end
        end
        posi_a = [posi_a_ul; posi_a_dl];
        % OFDMA
    case 'OFDMA'
        k_count = 0;
        for i = 1:noSearchAgents
            k_count = k_count+1;
            k = k_count;
            for m = 1:noBS_dl
                if (SBS_busy(1,noBS_ul + m)) % just considering DL SBSs
                    k = rem(k+1, noSubcs);      % if no. DL SBSs < no. subchannels, return to subchannel 1
                    posi_a_dl(m,k+1,i) = 1;     % allocate 1 subchannel for 1 DL SBS
                end
            end

            sub_remain = noSubcs - sum(SBS_busy(noBS_ul + 1 : end));
            % == number of sabcarriers remaining after allocated to DL SBSs
            if sub_remain > 0
                UE_BS_ul = UE_BS(1:noUE_ul, 1:noBS_ul);
                [xx,~]  = find(UE_BS_ul==1); % find UE-BS associations
                temp    = min(length(xx), sub_remain);
                subc_   = rem(k + randperm(sub_remain,temp), noSubcs)+1; % random occupied subcarriers
                % == 1 x temp vector
                u_      = randperm(length(xx),temp); % random associations
                % UE xx(u_(j)) offloads
                % to BS yy(u_(j))
                for j = 1:length(subc_)
                    posi_a_ul(xx(u_(j)), subc_(j), i) = 1;
                end
            end
        end
        posi_a = [posi_a_ul; posi_a_dl];
end
end


