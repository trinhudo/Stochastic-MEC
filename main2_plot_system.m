% Script to plot the system
% change noUsers and noSBS

noUsers = 13;
M_ul = 2;
M_dl = 2;
flag_plot = 1;

[UE_BS, UEs, BS] = location_voronoi(noUsers, M_ul, M_dl, flag_plot);

xlim([-500 500]);
ylim([-500 500]);

save('pos_BS_UEs.mat', 'UEs', 'BS', "UE_BS")