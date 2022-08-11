%5/28/19 Script developed to split the figures into more logical groups.
%State variables will still be plotted as Figs 1&2 on a separate script.
%This particular script will only plot the comparison of tracking error,
%cost, etc. for the different treatments

%% We the People, in Order to form a perfect Union,...
clc;
clearvars;
close all;
disp(['Time at start: ',datestr(now)])

%% Assign moth variables IF NECESSARY
% disp('Assign moth parameters')
L1 = 0.908;
L2 = 1.7475;

mothL = 2*L1 + 2*L2; %in cm
        %The body length of the organism in cm - THIS WILL CHANGE
        %FOR DAUBER (2.31), BEE (1.153), MOTH (5.311).

LSF_1 = 1;
          
%% List the directory stuff
listOf_fa_LSF_1_LE_0_Files = dir('../SimData_MPC/fa_Sims/LSF_1/Winstore_*_fa_LSF_1_LE_0.mat'); 
listOf_fs_LSF_1_LE_0_Files = dir('../SimData_MPC/fs_Sims/LSF_1/Winstore_*_fs_LSF_1_LE_0.mat'); 
listOf_ua_LSF_1_LE_0_Files = dir('../SimData_MPC/ua_Sims/LSF_1/Winstore_*_ua_LSF_1_LE_0.mat'); 
listOf_us_LSF_1_LE_0_Files = dir('../SimData_MPC/us_Sims/LSF_1/Winstore_*_us_LSF_1_LE_0.mat');  
listOf_uaW_LSF_1_LE_0_Files = dir('../SimData_MPC/ua_SimsW/LSF_1/Winstore_*_ua_LSF_1_LE_0.mat'); 
listOf_usW_LSF_1_LE_0_Files = dir('../SimData_MPC/us_SimsW/LSF_1/Winstore_*_us_LSF_1_LE_0.mat'); 

%the asterisk is a wildcard
%The dir function returns a "listing" of an M x 1 "structure." The 
%structure has five fields in this case listing: name, date, byte, isdir, 
%datenum.
%I used the wildcard because I know the number of text files will
%definitely increase as we gather more data.
%For more information enter   help dir   into MATLAB mainframe
disp('Load time vector')
load('../SimData_MPC/Tstore_MPC_hws_sp.mat')

%% Assign the imported numbers to internal variables
disp('Assign internal variables')

hws = 500; 
timestep = numel(Tstore(1:(end-100)))/hws;
% tderiv = Tstore(2);
% t_end = Tstore((end-100+1));

NumOf_fa_LSF_1_LE_0_Files = numel(listOf_fa_LSF_1_LE_0_Files);
NumOf_fs_LSF_1_LE_0_Files = numel(listOf_fs_LSF_1_LE_0_Files);
NumOf_ua_LSF_1_LE_0_Files = numel(listOf_ua_LSF_1_LE_0_Files);
NumOf_us_LSF_1_LE_0_Files = numel(listOf_us_LSF_1_LE_0_Files);
NumOf_uaW_LSF_1_LE_0_Files = numel(listOf_uaW_LSF_1_LE_0_Files);
NumOf_usW_LSF_1_LE_0_Files = numel(listOf_usW_LSF_1_LE_0_Files);

%% Import all the variables
disp('Importing curated data')

%imported cost
load('../CuratedData_MPC/fa/LSF_1/mean_impcost_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/mean_impcost_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/mean_impcost_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/mean_impcost_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/mean_impcost_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/mean_impcost_usW_LSF_1_LE_0.mat')

load('../CuratedData_MPC/fa/LSF_1/std_impcost_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/std_impcost_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/std_impcost_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/std_impcost_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/std_impcost_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/std_impcost_usW_LSF_1_LE_0.mat')

%tepfr -- tracking error per full run
load('../CuratedData_MPC/fa/LSF_1/mean_tepfr_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/mean_tepfr_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/mean_tepfr_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/mean_tepfr_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/mean_tepfr_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/mean_tepfr_usW_LSF_1_LE_0.mat')

%Mean of cost per full run -- for stats purposes
load('../CuratedData_MPC/fa/LSF_1/mean_cost_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/mean_cost_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/mean_cost_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/mean_cost_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/mean_cost_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/mean_cost_usW_LSF_1_LE_0.mat')

%Tracking error w.r.t. time
load('../CuratedData_MPC/fa/LSF_1/trackingerror_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/trackingerror_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/trackingerror_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/trackingerror_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/trackingerror_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/trackingerror_usW_LSF_1_LE_0.mat')

disp('Curated data is loaded')

%% Find the min & max values for...

%tepfr
tepfrMax_vec = [max(max(mean_tepfr_fa_LSF_1_LE_0)),...
    max(max(mean_tepfr_fs_LSF_1_LE_0)),...
    max(max(mean_tepfr_ua_LSF_1_LE_0)),...
    max(max(mean_tepfr_us_LSF_1_LE_0)),...
    max(max(mean_tepfr_uaW_LSF_1_LE_0)),...
    max(max(mean_tepfr_usW_LSF_1_LE_0))];
tepfr_max = max(tepfrMax_vec);

tepfrMin_vec = [min(min(mean_tepfr_fa_LSF_1_LE_0)),...
    min(min(mean_tepfr_fs_LSF_1_LE_0)),...
    min(min(mean_tepfr_ua_LSF_1_LE_0)),...
    min(min(mean_tepfr_us_LSF_1_LE_0)),...
    min(min(mean_tepfr_uaW_LSF_1_LE_0)),...
    min(min(mean_tepfr_usW_LSF_1_LE_0))];
tepfr_min = min(tepfrMin_vec);

%cost
costMax_vec = [max(max(mean_cost_fa_LSF_1_LE_0)),...
    max(max(mean_cost_fs_LSF_1_LE_0)),...
    max(max(mean_cost_ua_LSF_1_LE_0)),...
    max(max(mean_cost_us_LSF_1_LE_0)),...
    max(max(mean_cost_uaW_LSF_1_LE_0)),...
    max(max(mean_cost_usW_LSF_1_LE_0))];
cost_max = max(costMax_vec);

costMin_vec = [min(min(mean_cost_fa_LSF_1_LE_0)),...
    min(min(mean_cost_fs_LSF_1_LE_0)),...
    min(min(mean_cost_ua_LSF_1_LE_0)),...
    min(min(mean_cost_us_LSF_1_LE_0)),...
    min(min(mean_cost_uaW_LSF_1_LE_0)),...
    min(min(mean_cost_usW_LSF_1_LE_0))];
cost_min = min(costMin_vec);

%% Tracking error violin plots unscaled

fig3 = figure(3);
%Tracking error - fa LSF_1_LE_0
subplot(1,6,1)
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOf_fa_LSF_1_LE_0_Files),')'])
if NumOf_fa_LSF_1_LE_0_Files > 0
    distributionPlot(mean_tepfr_fa_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, tepfr_min, tepfr_max])
    % set(gca,'ytick',0:0.5:1);
%     axis([0, 2, 0.125, 0.2])
end
ylabel('Mean of tracking error (cm)')

%Tracking error - fs LSF_1_LE_0
subplot(1,6,2)
title({'fs';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOf_fs_LSF_1_LE_0_Files),')'])
if NumOf_fs_LSF_1_LE_0_Files > 0
    distributionPlot(mean_tepfr_fs_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, tepfr_min, tepfr_max])
    % set(gca,'ytick',0:0.5:1);
%     axis([0, 2, 0.125, 0.2])
end

%Tracking error - ua LSF_1_LE_0
subplot(1,6,3)
title({'ua';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOf_ua_LSF_1_LE_0_Files),')'])
if NumOf_ua_LSF_1_LE_0_Files > 0
    distributionPlot(mean_tepfr_ua_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, tepfr_min, tepfr_max])
    % set(gca,'ytick',0:0.5:1);
%     axis([0, 2, 0.125, 0.2])
end

%Tracking error - us LSF_1_LE_0
subplot(1,6,4)
title({'us';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOf_us_LSF_1_LE_0_Files),')'])
if NumOf_us_LSF_1_LE_0_Files > 0
    distributionPlot(mean_tepfr_us_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, tepfr_min, tepfr_max])
    % set(gca,'ytick',0:10:80);
%     axis([0, 2, 0.125, 0.2])
end

%Tracking error - uaW LSF_1_LE_0
subplot(1,6,5)
title({'uaW';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOf_uaW_LSF_1_LE_0_Files),')'])
if NumOf_uaW_LSF_1_LE_0_Files > 0
    distributionPlot(mean_tepfr_uaW_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, tepfr_min, tepfr_max])
    % set(gca,'ytick',0:0.5:1);
%     axis([0, 2, 0.125, 0.2])
end

%Tracking error - usW LSF_1_LE_0
subplot(1,6,6)
title({'usW';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOf_usW_LSF_1_LE_0_Files),')'])
if NumOf_usW_LSF_1_LE_0_Files > 0
    distributionPlot(mean_tepfr_usW_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, tepfr_min, tepfr_max])
    % set(gca,'ytick',0:10:80);
%     axis([0, 2, 0.125, 0.2])
end

set(fig3, 'Position', [120,120,850,525])

saveas(gcf,'../Figures_MPC/compTmt/LSF_1/Fig3_compTmt_TrackErrorViolin_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/compTmt/LSF_1/Fig3_compTmt_TrackErrorVoilin_LSF_1_v12h.jpg');

disp('Figure 3 done') 

%% Tracking error normalized to body length

fig4 = figure(4);
%Tracking error - fa LSF_1_LE_0
subplot(1,6,1)
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOf_fa_LSF_1_LE_0_Files),')'])
if NumOf_fa_LSF_1_LE_0_Files > 0
    distributionPlot(mean_tepfr_fa_LSF_1_LE_0./(mothL*LSF_1),'addSpread',true)
    axis([0, 2, tepfr_min/(mothL*LSF_1), tepfr_max/(mothL*LSF_1)])
    % set(gca,'ytick',0:0.05:0.2);
end
ylabel({'Mean of tracking error';'normalized to body length'})

%Tracking error - fs LSF_1_LE_0
subplot(1,6,2)
title({'fs';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOf_fs_LSF_1_LE_0_Files),')'])
if NumOf_fs_LSF_1_LE_0_Files > 0
    distributionPlot(mean_tepfr_fs_LSF_1_LE_0./(mothL*LSF_1),'addSpread',true)
    axis([0, 2, tepfr_min/(mothL*LSF_1), tepfr_max/(mothL*LSF_1)])
    % set(gca,'ytick',0:0.05:0.2);
end

%Tracking error - ua LSF_1_LE_0
subplot(1,6,3)
title({'ua';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOf_ua_LSF_1_LE_0_Files),')'])
if NumOf_ua_LSF_1_LE_0_Files > 0
    distributionPlot(mean_tepfr_ua_LSF_1_LE_0./(mothL*LSF_1),'addSpread',true)
    axis([0, 2, tepfr_min/(mothL*LSF_1), tepfr_max/(mothL*LSF_1)])
    % set(gca,'ytick',0:0.05:0.2);
end

%Tracking error - us LSF_1_LE_0
subplot(1,6,4)
title({'us';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOf_us_LSF_1_LE_0_Files),')'])
if NumOf_us_LSF_1_LE_0_Files > 0
    distributionPlot(mean_tepfr_us_LSF_1_LE_0./(mothL*LSF_1),'addSpread',true)
    axis([0, 2, tepfr_min/(mothL*LSF_1), tepfr_max/(mothL*LSF_1)])
    % set(gca,'ytick',0:2:16);
end

%Tracking error - uaW LSF_1_LE_0
subplot(1,6,5)
title({'uaW';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOf_uaW_LSF_1_LE_0_Files),')'])
if NumOf_uaW_LSF_1_LE_0_Files > 0
    distributionPlot(mean_tepfr_uaW_LSF_1_LE_0./(mothL*LSF_1),'addSpread',true)
    axis([0, 2, tepfr_min/(mothL*LSF_1), tepfr_max/(mothL*LSF_1)])
    % set(gca,'ytick',0:0.05:0.2);
end

%Tracking error - usW LSF_1_LE_0
subplot(1,6,6)
title({'usW';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOf_usW_LSF_1_LE_0_Files),')'])
if NumOf_usW_LSF_1_LE_0_Files > 0
    distributionPlot(mean_tepfr_usW_LSF_1_LE_0./(mothL*LSF_1),'addSpread',true)
    axis([0, 2, tepfr_min/(mothL*LSF_1), tepfr_max/(mothL*LSF_1)])
    % set(gca,'ytick',0:2:16);
end

set(fig4, 'Position', [120,120,850,525])

saveas(gcf,'../Figures_MPC/compTmt/LSF_1/Fig4_compTmt_TrackErrorViolin_lengthscaled_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/compTmt/LSF_1/Fig4_compTmt_TrackErrorVoilin_lengthscaled_LSF_1_v12h.jpg');

disp('Figure 4 done')

%% Cost function violin plots

fig5 = figure(5);
%Cost - fa LSF_1_LE_0
subplot(1,6,1)
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOf_fa_LSF_1_LE_0_Files),')'])
if NumOf_fa_LSF_1_LE_0_Files > 0
    distributionPlot(mean_cost_fa_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, cost_min, cost_max])
    % set(gca,'ytick',0:2e6:13e6);
end
ylabel('Mean of cost function value')

%Cost - fs LSF_1_LE_0
subplot(1,6,2)
title({'fs';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOf_fs_LSF_1_LE_0_Files),')'])
if NumOf_fs_LSF_1_LE_0_Files > 0
    distributionPlot(mean_cost_fs_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, cost_min, cost_max])
    % set(gca,'ytick',0:2e6:13e6);
end

%Cost - ua LSF_1_LE_0
subplot(1,6,3)
title({'ua';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOf_ua_LSF_1_LE_0_Files),')'])
if NumOf_ua_LSF_1_LE_0_Files > 0
    distributionPlot(mean_cost_ua_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, cost_min, cost_max])
    % set(gca,'ytick',0:2e6:13e6);
end

%Cost - us LSF_1_LE_0
subplot(1,6,4)
xlabel(['(n=',num2str(NumOf_us_LSF_1_LE_0_Files),')'])
title({'us';['LSF: ',num2str(LSF_1)];'LE: 0'})
if NumOf_us_LSF_1_LE_0_Files > 0
    distributionPlot(mean_cost_us_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, cost_min, cost_max])
    % set(gca,'ytick',0:2e11:2.5e11);
end

%Cost - uaW LSF_1_LE_0
subplot(1,6,5)
title({'uaW';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOf_uaW_LSF_1_LE_0_Files),')'])
if NumOf_uaW_LSF_1_LE_0_Files > 0
    distributionPlot(mean_cost_uaW_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, cost_min, cost_max])
    % set(gca,'ytick',0:2e6:13e6);
end

%Cost - usW LSF_1_LE_0
subplot(1,6,6)
xlabel(['(n=',num2str(NumOf_usW_LSF_1_LE_0_Files),')'])
title({'usW';['LSF: ',num2str(LSF_1)];'LE: 0'})
if NumOf_usW_LSF_1_LE_0_Files > 0
    distributionPlot(mean_cost_usW_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, cost_min, cost_max])
    % set(gca,'ytick',0:2e11:2.5e11);
end

saveas(gcf,'../Figures_MPC/compTmt/LSF_1/Fig5_compTmt_CostViolin_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/compTmt/LSF_1/Fig5_compTmt_CostViolin_LSF_1_v12h.jpg');

set(fig5, 'Position', [120,120,850,525])

disp('Figure 5 done') %p-value = 5.323e-11

%% log plot of cost

figure(6);
semilogy(Tstore(1:(timestep*0.25):(end-100)),mean_impcost_fa_LSF_1_LE_0,'Linewidth',2)
hold on;
semilogy(Tstore(1:(timestep*0.25):(end-100)),mean_impcost_fs_LSF_1_LE_0,'Linewidth',2)
hold on;
semilogy(Tstore(1:(timestep*0.25):(end-100)),mean_impcost_ua_LSF_1_LE_0,'Linewidth',2)
hold on;
semilogy(Tstore(1:(timestep*0.25):(end-100)),mean_impcost_us_LSF_1_LE_0,'Linewidth',2)
hold on;
semilogy(Tstore(1:(timestep*0.25):(end-100)),mean_impcost_uaW_LSF_1_LE_0,'Linewidth',2)
hold on;
semilogy(Tstore(1:(timestep*0.25):(end-100)),mean_impcost_usW_LSF_1_LE_0,'Linewidth',2)
hold on;
xlabel('Time (seconds)'); ylabel('Cost function value')
title('Comparing treatments')
legend(['fa, ','LSF: ',num2str(LSF_1),' LE: 0 - ', num2str(NumOf_fa_LSF_1_LE_0_Files)],...
    ['fs, ','LSF: ',num2str(LSF_1),' LE: 0 - ', num2str(NumOf_fs_LSF_1_LE_0_Files)],...
    ['ua, ','LSF: ',num2str(LSF_1),' LE: 0 - ', num2str(NumOf_ua_LSF_1_LE_0_Files)],...
    ['us, ','LSF: ',num2str(LSF_1),' LE: 0 - ', num2str(NumOf_us_LSF_1_LE_0_Files)],...
    ['uaW, ','LSF: ',num2str(LSF_1),' LE: 0 - ', num2str(NumOf_uaW_LSF_1_LE_0_Files)],...
    ['usW, ','LSF: ',num2str(LSF_1),' LE: 0 - ', num2str(NumOf_usW_LSF_1_LE_0_Files)])

saveas(gcf,'../Figures_MPC/compTmt/LSF_1/Fig6_compTmt_LogCost_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/compTmt/LSF_1/Fig6_compTmt_LogCost_LSF_1_v12h.jpg');

disp('Figure 6 done') 

%% Tracking error figure for SACNAS 2019

% figSACNAS = figure(7);
% %Tracking error - fs LSF_1_LE_0
% subplot(1,6,1)
% title({'fs'})
% xlabel(['(n=',num2str(NumOf_fs_LSF_1_LE_0_Files),')'])
% if NumOf_fs_LSF_1_LE_0_Files > 0
%     distributionPlot(mean_tepfr_fs_LSF_1_LE_0,'addSpread',true)
%     axis([0, 2, 0, 2])
%     % set(gca,'ytick',0:0.5:1);
% %     axis([0, 2, 0.125, 0.2])
% end
% ylabel('Mean of tracking error (cm)')
% 
% %Tracking error - fa LSF_1_LE_0
% subplot(1,6,2)
% 
% title({'fa'})
% xlabel(['(n=',num2str(NumOf_fa_LSF_1_LE_0_Files),')'])
% if NumOf_fa_LSF_1_LE_0_Files > 0
%     distributionPlot(mean_tepfr_fa_LSF_1_LE_0,'addSpread',true)
%     axis([0, 2, 0, 2])
%     % set(gca,'ytick',0:0.5:1);
% %     axis([0, 2, 0.125, 0.2])
% end
% 
% %Tracking error - us LSF_1_LE_0
% subplot(1,6,3)
% title({'us'})
% xlabel(['(n=',num2str(NumOf_us_LSF_1_LE_0_Files),')'])
% if NumOf_us_LSF_1_LE_0_Files > 0
%     distributionPlot(mean_tepfr_us_LSF_1_LE_0,'addSpread',true)
%     axis([0, 2, 0, 2])
%     % set(gca,'ytick',0:10:80);
% %     axis([0, 2, 0.125, 0.2])
% end
% 
% %Tracking error - ua LSF_1_LE_0
% subplot(1,6,4)
% title({'ua'})
% xlabel(['(n=',num2str(NumOf_ua_LSF_1_LE_0_Files),')'])
% if NumOf_ua_LSF_1_LE_0_Files > 0
%     distributionPlot(mean_tepfr_ua_LSF_1_LE_0,'addSpread',true)
%     axis([0, 2, 0, 2])
%     % set(gca,'ytick',0:0.5:1);
% %     axis([0, 2, 0.125, 0.2])
% end
% 
% set(figSACNAS, 'Position', [120,120,850,525])
% 
% saveas(gcf,'../Figures_MPC/compTmt/LSF_1/FigSACNAS2019_compTmt_TrackErrorViolin_LSF_1_v12h.fig');
% saveas(gcf,'../Figures_MPC/compTmt/LSF_1/FigSACNAS2019_compTmt_TrackErrorVoilin_LSF_1_v12h.jpg');

% disp('Figure for SACNAS done') 

%% Final window
disp('Completely done, son')
disp(['Time at end: ',datestr(now)])