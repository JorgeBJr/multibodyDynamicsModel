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
listOfLSF_1_LE_0_Files = dir('../SimData_MPC/fa_Sims/LSF_1/Winstore_*_fa_LSF_1_LE_0.mat'); 
listOfLSF_1_LE_p2_Files = dir('../SimData_MPC/fa_Sims/LSF_1/Winstore_*_fa_LSF_1_LE_p2.mat'); 
listOfLSF_1_LE_p4_Files = dir('../SimData_MPC/fa_Sims/LSF_1/Winstore_*_fa_LSF_1_LE_p4.mat'); 

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

NumOfLSF_1_LE_0_Files = numel(listOfLSF_1_LE_0_Files);
NumOfLSF_1_LE_p2_Files = numel(listOfLSF_1_LE_p2_Files);
NumOfLSF_1_LE_p4_Files = numel(listOfLSF_1_LE_p4_Files);

%% Import all the variables
disp('Importing curated data')

%imported cost
load('../CuratedData_MPC/fa/LSF_1/mean_impcost_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_impcost_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_impcost_fa_LSF_1_LE_p4.mat')

load('../CuratedData_MPC/fa/LSF_1/std_impcost_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/std_impcost_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/std_impcost_fa_LSF_1_LE_p4.mat')

%tepfr -- tracking error per full run
load('../CuratedData_MPC/fa/LSF_1/mean_tepfr_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_tepfr_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_tepfr_fa_LSF_1_LE_p4.mat')

%Mean of cost per full run -- for stats purposes
load('../CuratedData_MPC/fa/LSF_1/mean_cost_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_cost_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_cost_fa_LSF_1_LE_p4.mat')

%Tracking error w.r.t. time
load('../CuratedData_MPC/fa/LSF_1/trackingerror_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/trackingerror_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/trackingerror_fa_LSF_1_LE_p4.mat')

disp('Curated data is loaded')

%% Find the min & max values for...

%tepfr
tepfrMax_vec = [max(max(mean_tepfr_fa_LSF_1_LE_0)),...
    max(max(mean_tepfr_fa_LSF_1_LE_p2)),...
    max(max(mean_tepfr_fa_LSF_1_LE_p4))];
tepfr_max = max(tepfrMax_vec);

tepfrMin_vec = [min(min(mean_tepfr_fa_LSF_1_LE_0)),...
    min(min(mean_tepfr_fa_LSF_1_LE_p2)),...
    min(min(mean_tepfr_fa_LSF_1_LE_p4))];
tepfr_min = min(tepfrMin_vec);

%cost
costMax_vec = [max(max(mean_cost_fa_LSF_1_LE_0)),...
    max(max(mean_cost_fa_LSF_1_LE_p2)),...
    max(max(mean_cost_fa_LSF_1_LE_p4))];
cost_max = max(costMax_vec);

costMin_vec = [min(min(mean_cost_fa_LSF_1_LE_0)),...
    min(min(mean_cost_fa_LSF_1_LE_p2)),...
    min(min(mean_cost_fa_LSF_1_LE_p4))];
cost_min = min(costMin_vec);

%% Tracking error violin plots unscaled

fig3 = figure(3);
%Tracking error - con LSF_1_LE_0
subplot(1,3,1)
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOfLSF_1_LE_0_Files),')'])
if NumOfLSF_1_LE_0_Files > 0
    distributionPlot(mean_tepfr_fa_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, tepfr_min, tepfr_max])
    % set(gca,'ytick',0:0.5:1);
%     axis([0, 2, 0.125, 0.2])
end
ylabel('Mean of tracking error (cm)')

%Tracking error - con LSF_1_LE_p2
subplot(1,3,2)
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0.2'})
xlabel(['(n=',num2str(NumOfLSF_1_LE_p2_Files),')'])
if NumOfLSF_1_LE_p2_Files > 0
    distributionPlot(mean_tepfr_fa_LSF_1_LE_p2,'addSpread',true)
    axis([0, 2, tepfr_min, tepfr_max])
    % set(gca,'ytick',0:0.5:1);
%     axis([0, 2, 0.125, 0.2])
end

%Tracking error - con LSF_1_LE_p4
subplot(1,3,3)
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0.4'})
xlabel(['(n=',num2str(NumOfLSF_1_LE_p4_Files),')'])
if NumOfLSF_1_LE_p4_Files > 0
    distributionPlot(mean_tepfr_fa_LSF_1_LE_p4,'addSpread',true)
    axis([0, 2, tepfr_min, tepfr_max])
    % set(gca,'ytick',0:0.5:1);
%     axis([0, 2, 0.125, 0.2])
end

set(fig3, 'Position', [120,120,850,525])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig3_fa_TrackErrorViolin_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig3_fa_TrackErrorVoilin_LSF_1_v12h.jpg');

disp('Figure 3 done') 

%% Tracking error scaled

fig4 = figure(4);
%Tracking error - con LSF_1_LE_0
subplot(1,3,1)
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOfLSF_1_LE_0_Files),')'])
if NumOfLSF_1_LE_0_Files > 0
    distributionPlot(mean_tepfr_fa_LSF_1_LE_0./(mothL*LSF_1),'addSpread',true)
    axis([0, 2, tepfr_min/(mothL*LSF_1), tepfr_max/(mothL*LSF_1)])
    % set(gca,'ytick',0:0.05:0.2);
end
ylabel({'Mean of tracking error';'normalized to body length'})

%Tracking error - con LSF_1_LE_p2
subplot(1,3,2)
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0.2'})
xlabel(['(n=',num2str(NumOfLSF_1_LE_p2_Files),')'])
if NumOfLSF_1_LE_p2_Files > 0
    distributionPlot(mean_tepfr_fa_LSF_1_LE_p2./(mothL*LSF_1),'addSpread',true)
    axis([0, 2, tepfr_min/(mothL*LSF_1), tepfr_max/(mothL*LSF_1)])
    % set(gca,'ytick',0:0.05:0.2);
end

%Tracking error - con LSF_1_LE_p4
subplot(1,3,3)
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0.4'})
xlabel(['(n=',num2str(NumOfLSF_1_LE_p4_Files),')'])
if NumOfLSF_1_LE_p4_Files > 0
    distributionPlot(mean_tepfr_fa_LSF_1_LE_p4./(mothL*LSF_1),'addSpread',true)
    axis([0, 2, tepfr_min/(mothL*LSF_1), tepfr_max/(mothL*LSF_1)])
    % set(gca,'ytick',0:0.05:0.2);
end

set(fig4, 'Position', [120,120,850,525])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig4_fa_TrackErrorViolin_lengthscaled_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig4_fa_TrackErrorVoilin_lengthscaled_LSF_1_v12h.jpg');

disp('Figure 4 done')

%% Cost function violin plots

fig6 = figure(5);
%Cost - con LSF_1_LE_0
subplot(1,3,1)
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0'})
xlabel(['(n=',num2str(NumOfLSF_1_LE_0_Files),')'])
if NumOfLSF_1_LE_0_Files > 0
    distributionPlot(mean_cost_fa_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, cost_min, cost_max])
    % set(gca,'ytick',0:2e6:13e6);
end
ylabel('Mean of cost function value')

%Cost - con LSF_1_LE_p2
subplot(1,3,2)
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0.2'})
xlabel(['(n=',num2str(NumOfLSF_1_LE_p2_Files),')'])
if NumOfLSF_1_LE_p2_Files > 0
    distributionPlot(mean_cost_fa_LSF_1_LE_p2,'addSpread',true)
    axis([0, 2, cost_min, cost_max])
    % set(gca,'ytick',0:2e6:13e6);
end

%Cost - con LSF_1_LE_p4
subplot(1,3,3)
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0.4'})
xlabel(['(n=',num2str(NumOfLSF_1_LE_p4_Files),')'])
if NumOfLSF_1_LE_p4_Files > 0
    distributionPlot(mean_cost_fa_LSF_1_LE_p4,'addSpread',true)
    axis([0, 2, cost_min, cost_max])
    % set(gca,'ytick',0:2e6:13e6);
end

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig5_fa_CostViolin_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig5_fa_CostViolin_LSF_1_v12h.jpg');

set(fig6, 'Position', [120,120,850,525])

disp('Figure 5 done') %p-value = 5.323e-11

%% log plot of cost

figure(6);
semilogy(Tstore(1:(timestep*0.25):(end-100)),mean_impcost_fa_LSF_1_LE_0,'k','Linewidth',2)
hold on;
semilogy(Tstore(1:(timestep*0.25):(end-100)),mean_impcost_fa_LSF_1_LE_p2,'Linewidth',2)
hold on;
semilogy(Tstore(1:(timestep*0.25):(end-100)),mean_impcost_fa_LSF_1_LE_p4,'Linewidth',2)
hold on;
xlabel('Time (seconds)'); ylabel('Cost function value')
title('fa')
legend(['LSF: ',num2str(LSF_1),' LE: 0 - ', num2str(NumOfLSF_1_LE_0_Files)],...
    ['LSF: ',num2str(LSF_1),' LE: 0.2 - ', num2str(NumOfLSF_1_LE_p2_Files)],...
    ['LSF: ',num2str(LSF_1),' LE: 0.4 - ', num2str(NumOfLSF_1_LE_p4_Files)])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig6_fa_LogCost_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig6_fa_LogCost_LSF_1_v12h.jpg');

disp('Figure 6 done') 

%% Final window
disp('Completely done, son')
disp(['Time at end: ',datestr(now)])