%11/15/2018 Script developed to plot the tracking error per full run
%(tepfr) and cost function value with respect to Length Size factor (LSF)
%and Length Extension (LE).

%% We the People, in Order to form a perfect Union,...
clc;
clearvars;
close all;
disp(['Time at start: ',datestr(now)])

%% Assign moth variables IF NECESSARY
disp('Assign moth lengths, and LSF values')
L1 = 0.908;
L2 = 1.7475;

%Length Scale Factors
LSF_p5 = 0.5;
LSF_1 = 1;
% LSF_2 = 2;

bl = 5.311; %The body length of the organism in cm - THIS WILL CHANGE
          %FOR DAUBER (2.31), BEE (1.153), MOTH (5.311).

%The spread of values for the figure axes
LSF_spread = 1:2; %This is 1:3 when LSF_2 is included
LE_spread = 1:3;

%Assign the tepfr and cost matrices
mean_tepfr_fa_LSF_1 = nan(40,3);
% mean_tepfr_fa_LSF_2 = nan(40,3);
mean_tepfr_fa_LSF_p5 = nan(40,3);

mean_cost_fa_LSF_1 = nan(40,3);
% mean_cost_fa_LSF_2 = nan(40,3);
mean_cost_fa_LSF_p5 = nan(40,3);
          
%% Import all the variables
disp('Importing curated data')

%LSF_p5
%tepfr -- tracking error per full run
load('../CuratedData_MPC/fa/LSF_p5/mean_tepfr_fa_LSF_p5_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_p5/mean_tepfr_fa_LSF_p5_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_p5/mean_tepfr_fa_LSF_p5_LE_p4.mat')

for i = 1:numel(mean_tepfr_fa_LSF_p5_LE_0)
    mean_tepfr_fa_LSF_p5(i,1) = mean_tepfr_fa_LSF_p5_LE_0(i,1);
end
for i = 1:numel(mean_tepfr_fa_LSF_p5_LE_p2)
    mean_tepfr_fa_LSF_p5(i,2) = mean_tepfr_fa_LSF_p5_LE_p2(i,1);
end
for i = 1:numel(mean_tepfr_fa_LSF_p5_LE_p4)
    mean_tepfr_fa_LSF_p5(i,3) = mean_tepfr_fa_LSF_p5_LE_p4(i,1);
end

clear mean_tepfr_fa_LSF_p5_LE_0
clear mean_tepfr_fa_LSF_p5_LE_p2
clear mean_tepfr_fa_LSF_p5_LE_p4

%Mean of cost per full run -- for stats purposes
load('../CuratedData_MPC/fa/LSF_p5/mean_cost_fa_LSF_p5_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_p5/mean_cost_fa_LSF_p5_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_p5/mean_cost_fa_LSF_p5_LE_p4.mat')

for i = 1:numel(mean_cost_fa_LSF_p5_LE_0)
    mean_cost_fa_LSF_p5(i,1) = mean_cost_fa_LSF_p5_LE_0(i,1);
end
for i = 1:numel(mean_cost_fa_LSF_p5_LE_p2)
    mean_cost_fa_LSF_p5(i,2) = mean_cost_fa_LSF_p5_LE_p2(i,1);
end
for i = 1:numel(mean_cost_fa_LSF_p5_LE_p4)
    mean_cost_fa_LSF_p5(i,3) = mean_cost_fa_LSF_p5_LE_p4(i,1);
end

clear mean_cost_fa_LSF_p5_LE_0
clear mean_cost_fa_LSF_p5_LE_p2
clear mean_cost_fa_LSF_p5_LE_p4

%LSF_1
%tepfr -- tracking error per full run
load('../CuratedData_MPC/fa/LSF_1/mean_tepfr_fa_LSF_1_LE_0.mat');
load('../CuratedData_MPC/fa/LSF_1/mean_tepfr_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_tepfr_fa_LSF_1_LE_p4.mat')

for i = 1:numel(mean_tepfr_fa_LSF_1_LE_0)
    mean_tepfr_fa_LSF_1(i,1) = mean_tepfr_fa_LSF_1_LE_0(i,1);
end
for i = 1:numel(mean_tepfr_fa_LSF_1_LE_p2)
    mean_tepfr_fa_LSF_1(i,2) = mean_tepfr_fa_LSF_1_LE_p2(i,1);
end
for i = 1:numel(mean_tepfr_fa_LSF_1_LE_p4)
    mean_tepfr_fa_LSF_1(i,3) = mean_tepfr_fa_LSF_1_LE_p4(i,1);
end

clear mean_tepfr_fa_LSF_1_LE_0
clear mean_tepfr_fa_LSF_1_LE_p2
clear mean_tepfr_fa_LSF_1_LE_p4

%Mean of cost per full run -- for stats purposes
load('../CuratedData_MPC/fa/LSF_1/mean_cost_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_cost_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_cost_fa_LSF_1_LE_p4.mat')

for i = 1:numel(mean_cost_fa_LSF_1_LE_0)
    mean_cost_fa_LSF_1(i,1) = mean_cost_fa_LSF_1_LE_0(i,1);
end
for i = 1:numel(mean_cost_fa_LSF_1_LE_p2)
    mean_cost_fa_LSF_1(i,2) = mean_cost_fa_LSF_1_LE_p2(i,1);
end
for i = 1:numel(mean_cost_fa_LSF_1_LE_p4)
    mean_cost_fa_LSF_1(i,3) = mean_cost_fa_LSF_1_LE_p4(i,1);
end

clear mean_cost_fa_LSF_1_LE_0
clear mean_cost_fa_LSF_1_LE_p2
clear mean_cost_fa_LSF_1_LE_p4

% %LSF_2
% %tepfr -- tracking error per full run
% load('../CuratedData_MPC/fa/LSF_2/mean_tepfr_fa_LSF_2_LE_0.mat')
% load('../CuratedData_MPC/fa/LSF_2/mean_tepfr_fa_LSF_2_LE_p2.mat')
% load('../CuratedData_MPC/fa/LSF_2/mean_tepfr_fa_LSF_2_LE_p4.mat')
% 
% for i = 1:numel(mean_tepfr_fa_LSF_2_LE_0)
%     mean_tepfr_fa_LSF_2(i,1) = mean_tepfr_fa_LSF_2_LE_0(i,1);
% end
% for i = 1:numel(mean_tepfr_fa_LSF_2_LE_p2)
%     mean_tepfr_fa_LSF_2(i,2) = mean_tepfr_fa_LSF_2_LE_p2(i,1);
% end
% for i = 1:numel(mean_tepfr_fa_LSF_2_LE_p4)
%     mean_tepfr_fa_LSF_2(i,3) = mean_tepfr_fa_LSF_2_LE_p4(i,1);
% end
% 
% clear mean_tepfr_fa_LSF_2_LE_0
% clear mean_tepfr_fa_LSF_2_LE_p2
% clear mean_tepfr_fa_LSF_2_LE_p4
% 
% %Mean of cost per full run -- for stats purposes
% load('../CuratedData_MPC/fa/LSF_2/mean_cost_fa_LSF_2_LE_0.mat')
% load('../CuratedData_MPC/fa/LSF_2/mean_cost_fa_LSF_2_LE_p2.mat')
% load('../CuratedData_MPC/fa/LSF_2/mean_cost_fa_LSF_2_LE_p4.mat')
% 
% for i = 1:numel(mean_cost_fa_LSF_2_LE_0)
%     mean_cost_fa_LSF_2(i,1) = mean_cost_fa_LSF_2_LE_0(i,1);
% end
% for i = 1:numel(mean_cost_fa_LSF_2_LE_p2)
%     mean_cost_fa_LSF_2(i,2) = mean_cost_fa_LSF_2_LE_p2(i,1);
% end
% for i = 1:numel(mean_cost_fa_LSF_2_LE_p4)
%     mean_cost_fa_LSF_2(i,3) = mean_cost_fa_LSF_2_LE_p4(i,1);
% end
% 
% clear mean_cost_fa_LSF_2_LE_0
% clear mean_cost_fa_LSF_2_LE_p2
% clear mean_cost_fa_LSF_2_LE_p4

disp('Curated data is loaded')

%% Calculate the mean and standard deviations for each LE/LSF pair
tepfr_meanvals = nan(numel(LSF_spread), numel(LE_spread),1);
cost_meanvals = nan(numel(LSF_spread), numel(LE_spread),1);

tepfr_stdevvals = nan(numel(LSF_spread), numel(LE_spread),1);
cost_stdevvals = nan(numel(LSF_spread), numel(LE_spread),1);

%Temporary mean of mean_cost_con vectors
mean_mtc_LSF_p5 = nan(1,numel(LE_spread));
mean_mtc_LSF_1 = nan(1,numel(LE_spread));
% mean_mtc_LSF_2 = nan(1,numel(LE_spread));

mean_mcc_LSF_p5 = nan(1,numel(LE_spread));
mean_mcc_LSF_1 = nan(1,numel(LE_spread));
% mean_mcc_LSF_2 = nan(1,numel(LE_spread));

%Temporary standard deviation of std_cost_con vectors
std_mtc_LSF_p5 = nan(1,numel(LE_spread));
std_mtc_LSF_1 = nan(1,numel(LE_spread));
% std_mtc_LSF_2 = nan(1,numel(LE_spread));

std_mcc_LSF_p5 = nan(1,numel(LE_spread));
std_mcc_LSF_1 = nan(1,numel(LE_spread));
% std_mcc_LSF_2 = nan(1,numel(LE_spread));


for i = 1:numel(LE_spread)
    %Means
    mean_mtc_LSF_p5(1,i) = nanmean(mean_tepfr_fa_LSF_p5(:,i));
    mean_mcc_LSF_p5(1,i) = nanmean(mean_cost_fa_LSF_p5(:,i));
    mean_mtc_LSF_1(1,i) = nanmean(mean_tepfr_fa_LSF_1(:,i));
    mean_mcc_LSF_1(1,i) = nanmean(mean_cost_fa_LSF_1(:,i));
%     mean_mtc_LSF_2(1,i) = nanmean(mean_tepfr_fa_LSF_2(:,i));
%     mean_mcc_LSF_2(1,i) = nanmean(mean_cost_fa_LSF_2(:,i));

    %Standard deviations
    std_mtc_LSF_p5(1,i) = nanstd(mean_tepfr_fa_LSF_p5(:,i));
    std_mcc_LSF_p5(1,i) = nanstd(mean_cost_fa_LSF_p5(:,i));
    std_mtc_LSF_1(1,i) = nanstd(mean_tepfr_fa_LSF_1(:,i));
    std_mcc_LSF_1(1,i) = nanstd(mean_cost_fa_LSF_1(:,i));
%     std_mtc_LSF_2(1,i) = nanstd(mean_tepfr_fa_LSF_2(:,i));
%     std_mcc_LSF_2(1,i) = nanstd(mean_cost_fa_LSF_2(:,i));  
end

%Assign the calculated means and standard deviations into their appropriate
%matrices

%Tracking error per full run 
tepfr_meanvals(1,:) = mean_mtc_LSF_p5;
tepfr_meanvals(2,:) = mean_mtc_LSF_1;
% tepfr_meanvals(3,:) = mean_mtc_LSF_2;

tepfr_stdevvals(1,:) = std_mtc_LSF_p5;
tepfr_stdevvals(2,:) = std_mtc_LSF_1;
% tepfr_stdevvals(3,:) = std_mtc_LSF_2;

%Cost function value
cost_meanvals(1,:) = mean_mcc_LSF_p5;
cost_meanvals(2,:) = mean_mcc_LSF_1;
% cost_meanvals(3,:) = mean_mcc_LSF_2;

cost_stdevvals(1,:) = std_mcc_LSF_p5;
cost_stdevvals(2,:) = std_mcc_LSF_1;
% cost_stdevvals(3,:) = std_mcc_LSF_2;

%% Tracking error per full run surface plot

fig1 = figure(1);
subplot(1,2,1)
surfc(tepfr_meanvals);
title('fa')
zlabel('Mean of tracking error (cm)'); 
xlabel('Petiole length extention');
ylabel('Length scale factor');
xticks(LE_spread)
xticklabels({'0','0.2','0.4'})
yticks(LSF_spread)
% yticklabels({'0.5','1','2'})
yticklabels({'0.5','1'})

subplot(1,2,2)
contour(tepfr_meanvals,'ShowText','on','LineWidth',3)
title('fa')
xlabel('Petiole length extention');
ylabel('Length scale factor');
xticks(LE_spread)
xticklabels({'0','0.2','0.4'})
yticks(LSF_spread)
% yticklabels({'0.5','1','2'})
yticklabels({'0.5','1'})

set(fig1, 'Position', [0,0, 1200, 2400])

saveas(gcf,'../Figures_MPC/fa/compSize/Fig18_fa_tepfr_surf_LSFp5-1_12h.fig');
saveas(gcf,'../Figures_MPC/fa/compSize/Fig18_fa_tepfr_surf_LSFp5-1_12h.jpg');

disp('Figure 1 done')

%% Tracking error per full run surface plot

fig2 = figure(2);
subplot(1,2,1)
surfc(tepfr_stdevvals);
title('fa')
zlabel('Standard deviation of tracking error (cm)'); 
xlabel('Petiole length extention');
ylabel('Length scale factor');
xticks(LE_spread)
xticklabels({'0','0.2','0.4'})
yticks(LSF_spread)
% yticklabels({'0.5','1','2'})
yticklabels({'0.5','1'})

subplot(1,2,2)
contour(tepfr_stdevvals,'ShowText','on','LineWidth',3)
title('fa')
xlabel('Petiole length extention');
ylabel('Length scale factor');
xticks(LE_spread)
xticklabels({'0','0.2','0.4'})
yticks(LSF_spread)
% yticklabels({'0.5','1','2'})
yticklabels({'0.5','1'})

set(fig2, 'Position', [0,0, 1200, 2400])

saveas(gcf,'../Figures_MPC/fa/compSize/Fig19_fa_tepfr_surf_std_LSFp5-1_12h.fig');
saveas(gcf,'../Figures_MPC/fa/compSize/Fig19_fa_tepfr_surf_std_LSFp5-1_12h.jpg');

disp('Figure 2 done')

%% Tracking error per full run surface plot

fig3 = figure(3);
subplot(1,2,1)
surfc(cost_meanvals);
title('fa')
zlabel('Mean of cost function value'); 
xlabel('Petiole length extention');
ylabel('Length scale factor');
xticks(LE_spread)
xticklabels({'0','0.2','0.4'})
yticks(LSF_spread)
% yticklabels({'0.5','1','2'})
yticklabels({'0.5','1'})

subplot(1,2,2)
contour(cost_meanvals,'ShowText','on','LineWidth',3)
title('fa')
xlabel('Petiole length extention');
ylabel('Length scale factor');
xticks(LE_spread)
xticklabels({'0','0.2','0.4'})
yticks(LSF_spread)
% yticklabels({'0.5','1','2'})
yticklabels({'0.5','1'})

set(fig3, 'Position', [0,0, 1200, 2400])

saveas(gcf,'../Figures_MPC/fa/compSize/Fig20_fa_cost_surf_LSFp5-1_12h.fig');
saveas(gcf,'../Figures_MPC/fa/compSize/Fig20_fa_cost_surf_LSFp5-1_12h.jpg');

disp('Figure 3 done')

%% Tracking error per full run surface plot

fig4 = figure(4);
subplot(1,2,1)
surfc(cost_stdevvals);
title('fa')
zlabel('Standard deviation of cost function value'); 
xlabel('Petiole length extention');
ylabel('Length scale factor');
xticks(LE_spread)
xticklabels({'0','0.2','0.4'})
yticks(LSF_spread)
% yticklabels({'0.5','1','2'})
yticklabels({'0.5','1'})

subplot(1,2,2)
contour(cost_stdevvals,'ShowText','on','LineWidth',3)
title('fa')
xlabel('Petiole length extention');
ylabel('Length scale factor');
xticks(LE_spread)
xticklabels({'0','0.2','0.4'})
yticks(LSF_spread)
% yticklabels({'0.5','1','2'})
yticklabels({'0.5','1'})

set(fig4, 'Position', [0,0, 1200, 2400])

saveas(gcf,'../Figures_MPC/fa/compSize/Fig21_fa_cost_surf_std_LSFp5-1_12h.fig');
saveas(gcf,'../Figures_MPC/fa/compSize/Fig21_fa_cost_surf_std_LSFp5-1_12h.jpg');

disp('Figure 4 done')

%% Final window
disp('Completely done, son')
disp(['Time at end: ',datestr(now)])