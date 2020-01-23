%8/8/17 Script developed to curate the data into a form more usable by
%error_11h
%8/16/18 Script modified for LengthScaleFactor (LSF)
%10/17/18 Script modified for size extension (LE)
%1/11/19 Script modified for Model Predictive Control (MPC)
    %1/16/19 Script corrected previous indexing error.
    
%v11b is for horizontal aggressive maneuver
%v11c is for aggressive maneuver
%v12h is for sum of prime number sines

%% We the People, in Order to form a perfect Union,...
clc;
clearvars;
close all;
disp(['Time at start: ',datestr(now)])

%% Assign moth variables IF NECESSARY
% disp('Assign moth parameters')

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

%Create the sum of primes signal
signal_amp = 5; %in cm 
prime_f = [0.2, 0.3, 0.5, 0.7, 1.1, 1.7, 2.9, 4.3, 7.9, 13.7, 19.9]; %in Hz
prime_a = (signal_amp./(2*pi.*prime_f)).*(2*pi.*prime_f(1)); %in cm
prime_ph = zeros(1,numel(prime_f)); %Not sure if we'll require a phase

%Goal criteria (for the cost function)
y_g = prime_a(1)*sin(2*pi*prime_f(1)*Tstore(1:(end-100+1)) + prime_ph(1)) +...
    prime_a(2)*sin(2*pi*prime_f(2)*Tstore(1:(end-100+1)) + prime_ph(2)) +...
    prime_a(3)*sin(2*pi*prime_f(3)*Tstore(1:(end-100+1)) + prime_ph(3)) +...
    prime_a(4)*sin(2*pi*prime_f(4)*Tstore(1:(end-100+1)) + prime_ph(4)) +...
    prime_a(5)*sin(2*pi*prime_f(5)*Tstore(1:(end-100+1)) + prime_ph(5)) +...
    prime_a(6)*sin(2*pi*prime_f(6)*Tstore(1:(end-100+1)) + prime_ph(6)) +...
    prime_a(7)*sin(2*pi*prime_f(7)*Tstore(1:(end-100+1)) + prime_ph(7)) +...
    prime_a(8)*sin(2*pi*prime_f(8)*Tstore(1:(end-100+1)) + prime_ph(8)) +...
    prime_a(9)*sin(2*pi*prime_f(9)*Tstore(1:(end-100+1)) + prime_ph(9)) +...
    prime_a(10)*sin(2*pi*prime_f(10)*Tstore(1:(end-100+1)) + prime_ph(10)) +...
    prime_a(11)*sin(2*pi*prime_f(11)*Tstore(1:(end-100+1)) + prime_ph(11));

ydot_g = 2*pi*prime_a(1)*prime_f(1)*cos(2*pi*prime_f(1)*Tstore(1:(end-100+1)) + prime_ph(1)) +...
    2*pi*prime_a(2)*prime_f(2)*cos(2*pi*prime_f(2)*Tstore(1:(end-100+1)) + prime_ph(2)) +...
    2*pi*prime_a(3)*prime_f(3)*cos(2*pi*prime_f(3)*Tstore(1:(end-100+1)) + prime_ph(3)) +...
    2*pi*prime_a(4)*prime_f(4)*cos(2*pi*prime_f(4)*Tstore(1:(end-100+1)) + prime_ph(4)) +...
    2*pi*prime_a(5)*prime_f(5)*cos(2*pi*prime_f(5)*Tstore(1:(end-100+1)) + prime_ph(5)) +...
    2*pi*prime_a(6)*prime_f(6)*cos(2*pi*prime_f(6)*Tstore(1:(end-100+1)) + prime_ph(6)) +...
    2*pi*prime_a(7)*prime_f(7)*cos(2*pi*prime_f(7)*Tstore(1:(end-100+1)) + prime_ph(7)) +...
    2*pi*prime_a(8)*prime_f(8)*cos(2*pi*prime_f(8)*Tstore(1:(end-100+1)) + prime_ph(8)) +...
    2*pi*prime_a(9)*prime_f(9)*cos(2*pi*prime_f(9)*Tstore(1:(end-100+1)) + prime_ph(9)) +...
    2*pi*prime_a(10)*prime_f(10)*cos(2*pi*prime_f(10)*Tstore(1:(end-100+1)) + prime_ph(10)) +...
    2*pi*prime_a(11)*prime_f(11)*cos(2*pi*prime_f(11)*Tstore(1:(end-100+1)) + prime_ph(11));

y_g = y_g';
x_g = zeros(1,numel(y_g));
theta_g = zeros(1,numel(y_g));
theta_g(1,:) = pi/4;

ydot_g = ydot_g';
xdot_g = zeros(1,numel(y_g));
thetadot_g = zeros(1,numel(y_g));

NumOfLSF_1_LE_0_Files = numel(listOfLSF_1_LE_0_Files);
NumOfLSF_1_LE_p2_Files = numel(listOfLSF_1_LE_p2_Files);
NumOfLSF_1_LE_p4_Files = numel(listOfLSF_1_LE_p4_Files);

%% Import all the variables
disp('Importing curated data')

%con - LSF_1_LE_0
load('../CuratedData_MPC/fa/LSF_1/mean_x_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_y_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_theta_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_phi_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_xdot_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_ydot_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_thetadot_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_phidot_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_beta_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_dist_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/std_x_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/std_y_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/std_theta_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/std_phi_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/std_xdot_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/std_ydot_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/std_thetadot_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/std_phidot_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/std_beta_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/std_dist_fa_LSF_1_LE_0.mat')

load('../CuratedData_MPC/fa/LSF_1/Winstore_x_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_y_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_theta_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_phi_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_beta_fa_LSF_1_LE_0.mat')

%con - LSF_1_LE_p2
load('../CuratedData_MPC/fa/LSF_1/mean_x_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_y_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_theta_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_phi_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_xdot_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_ydot_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_thetadot_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_phidot_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_beta_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_dist_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/std_x_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/std_y_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/std_theta_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/std_phi_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/std_xdot_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/std_ydot_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/std_thetadot_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/std_phidot_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/std_beta_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/std_dist_fa_LSF_1_LE_p2.mat')

load('../CuratedData_MPC/fa/LSF_1/Winstore_x_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_y_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_theta_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_phi_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_beta_fa_LSF_1_LE_p2.mat')

%con - LSF_1_LE_p4
load('../CuratedData_MPC/fa/LSF_1/mean_x_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_y_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_theta_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_phi_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_xdot_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_ydot_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_thetadot_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_phidot_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_beta_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/mean_dist_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/std_x_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/std_y_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/std_theta_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/std_phi_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/std_xdot_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/std_ydot_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/std_thetadot_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/std_phidot_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/std_beta_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/std_dist_fa_LSF_1_LE_p4.mat')

load('../CuratedData_MPC/fa/LSF_1/Winstore_x_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_y_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_theta_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_phi_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_beta_fa_LSF_1_LE_p4.mat')

disp('Curated data is loaded')

%% Find the min & max values for...
%x
xMax_vec = [max(max(mean_x_fa_LSF_1_LE_0 + std_x_fa_LSF_1_LE_0)),...
    max(max(mean_x_fa_LSF_1_LE_p2 + std_x_fa_LSF_1_LE_p2)),...
    max(max(mean_x_fa_LSF_1_LE_p4 + std_x_fa_LSF_1_LE_p4))];
x_max = max(xMax_vec);

xMin_vec = [min(min(mean_x_fa_LSF_1_LE_0 - std_x_fa_LSF_1_LE_0)),...
    min(min(mean_x_fa_LSF_1_LE_p2 - std_x_fa_LSF_1_LE_p2)),...
    min(min(mean_x_fa_LSF_1_LE_p4 - std_x_fa_LSF_1_LE_p4))];
x_min = min(xMin_vec);

%y
yMax_vec = [max(max(mean_y_fa_LSF_1_LE_0 + std_y_fa_LSF_1_LE_0)),...
    max(max(mean_y_fa_LSF_1_LE_p2 + std_y_fa_LSF_1_LE_p2)),...
    max(max(mean_y_fa_LSF_1_LE_p4 + std_y_fa_LSF_1_LE_p4))];
y_max = max(yMax_vec);

yMin_vec = [min(min(mean_y_fa_LSF_1_LE_0 - std_y_fa_LSF_1_LE_0)),...
    min(min(mean_y_fa_LSF_1_LE_p2 - std_y_fa_LSF_1_LE_p2)),...
    min(min(mean_y_fa_LSF_1_LE_p4 - std_y_fa_LSF_1_LE_p4))];
y_min = min(yMin_vec);

%theta
thetaMax_vec = [max(max(mean_theta_fa_LSF_1_LE_0 + std_theta_fa_LSF_1_LE_0)),...
    max(max(mean_theta_fa_LSF_1_LE_p2 + std_theta_fa_LSF_1_LE_p2)),...
    max(max(mean_theta_fa_LSF_1_LE_p4 + std_theta_fa_LSF_1_LE_p4))];
theta_max = max(thetaMax_vec);

thetaMin_vec = [min(min(mean_theta_fa_LSF_1_LE_0 - std_theta_fa_LSF_1_LE_0)),...
    min(min(mean_theta_fa_LSF_1_LE_p2 - std_theta_fa_LSF_1_LE_p2)),...
    min(min(mean_theta_fa_LSF_1_LE_p4 - std_theta_fa_LSF_1_LE_p4))];
theta_min = min(thetaMin_vec);

%phi
phiMax_vec = [max(max(mean_phi_fa_LSF_1_LE_0 + std_phi_fa_LSF_1_LE_0)),...
    max(max(mean_phi_fa_LSF_1_LE_p2 + std_phi_fa_LSF_1_LE_p2)),...
    max(max(mean_phi_fa_LSF_1_LE_p4 + std_phi_fa_LSF_1_LE_p4))];
phi_max = max(phiMax_vec);

phiMin_vec = [min(min(mean_phi_fa_LSF_1_LE_0 - std_phi_fa_LSF_1_LE_0)),...
    min(min(mean_phi_fa_LSF_1_LE_p2 - std_phi_fa_LSF_1_LE_p2)),...
    min(min(mean_phi_fa_LSF_1_LE_p4 - std_phi_fa_LSF_1_LE_p4))];
phi_min = min(phiMin_vec);

%beta
betaMax_vec = [max(max(mean_beta_fa_LSF_1_LE_0 + std_beta_fa_LSF_1_LE_0)),...
    max(max(mean_beta_fa_LSF_1_LE_p2 + std_beta_fa_LSF_1_LE_p2)),...
    max(max(mean_beta_fa_LSF_1_LE_p4 + std_beta_fa_LSF_1_LE_p4))];
beta_max = max(betaMax_vec);

betaMin_vec = [min(min(mean_beta_fa_LSF_1_LE_0 - std_beta_fa_LSF_1_LE_0)),...
    min(min(mean_beta_fa_LSF_1_LE_p2 - std_beta_fa_LSF_1_LE_p2)),...
    min(min(mean_beta_fa_LSF_1_LE_p4 - std_beta_fa_LSF_1_LE_p4))];
beta_min = min(betaMin_vec);

%% Set up internal values for Figure 1

abdoflip(1:numel(Tstore)) = 2*pi;
livebetarange(1:numel(Tstore)) = 20*(pi/180); 

%% State variables w.r.t. time - individfally
fig1 = figure(1);

%x - con LSF_1_LE_0
subplot(5,3,1)
if NumOfLSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),Winstore_x_fa_LSF_1_LE_0,'LineWidth',2);
end
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0'});
ylabel({'x';'(cm)'});
% axis([0, Tstore(end-100+1), x_min, x_max])
axis([0, Tstore(end-100+1), -5, 5])

%x - con LSF_1_LE_p2
subplot(5,3,2)
if NumOfLSF_1_LE_p2_Files > 0
    plot(Tstore(1:(end-100+1)),Winstore_x_fa_LSF_1_LE_p2,'LineWidth',2);
end
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0.2'});
% axis([0, Tstore(end-100+1), x_min, x_max])
axis([0, Tstore(end-100+1), -5, 5])

%x - con LSF_1_LE_p4
subplot(5,3,3)
if NumOfLSF_1_LE_p4_Files > 0
    plot(Tstore(1:(end-100+1)),Winstore_x_fa_LSF_1_LE_p4,'LineWidth',2);
end
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0.4'});
% axis([0, Tstore(end-100+1), x_min, x_max])
axis([0, Tstore(end-100+1), -5, 5])

%y - con LSF_1_LE_0
subplot(5,3,4)
if NumOfLSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),Winstore_y_fa_LSF_1_LE_0,'LineWidth',2);
    hold on;
end
ylabel({'y';'(cm)'});
plot(Tstore(1:(end-100+1)), y_g,'k--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])
axis([0, Tstore(end-100+1), -11, 11])

%y - con LSF_1_LE_p2
subplot(5,3,5)
if NumOfLSF_1_LE_p2_Files > 0
    plot(Tstore(1:(end-100+1)),Winstore_y_fa_LSF_1_LE_p2,'LineWidth',2);
    hold on;
end
plot(Tstore(1:(end-100+1)), y_g,'k--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])
axis([0, Tstore(end-100+1), -11, 11])

%y - con LSF_1_LE_p4
subplot(5,3,6)
if NumOfLSF_1_LE_p4_Files > 0
    plot(Tstore(1:(end-100+1)),Winstore_y_fa_LSF_1_LE_p4,'LineWidth',2);
    hold on;
end
plot(Tstore(1:(end-100+1)), y_g,'k--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])
axis([0, Tstore(end-100+1), -11, 11])

%theta - con LSF_1_LE_0
subplot(5,3,7)
if NumOfLSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_theta_fa_LSF_1_LE_0,'LineWidth',2);
end
ylabel({'theta';'(degrees)'});
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
axis([0, Tstore(end-100+1), -20, 110])

%theta - con LSF_1_LE_p2
subplot(5,3,8)
if NumOfLSF_1_LE_p2_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_theta_fa_LSF_1_LE_p2,'LineWidth',2);
end
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
axis([0, Tstore(end-100+1), -20, 110])

%theta - con LSF_1_LE_p4
subplot(5,3,9)
if NumOfLSF_1_LE_p4_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_theta_fa_LSF_1_LE_p4,'LineWidth',2);
end
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
axis([0, Tstore(end-100+1), -inf, inf])

%phi - con LSF_1_LE_0
subplot(5,3,10)
if NumOfLSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_phi_fa_LSF_1_LE_0,'LineWidth',2);
    hold on;
end
plot(Tstore, (180/pi)*abdoflip,'k--','LineWidth',1)
ylabel({'phi';'(degrees)'});
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])
axis([0, Tstore(end-100+1), -270, 800])

%phi - con LSF_1_LE_p2
subplot(5,3,11)
if NumOfLSF_1_LE_p2_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_phi_fa_LSF_1_LE_p2,'LineWidth',2);
    hold on;
end
plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])
axis([0, Tstore(end-100+1), -270, 800])

%phi - con LSF_1_LE_p4
subplot(5,3,12)
if NumOfLSF_1_LE_p4_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_phi_fa_LSF_1_LE_p4,'LineWidth',2);
    hold on;
end
plot(Tstore, (180/pi)*abdoflip,'k--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])
axis([0, Tstore(end-100+1), -inf, inf])

%beta - con LSF_1_LE_0
subplot(5,3,13)
if NumOfLSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_beta_fa_LSF_1_LE_0,'LineWidth',2);
    hold on;
end
xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOfLSF_1_LE_0_Files),')']}); 
ylabel({'beta';'(degrees)'});
plot(Tstore, (180/pi)*livebetarange,'k--','LineWidth',1)
hold on;
plot(Tstore, -(180/pi)*livebetarange,'k--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',0:40:80);
axis([0, Tstore(end-100+1), -500, 500])

%beta - con LSF_1_LE_p2
subplot(5,3,14)
if NumOfLSF_1_LE_p2_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_beta_fa_LSF_1_LE_p2,'LineWidth',2);
    hold on;
end
xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOfLSF_1_LE_p2_Files),')']});
plot(Tstore, (180/pi)*livebetarange,'k--','LineWidth',1)
hold on;
plot(Tstore, -(180/pi)*livebetarange,'k--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',0:40:80);
axis([0, Tstore(end-100+1), -500, 500])

%beta - con LSF_1_LE_p4
subplot(5,3,15)
if NumOfLSF_1_LE_p4_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_beta_fa_LSF_1_LE_p4,'LineWidth',2);
    hold on;
end
xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOfLSF_1_LE_p4_Files),')']}); 
plot(Tstore, (180/pi)*livebetarange,'k--','LineWidth',1)
hold on;
plot(Tstore, -(180/pi)*livebetarange,'k--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',0:40:80);
axis([0, Tstore(end-100+1), -inf, inf])

set(fig1, 'Position', [0,0, 1000, 600])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig1_fa_stateVarsIndiv_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig1_fa_stateVarsIndiv_LSF_1_v12h.jpg');

disp('Figure 1 done')

%% State variables w.r.t. time - mean
fig2 = figure(2);

%x - con LSF_1_LE_0
subplot(5,3,1)
Fig2_1 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_fa_LSF_1_LE_0,std_x_fa_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_1.mainLine.LineWidth = 2;
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0'})
ylabel({'x';'(cm)'});
% axis([0, Tstore(end-100+1), x_min, x_max])
axis([0, Tstore(end-100+1), -2, 2])

%x - con LSF_1_LE_p2
subplot(5,3,2)
Fig2_2 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_fa_LSF_1_LE_p2,...
    std_x_fa_LSF_1_LE_p2,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_2.mainLine.LineWidth = 2;
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0.2'})
% axis([0, Tstore(end-100+1), x_min, x_max])
axis([0, Tstore(end-100+1), -2, 2])

%x - con LSF_1_LE_p4
subplot(5,3,3)
Fig2_3 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_fa_LSF_1_LE_p4,...
    std_x_fa_LSF_1_LE_p4,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_3.mainLine.LineWidth = 2;
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0.4'})
% axis([0, Tstore(end-100+1), x_min, x_max])
axis([0, Tstore(end-100+1), -2, 2])

%y - con LSF_1_LE_0
subplot(5,3,4)
Fig2_8 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_fa_LSF_1_LE_0,...
    std_y_fa_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_8.mainLine.LineWidth = 1.2; 
hold on;
ylabel({'y';'(cm)'});
plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])
axis([0, Tstore(end-100+1), -11, 11])

%y - con LSF_1_LE_p2
subplot(5,3,5)
Fig2_9 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_fa_LSF_1_LE_p2,...
    std_y_fa_LSF_1_LE_p2,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_9.mainLine.LineWidth = 1.2; 
hold on;
plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])
axis([0, Tstore(end-100+1), -11, 11])

%y - con LSF_1_LE_p4
subplot(5,3,6)
Fig2_10 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_fa_LSF_1_LE_p4,...
    std_y_fa_LSF_1_LE_p4,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_10.mainLine.LineWidth = 1.2;
hold on;
plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])
axis([0, Tstore(end-100+1), -11, 11])

%theta - con LSF_1_LE_0
subplot(5,3,7)
Fig2_15 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_fa_LSF_1_LE_0,...
    (180/pi)*std_theta_fa_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
ylabel({'theta';'(degrees)'});
Fig2_15.mainLine.LineWidth = 2;
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
axis([0, Tstore(end-100+1), 20, 70])

%theta - con LSF_1_LE_p2
subplot(5,3,8)
Fig2_16 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_fa_LSF_1_LE_p2,...
    (180/pi)*std_theta_fa_LSF_1_LE_p2,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_16.mainLine.LineWidth = 2;
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
axis([0, Tstore(end-100+1), 20, 70])

%theta - con LSF_1_LE_p4
subplot(5,3,9)
Fig2_17 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_fa_LSF_1_LE_p4,...
    (180/pi)*std_theta_fa_LSF_1_LE_p4,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_17.mainLine.LineWidth = 2;
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
axis([0, Tstore(end-100+1), 20, 70])

%phi - con LSF_1_LE_0
subplot(5,3,10)
Fig2_22 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_phi_fa_LSF_1_LE_0,...
    (180/pi)*std_phi_fa_LSF_1_LE_0,...
    'lineprops', '-k','transparent',false,'patchSaturation',0.075);
Fig2_22.mainLine.LineWidth = 2;
hold on;
plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
ylabel({'phi';'(degrees)'});
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])
axis([0, Tstore(end-100+1), 100, 420])

%phi - con LSF_1_LE_p2
subplot(5,3,11)
Fig2_23 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_phi_fa_LSF_1_LE_p2,...
    (180/pi)*std_phi_fa_LSF_1_LE_p2,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_23.mainLine.LineWidth = 2;
hold on;
plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])
axis([0, Tstore(end-100+1), 100, 420])

%phi - con LSF_1_LE_p4
subplot(5,3,12)
Fig2_24 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_phi_fa_LSF_1_LE_p4,...
    (180/pi)*std_phi_fa_LSF_1_LE_p4,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_24.mainLine.LineWidth = 2;
hold on;
plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])
axis([0, Tstore(end-100+1), -inf, inf])

%beta - con LSF_1_LE_0
subplot(5,3,13)
Fig2_29 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_fa_LSF_1_LE_0,...
    (180/pi)*std_beta_fa_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_29.mainLine.LineWidth = 2;
xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOfLSF_1_LE_0_Files),')']}); 
ylabel({'beta';'(degrees)'});
hold on;
plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
hold on;
plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',0:40:80);
axis([0, Tstore(end-100+1), -200, 200])

%beta - con LSF_1_LE_p2
subplot(5,3,14)
Fig2_30 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_fa_LSF_1_LE_p2,...
    (180/pi)*std_beta_fa_LSF_1_LE_p2,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_30.mainLine.LineWidth = 2;
xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOfLSF_1_LE_p2_Files),')']}); 
hold on;
plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
hold on;
plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',0:40:80);
axis([0, Tstore(end-100+1), -200, 200])

%beta - con LSF_1_LE_p4
subplot(5,3,15)
Fig2_31 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_fa_LSF_1_LE_p4,...
    (180/pi)*std_beta_fa_LSF_1_LE_p4,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);   
Fig2_31.mainLine.LineWidth = 2;
xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOfLSF_1_LE_p4_Files),')']}); 
hold on;
plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
hold on;
plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',0:40:80);
axis([0, Tstore(end-100+1), -inf, inf])

set(fig2, 'Position', [0,0, 1000, 600])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig2_fa_stateVarsMean_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig2_fa_stateVarsMean_LSF_1_v12h.jpg');

disp('Figure 2 done')

%% Tracking error w.r.t. time

% figure(3);
% %y motion- goal
% subplot(2,1,1)
% plot(Tstore(1:(end-100+1)),y_g,'-k','LineWidth',2);
% ylabel({'Goal y-motion'; '(cm)'}); 
% % xlabel({'0'; ['(n = ',num2str(NumOfconLSF_1_LE_0_Files),')']});
% axis([0, Tstore(end-100+1), -inf,inf])
% 
% %Tracking error - con LSF_1_LE_0
% subplot(2,1,2)
% Fig2_1 = shadedErrorBar(Tstore(1:(end-100+1)),mean_dist_fa_LSF_1_LE_0,...
%     std_dist_fa_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% Fig2_1.mainLine.LineWidth = 2;
% ylabel({'Tracking error'; '(cm)'}); 
% xlabel({'fa';'LE: 0'; ['(n = ',num2str(NumOfLSF_1_LE_0_Files),')']});
% axis([0, Tstore(end-100+1), -inf,inf])
% % set(gca,'ytick',0:0.1:0.6);
% 
% saveas(gcf,'../Figures_MPC/fa/LSF_1/FigOPT_fa_error_time_LSF_1_v12h.fig');
% saveas(gcf,'../Figures_MPC/fa/LSF_1/FigOPT_fa_error_time_LSF_1_v12h.jpg');
% 
% disp('OPTIONAL - Figure of tracking error w.r.t. time done')

%% Final window
disp('Completely done, son')
disp(['Time at end: ',datestr(now)])