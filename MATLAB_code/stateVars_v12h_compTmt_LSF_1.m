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

NumOf_fa_LSF_1_LE_0_Files = numel(listOf_fa_LSF_1_LE_0_Files);
NumOf_fs_LSF_1_LE_0_Files = numel(listOf_fs_LSF_1_LE_0_Files);
NumOf_ua_LSF_1_LE_0_Files = numel(listOf_ua_LSF_1_LE_0_Files);
NumOf_us_LSF_1_LE_0_Files = numel(listOf_us_LSF_1_LE_0_Files);
NumOf_uaW_LSF_1_LE_0_Files = numel(listOf_uaW_LSF_1_LE_0_Files);
NumOf_usW_LSF_1_LE_0_Files = numel(listOf_usW_LSF_1_LE_0_Files);

%% Import all the variables
disp('Importing curated data')

%fa - LSF_1_LE_0
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

%fs - LSF_1_LE_0
load('../CuratedData_MPC/fs/LSF_1/mean_x_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/mean_y_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/mean_theta_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/mean_phi_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/mean_xdot_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/mean_ydot_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/mean_thetadot_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/mean_phidot_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/mean_beta_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/mean_dist_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/std_x_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/std_y_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/std_theta_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/std_phi_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/std_xdot_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/std_ydot_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/std_thetadot_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/std_phidot_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/std_beta_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/std_dist_fs_LSF_1_LE_0.mat')

load('../CuratedData_MPC/fs/LSF_1/Winstore_x_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/Winstore_y_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/Winstore_theta_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/Winstore_phi_fs_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fs/LSF_1/Winstore_beta_fs_LSF_1_LE_0.mat')

%ua - LSF_1_LE_0
load('../CuratedData_MPC/ua/LSF_1/mean_x_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/mean_y_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/mean_theta_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/mean_phi_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/mean_xdot_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/mean_ydot_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/mean_thetadot_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/mean_phidot_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/mean_beta_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/mean_dist_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/std_x_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/std_y_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/std_theta_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/std_phi_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/std_xdot_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/std_ydot_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/std_thetadot_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/std_phidot_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/std_beta_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/std_dist_ua_LSF_1_LE_0.mat')

load('../CuratedData_MPC/ua/LSF_1/Winstore_x_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/Winstore_y_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/Winstore_theta_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/Winstore_phi_ua_LSF_1_LE_0.mat')
load('../CuratedData_MPC/ua/LSF_1/Winstore_beta_ua_LSF_1_LE_0.mat')

%us - LSF_1_LE_0
load('../CuratedData_MPC/us/LSF_1/mean_x_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/mean_y_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/mean_theta_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/mean_phi_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/mean_xdot_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/mean_ydot_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/mean_thetadot_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/mean_phidot_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/mean_beta_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/mean_dist_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/std_x_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/std_y_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/std_theta_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/std_phi_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/std_xdot_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/std_ydot_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/std_thetadot_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/std_phidot_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/std_beta_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/std_dist_us_LSF_1_LE_0.mat')

load('../CuratedData_MPC/us/LSF_1/Winstore_x_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/Winstore_y_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/Winstore_theta_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/Winstore_phi_us_LSF_1_LE_0.mat')
load('../CuratedData_MPC/us/LSF_1/Winstore_beta_us_LSF_1_LE_0.mat')

%uaW - LSF_1_LE_0
load('../CuratedData_MPC/uaW/LSF_1/mean_x_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/mean_y_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/mean_theta_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/mean_phi_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/mean_xdot_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/mean_ydot_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/mean_thetadot_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/mean_phidot_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/mean_beta_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/mean_dist_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/std_x_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/std_y_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/std_theta_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/std_phi_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/std_xdot_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/std_ydot_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/std_thetadot_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/std_phidot_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/std_beta_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/std_dist_uaW_LSF_1_LE_0.mat')

load('../CuratedData_MPC/uaW/LSF_1/Winstore_x_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/Winstore_y_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/Winstore_theta_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/Winstore_phi_uaW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/uaW/LSF_1/Winstore_beta_uaW_LSF_1_LE_0.mat')

%usW - LSF_1_LE_0
load('../CuratedData_MPC/usW/LSF_1/mean_x_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/mean_y_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/mean_theta_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/mean_phi_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/mean_xdot_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/mean_ydot_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/mean_thetadot_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/mean_phidot_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/mean_beta_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/mean_dist_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/std_x_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/std_y_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/std_theta_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/std_phi_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/std_xdot_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/std_ydot_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/std_thetadot_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/std_phidot_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/std_beta_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/std_dist_usW_LSF_1_LE_0.mat')

load('../CuratedData_MPC/usW/LSF_1/Winstore_x_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/Winstore_y_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/Winstore_theta_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/Winstore_phi_usW_LSF_1_LE_0.mat')
load('../CuratedData_MPC/usW/LSF_1/Winstore_beta_usW_LSF_1_LE_0.mat')

disp('Curated data is loaded')

%% Find the min & max values for...
%x
xMax_vec = [max(max(Winstore_x_fa_LSF_1_LE_0)),...
    max(max(Winstore_x_fs_LSF_1_LE_0)),...
    max(max(Winstore_x_ua_LSF_1_LE_0)),...
    max(max(Winstore_x_us_LSF_1_LE_0)),...
    max(max(Winstore_x_uaW_LSF_1_LE_0)),...
    max(max(Winstore_x_usW_LSF_1_LE_0))];
x_max = max(xMax_vec);

xMin_vec = [min(min(Winstore_x_fa_LSF_1_LE_0)),...
    min(min(Winstore_x_fs_LSF_1_LE_0)),...
    min(min(Winstore_x_ua_LSF_1_LE_0)),...
    min(min(Winstore_x_us_LSF_1_LE_0)),...
    min(min(Winstore_x_uaW_LSF_1_LE_0)),...
    min(min(Winstore_x_usW_LSF_1_LE_0))];
x_min = min(xMin_vec);

%y
yMax_vec = [max(max(Winstore_y_fa_LSF_1_LE_0)),...
    max(max(Winstore_y_fs_LSF_1_LE_0)),...
    max(max(Winstore_y_ua_LSF_1_LE_0)),...
    max(max(Winstore_y_us_LSF_1_LE_0)),...
    max(max(Winstore_y_uaW_LSF_1_LE_0)),...
    max(max(Winstore_y_usW_LSF_1_LE_0))];
y_max = max(yMax_vec);

yMin_vec = [min(min(Winstore_y_fa_LSF_1_LE_0)),...
    min(min(Winstore_y_fs_LSF_1_LE_0)),...
    min(min(Winstore_y_ua_LSF_1_LE_0)),...
    min(min(Winstore_y_us_LSF_1_LE_0)),...
    min(min(Winstore_y_uaW_LSF_1_LE_0)),...
    min(min(Winstore_y_usW_LSF_1_LE_0))];
y_min = min(yMin_vec);

%theta
thetaMax_vec = [max(max(Winstore_theta_fa_LSF_1_LE_0)),...
    max(max(Winstore_theta_fs_LSF_1_LE_0)),...
    max(max(Winstore_theta_ua_LSF_1_LE_0)),...
    max(max(Winstore_theta_us_LSF_1_LE_0)),...
    max(max(Winstore_theta_uaW_LSF_1_LE_0)),...
    max(max(Winstore_theta_usW_LSF_1_LE_0))];
theta_max = max(thetaMax_vec);

thetaMin_vec = [min(min(Winstore_theta_fa_LSF_1_LE_0)),...
    min(min(Winstore_theta_fs_LSF_1_LE_0)),...
    min(min(Winstore_theta_ua_LSF_1_LE_0)),...
    min(min(Winstore_theta_us_LSF_1_LE_0)),...
    min(min(Winstore_theta_uaW_LSF_1_LE_0)),...
    min(min(Winstore_theta_usW_LSF_1_LE_0))];
theta_min = min(thetaMin_vec);

%phi
phiMax_vec = [max(max(Winstore_phi_fa_LSF_1_LE_0)),...
    max(max(Winstore_phi_fs_LSF_1_LE_0)),...
    max(max(Winstore_phi_ua_LSF_1_LE_0)),...
    max(max(Winstore_phi_us_LSF_1_LE_0)),...
    max(max(Winstore_phi_uaW_LSF_1_LE_0)),...
    max(max(Winstore_phi_usW_LSF_1_LE_0))];
phi_max = max(phiMax_vec);

phiMin_vec = [min(min(Winstore_phi_fa_LSF_1_LE_0)),...
    min(min(Winstore_phi_fs_LSF_1_LE_0)),...
    min(min(Winstore_phi_ua_LSF_1_LE_0)),...
    min(min(Winstore_phi_us_LSF_1_LE_0)),...
    min(min(Winstore_phi_uaW_LSF_1_LE_0)),...
    min(min(Winstore_phi_usW_LSF_1_LE_0))];
phi_min = min(phiMin_vec);

%beta
betaMax_vec = [max(max(Winstore_beta_fa_LSF_1_LE_0)),...
    max(max(Winstore_beta_fs_LSF_1_LE_0)),...
    max(max(Winstore_beta_ua_LSF_1_LE_0)),...
    max(max(Winstore_beta_us_LSF_1_LE_0)),...
    max(max(Winstore_beta_uaW_LSF_1_LE_0)),...
    max(max(Winstore_beta_usW_LSF_1_LE_0))];
beta_max = max(betaMax_vec);

betaMin_vec = [min(min(Winstore_beta_fa_LSF_1_LE_0)),...
    min(min(Winstore_beta_fs_LSF_1_LE_0)),...
    min(min(Winstore_beta_ua_LSF_1_LE_0)),...
    min(min(Winstore_beta_us_LSF_1_LE_0)),...
    min(min(Winstore_beta_uaW_LSF_1_LE_0)),...
    min(min(Winstore_beta_usW_LSF_1_LE_0))];
beta_min = min(betaMin_vec);

%% Set up internal values for Figure 1

abdoflip(1:numel(Tstore)) = 2*pi;
livebetarange(1:numel(Tstore)) = 20*(pi/180); 

stop

%% State variables w.r.t. time - individfally
fig1 = figure(1);

%x - fa LSF_1_LE_0
subplot(5,4,1) %(5,6,1)
if NumOf_fa_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),Winstore_x_fa_LSF_1_LE_0,'LineWidth',2);
end
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0'});
ylabel({'x';'(cm)'});
% axis([0, Tstore(end-100+1), x_min, x_max])
axis([0, Tstore(end-100+1), -1, 1])

%x - fs LSF_1_LE_0
subplot(5,4,2)
if NumOf_fs_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),Winstore_x_fs_LSF_1_LE_0,'LineWidth',2);
end
title({'fs';['LSF: ',num2str(LSF_1)];'LE: 0'});
% axis([0, Tstore(end-100+1), x_min, x_max])
axis([0, Tstore(end-100+1), -1, 1])

% %x - ua LSF_1_LE_0
% subplot(5,4,3)
% if NumOf_ua_LSF_1_LE_0_Files > 0
%     plot(Tstore(1:(end-100+1)),Winstore_x_ua_LSF_1_LE_0,'LineWidth',2);
% end
% title({'ua';['LSF: ',num2str(LSF_1)];'LE: 0'});
% % axis([0, Tstore(end-100+1), x_min, x_max])
% axis([0, Tstore(end-100+1), -inf, inf])
% 
% %x - us LSF_1_LE_0
% subplot(5,4,4)
% if NumOf_us_LSF_1_LE_0_Files > 0
%     plot(Tstore(1:(end-100+1)),Winstore_x_us_LSF_1_LE_0,'LineWidth',2);
% end
% title({'us';['LSF: ',num2str(LSF_1)];'LE: 0'});
% % axis([0, Tstore(end-100+1), x_min, x_max])
% axis([0, Tstore(end-100+1), -inf, inf])

%x - uaW LSF_1_LE_0
subplot(5,4,3)
if NumOf_uaW_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),Winstore_x_uaW_LSF_1_LE_0,'LineWidth',2);
end
title({'uaW';['LSF: ',num2str(LSF_1)];'LE: 0'});
% axis([0, Tstore(end-100+1), x_min, x_max])
axis([0, Tstore(end-100+1), -inf, inf])

%x - usW LSF_1_LE_0
subplot(5,4,4)
if NumOf_usW_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),Winstore_x_usW_LSF_1_LE_0,'LineWidth',2);
end
title({'usW';['LSF: ',num2str(LSF_1)];'LE: 0'});
% axis([0, Tstore(end-100+1), x_min, x_max])
axis([0, Tstore(end-100+1), -inf, inf])

%y - fa LSF_1_LE_0
subplot(5,4,5)
if NumOf_fa_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),Winstore_y_fa_LSF_1_LE_0,'LineWidth',2);
    hold on;
end
ylabel({'y';'(cm)'});
plot(Tstore(1:(end-100+1)), y_g,'k--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])
axis([0, Tstore(end-100+1), -11, 11])

%y - fs LSF_1_LE_0
subplot(5,4,6)
if NumOf_fs_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),Winstore_y_fs_LSF_1_LE_0,'LineWidth',2);
    hold on;
end
plot(Tstore(1:(end-100+1)), y_g,'k--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])
axis([0, Tstore(end-100+1), -11, 11])

% %y - ua LSF_1_LE_0
% subplot(5,4,9)
% if NumOf_ua_LSF_1_LE_0_Files > 0
%     plot(Tstore(1:(end-100+1)),Winstore_y_ua_LSF_1_LE_0,'LineWidth',2);
%     hold on;
% end
% plot(Tstore(1:(end-100+1)), y_g,'k--','LineWidth',1);
% % axis([0, Tstore(end-100+1), y_min, y_max])
% axis([0, Tstore(end-100+1), -inf, inf])

% %y - us LSF_1_LE_0
% subplot(5,4,10)
% if NumOf_us_LSF_1_LE_0_Files > 0
%     plot(Tstore(1:(end-100+1)),Winstore_y_us_LSF_1_LE_0,'LineWidth',2);
%     hold on;
% end
% plot(Tstore(1:(end-100+1)), y_g,'k--','LineWidth',1);
% % axis([0, Tstore(end-100+1), y_min, y_max])
% axis([0, Tstore(end-100+1), -inf, inf])

%y - uaW LSF_1_LE_0
subplot(5,4,7)
if NumOf_uaW_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),Winstore_y_uaW_LSF_1_LE_0,'LineWidth',2);
    hold on;
end
plot(Tstore(1:(end-100+1)), y_g,'k--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])
axis([0, Tstore(end-100+1), -inf, inf])

%y - usW LSF_1_LE_0
subplot(5,4,8)
if NumOf_usW_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),Winstore_y_usW_LSF_1_LE_0,'LineWidth',2);
    hold on;
end
plot(Tstore(1:(end-100+1)), y_g,'k--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])
axis([0, Tstore(end-100+1), -inf, inf])


%theta - fa LSF_1_LE_0
subplot(5,4,9)
if NumOf_fa_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_theta_fa_LSF_1_LE_0,'LineWidth',2);
end
ylabel({'theta';'(degrees)'});
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
axis([0, Tstore(end-100+1), 25, 60])

%theta - fs LSF_1_LE_0
subplot(5,4,10)
if NumOf_fs_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_theta_fs_LSF_1_LE_0,'LineWidth',2);
end
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
axis([0, Tstore(end-100+1), 25, 60])

% %theta - ua LSF_1_LE_0
% subplot(5,4,15)
% if NumOf_ua_LSF_1_LE_0_Files > 0
%     plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_theta_ua_LSF_1_LE_0,'LineWidth',2);
% end
% % axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
% axis([0, Tstore(end-100+1), -inf, inf])
% 
% %theta - us LSF_1_LE_0
% subplot(5,4,16)
% if NumOf_us_LSF_1_LE_0_Files > 0
%     plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_theta_us_LSF_1_LE_0,'LineWidth',2);
% end
% % axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
% axis([0, Tstore(end-100+1), -inf, inf])

%theta - uaW LSF_1_LE_0
subplot(5,4,11)
if NumOf_uaW_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_theta_uaW_LSF_1_LE_0,'LineWidth',2);
end
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
axis([0, Tstore(end-100+1), -inf, inf])

%theta - usW LSF_1_LE_0
subplot(5,4,12)
if NumOf_usW_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_theta_usW_LSF_1_LE_0,'LineWidth',2);
end
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
axis([0, Tstore(end-100+1), -inf, inf])

%phi - fa LSF_1_LE_0
subplot(5,4,13)
if NumOf_fa_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_phi_fa_LSF_1_LE_0,'LineWidth',2);
    hold on;
end
% plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
ylabel({'phi';'(degrees)'});
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])
axis([0, Tstore(end-100+1), 190, 280])

%phi - fs LSF_1_LE_0
subplot(5,4,14)
if NumOf_fs_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_phi_fs_LSF_1_LE_0,'LineWidth',2);
    hold on;
end
% plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])
axis([0, Tstore(end-100+1), 190, 280])

% %phi - ua LSF_1_LE_0
% subplot(5,4,21)
% if NumOf_ua_LSF_1_LE_0_Files > 0
%     plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_phi_ua_LSF_1_LE_0,'LineWidth',2);
%     hold on;
% end
% % plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% % axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])
% axis([0, Tstore(end-100+1), -inf, inf])

% %phi - us LSF_1_LE_0
% subplot(5,4,22)
% if NumOf_us_LSF_1_LE_0_Files > 0
%     plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_phi_us_LSF_1_LE_0,'LineWidth',2);
%     hold on;
% end
% % plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% % axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])
% axis([0, Tstore(end-100+1), -inf, inf])

%phi - uaW LSF_1_LE_0
subplot(5,4,15)
if NumOf_uaW_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_phi_uaW_LSF_1_LE_0,'LineWidth',2);
    hold on;
end
% plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])
axis([0, Tstore(end-100+1), -inf, inf])

%phi - usW LSF_1_LE_0
subplot(5,4,16)
if NumOf_usW_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_phi_usW_LSF_1_LE_0,'LineWidth',2);
    hold on;
end
% plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])
axis([0, Tstore(end-100+1), -inf, inf])

%beta - fa LSF_1_LE_0
subplot(5,4,17)
if NumOf_fa_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_beta_fa_LSF_1_LE_0,'LineWidth',2);
    hold on;
end
xlabel({'Time (sec),'; ['(n = ',num2str(NumOf_fa_LSF_1_LE_0_Files),')']}); 
ylabel({'beta';'(degrees)'});
% plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% hold on;
% plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',0:40:80);
axis([0, Tstore(end-100+1), -20, 45])

%beta - fs LSF_1_LE_0
subplot(5,4,18)
if NumOf_fs_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_beta_fs_LSF_1_LE_0,'LineWidth',2);
    hold on;
end
xlabel({'Time (sec),'; ['(n = ',num2str(NumOf_fs_LSF_1_LE_0_Files),')']});
% plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% hold on;
% plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',0:40:80);
axis([0, Tstore(end-100+1), -20, 45])

% %beta - ua LSF_1_LE_0
% subplot(5,4,27)
% if NumOf_ua_LSF_1_LE_0_Files > 0
%     plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_beta_ua_LSF_1_LE_0,'LineWidth',2);
%     hold on;
% end
% xlabel({'Time (sec),'; ['(n = ',num2str(NumOf_ua_LSF_1_LE_0_Files),')']}); 
% % plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% % hold on;
% % plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% % axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% % set(gca,'ytick',0:40:80);
% axis([0, Tstore(end-100+1), -inf, inf])

% %beta - us LSF_1_LE_0
% subplot(5,4,28)
% if NumOf_us_LSF_1_LE_0_Files > 0
%     plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_beta_us_LSF_1_LE_0,'LineWidth',2);
%     hold on;
% end
% xlabel({'Time (sec),'; ['(n = ',num2str(NumOf_us_LSF_1_LE_0_Files),')']}); 
% % plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% % hold on;
% % plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% % axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% % % set(gca,'ytick',-120:40:180);
% axis([0, Tstore(end-100+1), -inf, inf])

%beta - uaW LSF_1_LE_0
subplot(5,4,19)
if NumOf_uaW_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_beta_uaW_LSF_1_LE_0,'LineWidth',2);
    hold on;
end
xlabel({'Time (sec),'; ['(n = ',num2str(NumOf_uaW_LSF_1_LE_0_Files),')']}); 
% plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% hold on;
% plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',0:40:80);
axis([0, Tstore(end-100+1), -inf, inf])

%beta - usW LSF_1_LE_0
subplot(5,4,20)
if NumOf_usW_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_beta_usW_LSF_1_LE_0,'LineWidth',2);
    hold on;
end
xlabel({'Time (sec),'; ['(n = ',num2str(NumOf_usW_LSF_1_LE_0_Files),')']}); 
% plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% hold on;
% plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% % set(gca,'ytick',-120:40:180);
axis([0, Tstore(end-100+1), -inf, inf])

set(fig1, 'Position', [0,0, 1250, 2500])

saveas(gcf,'../Figures_MPC/compTmt/LSF_1/Fig1_compTmt_stateVarsIndiv_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/compTmt/LSF_1/Fig1_compTmt_stateVarsIndiv_LSF_1_v12h.jpg');

disp('Figure 1 done')

%% Troubleshoot beta
% the cases where the abdomen spins over 500 degrees (unrealistic)


temp1 = (180/pi).*Winstore_beta_ua_LSF_1_LE_0;
f1 = temp1>500;
id1 = f1(:,25000:25002);

temp2 = (180/pi).*Winstore_beta_us_LSF_1_LE_0;
f2 = temp2>500;
id2 = f2(:,25000:25002);

f2b = temp2<0;
id2b = f2b(:,25000:25002);

%Note, the above filters helped weed out a majority of the problematic
%cases, but not all. The rest I had to identify at different time points
%manually.

%The problematic cases for "ua" are as follows:
    %1,3,6,7,8,9,10,11,15,17,18,19,22,25,28,31,37,38,39

%The problematic cases for "us" are as follows:
    %2,3,4,5,6,9,12,13,14,15,16,17,18,19,20,21,22,23,25,28,29,31,32,33,35,36,37,39

%These cases will be manually removed and the good ones will be renamed.

fig3 = figure(3);

%beta - ua LSF_1_LE_0
subplot(2,1,1)
if NumOf_ua_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_beta_ua_LSF_1_LE_0,'LineWidth',2);
    hold on;
end
ylabel('Beta (degrees), ua treatment')
xlabel({'Time (sec),'; ['(n = ',num2str(NumOf_ua_LSF_1_LE_0_Files-19),')']}); 
% plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% hold on;
% plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',0:40:80);
axis([0, Tstore(end-100+1), -inf, inf])

%beta - us LSF_1_LE_0
subplot(2,1,2)
if NumOf_us_LSF_1_LE_0_Files > 0
    plot(Tstore(1:(end-100+1)),(180/pi).*Winstore_beta_us_LSF_1_LE_0,'LineWidth',2);
    hold on;
end
ylabel('Beta (degrees), us treatment')
xlabel({'Time (sec),'; ['(n = ',num2str(NumOf_us_LSF_1_LE_0_Files-28),')']}); 
% plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% hold on;
% plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% % set(gca,'ytick',-120:40:180);
axis([0, Tstore(end-100+1), -inf, inf])

stop

%% State variables w.r.t. time - mean
fig2 = figure(2);

%x - fa LSF_1_LE_0
subplot(5,4,1)
Fig2_1 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_fa_LSF_1_LE_0,...
    std_x_fa_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_1.mainLine.LineWidth = 2;
% title({'Fully actuated';['LSF: ',num2str(LSF_1)];'LE: 0'})
title('Fully actuated')
ylabel({'x';'(cm)'});
% axis([0, Tstore(end-100+1), x_min, x_max])
axis([0, Tstore(end-100+1), -1, 1])

%x - fs LSF_1_LE_0
subplot(5,4,2)
Fig2_2 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_fs_LSF_1_LE_0,...
    std_x_fs_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_2.mainLine.LineWidth = 2;
% title({'fs';['LSF: ',num2str(LSF_1)];'LE: 0'})
title({'Fully actuated'; 'and location of';'applied force shifted'})
% axis([0, Tstore(end-100+1), x_min, x_max])
axis([0, Tstore(end-100+1), -1, 1])

% %x - ua LSF_1_LE_0
% subplot(5,4,3)
% Fig2_5 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_ua_LSF_1_LE_0,...
%     std_x_ua_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% Fig2_5.mainLine.LineWidth = 2;
% % title({'ua';['LSF: ',num2str(LSF_1)];'LE: 0'})
% title('Underactuated')
% % axis([0, Tstore(end-100+1), x_min, x_max])
% axis([0, Tstore(end-100+1), -35, 35])
% 
% %x - us LSF_1_LE_0
% subplot(5,4,4)
% Fig2_6 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_us_LSF_1_LE_0,...
%     std_x_us_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% Fig2_6.mainLine.LineWidth = 2;
% % title({'us';['LSF: ',num2str(LSF_1)];'LE: 0'})
% title({'Underactuated'; 'and location of';'applied force shifted'})
% % axis([0, Tstore(end-100+1), x_min, x_max])
% axis([0, Tstore(end-100+1), -35, 35])

%x - uaW LSF_1_LE_0
subplot(5,4,3)
Fig2_3 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_uaW_LSF_1_LE_0,...
    std_x_uaW_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_3.mainLine.LineWidth = 2;
% title({'ua';['LSF: ',num2str(LSF_1)];'LE: 0'})
title({'Underactuated'; '(physically constrained rotations)'})
% axis([0, Tstore(end-100+1), x_min, x_max])
axis([0, Tstore(end-100+1), -5, 5])

%x - usW LSF_1_LE_0
subplot(5,4,4)
Fig2_4 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_usW_LSF_1_LE_0,...
    std_x_usW_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_4.mainLine.LineWidth = 2;
% title({'us';['LSF: ',num2str(LSF_1)];'LE: 0'})
title({'Underactuated'; 'and location of';'applied force shifted';...
    '(physically constrained rotations)'})
% axis([0, Tstore(end-100+1), x_min, x_max])
axis([0, Tstore(end-100+1), -5, 5])

%y - fa LSF_1_LE_0
subplot(5,4,5)
Fig2_5 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_fa_LSF_1_LE_0,...
    std_y_fa_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_5.mainLine.LineWidth = 2; 
hold on;
ylabel({'y';'(cm)'});
plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])
axis([0, Tstore(end-100+1), -11, 11])

%y - fs LSF_1_LE_0
subplot(5,4,6)
Fig2_6 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_fs_LSF_1_LE_0,...
    std_y_fs_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_6.mainLine.LineWidth = 2; 
hold on;
plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])
axis([0, Tstore(end-100+1), -11, 11])

% %y - ua LSF_1_LE_0
% subplot(5,4,9)
% Fig2_9 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_ua_LSF_1_LE_0,...
%     std_y_ua_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% Fig2_9.mainLine.LineWidth = 2;
% hold on;
% plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% % axis([0, Tstore(end-100+1), y_min, y_max])
% axis([0, Tstore(end-100+1), -11, 11])
% 
% %y - us LSF_1_LE_0
% subplot(5,4,10)
% Fig2_10 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_us_LSF_1_LE_0,...
%     std_y_us_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% Fig2_10.mainLine.LineWidth = 2; 
% hold on;
% plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% % axis([0, Tstore(end-100+1), y_min, y_max])
% axis([0, Tstore(end-100+1), -11, 11])

%y - uaW LSF_1_LE_0
subplot(5,4,7)
Fig2_7 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_uaW_LSF_1_LE_0,...
    std_y_uaW_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_7.mainLine.LineWidth = 2;
hold on;
plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])
axis([0, Tstore(end-100+1), -11, 11])

%y - usW LSF_1_LE_0
subplot(5,4,8)
Fig2_8 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_usW_LSF_1_LE_0,...
    std_y_usW_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_8.mainLine.LineWidth = 2; 
hold on;
plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% axis([0, Tstore(end-100+1), y_min, y_max])
axis([0, Tstore(end-100+1), -11, 11])


%theta - fa LSF_1_LE_0
subplot(5,4,9)
Fig2_9 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_fa_LSF_1_LE_0,...
    (180/pi)*std_theta_fa_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
% ylabel({'theta';'(degrees)'});
ylabel({'Head-thorax';'motion';'(degrees)'});
Fig2_9.mainLine.LineWidth = 2;
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
axis([0, Tstore(end-100+1), 35, 50])

%theta - fs LSF_1_LE_0
subplot(5,4,10)
Fig2_10 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_fs_LSF_1_LE_0,...
    (180/pi)*std_theta_fs_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_10.mainLine.LineWidth = 2;
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
axis([0, Tstore(end-100+1), 35, 50])

% %theta - ua LSF_1_LE_0
% subplot(5,4,15)
% Fig2_15 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_ua_LSF_1_LE_0,...
%     (180/pi)*std_theta_ua_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% Fig2_15.mainLine.LineWidth = 2;
% % axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
% axis([0, Tstore(end-100+1), -110, 230])
% 
% %theta - us LSF_1_LE_0
% subplot(5,4,16)
% Fig2_16 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_us_LSF_1_LE_0,...
%     (180/pi)*std_theta_us_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% Fig2_16.mainLine.LineWidth = 2;
% % axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
% axis([0, Tstore(end-100+1), -110, 230])

%theta - uaW LSF_1_LE_0
subplot(5,4,11)
Fig2_11 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_uaW_LSF_1_LE_0,...
    (180/pi)*std_theta_uaW_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_11.mainLine.LineWidth = 2;
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
axis([0, Tstore(end-100+1), -100, 175])

%theta - usW LSF_1_LE_0
subplot(5,4,12)
Fig2_12 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_usW_LSF_1_LE_0,...
    (180/pi)*std_theta_usW_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_12.mainLine.LineWidth = 2;
% axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
axis([0, Tstore(end-100+1), -100, 175])

%phi - fa LSF_1_LE_0
subplot(5,4,13)
Fig2_13 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_phi_fa_LSF_1_LE_0,...
    (180/pi)*std_phi_fa_LSF_1_LE_0,...
    'lineprops', '-k','transparent',false,'patchSaturation',0.075);
Fig2_13.mainLine.LineWidth = 2;
% hold on;
% plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% ylabel({'phi';'(degrees)'});
ylabel({'Abdomen motion';'(degrees)'});
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])
axis([0, Tstore(end-100+1), 210, 255])

%phi - fs LSF_1_LE_0
subplot(5,4,14)
Fig2_14 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_phi_fs_LSF_1_LE_0,...
    (180/pi)*std_phi_fs_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_14.mainLine.LineWidth = 2;
% hold on;
% plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])
axis([0, Tstore(end-100+1), 210, 255])

% %phi - ua LSF_1_LE_0
% subplot(5,4,21)
% Fig2_21 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_phi_ua_LSF_1_LE_0,...
%     (180/pi)*std_phi_ua_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% Fig2_21.mainLine.LineWidth = 2;
% % hold on;
% % plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% % axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])
% axis([0, Tstore(end-100+1), 120, 1550])
% 
% %phi - us LSF_1_LE_0
% subplot(5,4,22)
% Fig2_22 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_phi_us_LSF_1_LE_0,...
%     (180/pi)*std_phi_us_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);   
% Fig2_22.mainLine.LineWidth = 2;
% % hold on;
% % plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% % axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])
% axis([0, Tstore(end-100+1), 120, 1550])

%phi - uaW LSF_1_LE_0
subplot(5,4,15)
Fig2_15 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_phi_uaW_LSF_1_LE_0,...
    (180/pi)*std_phi_uaW_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_15.mainLine.LineWidth = 2;
% hold on;
% plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])
axis([0, Tstore(end-100+1), 130, 750])

%phi - usW LSF_1_LE_0
subplot(5,4,16)
Fig2_16 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_phi_usW_LSF_1_LE_0,...
    (180/pi)*std_phi_usW_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);   
Fig2_16.mainLine.LineWidth = 2;
% hold on;
% plot(Tstore, (180/pi)*abdoflip,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*phi_min, (180/pi)*phi_max])
axis([0, Tstore(end-100+1), 130, 750])

%beta - fa LSF_1_LE_0
subplot(5,4,17)
Fig2_17 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_fa_LSF_1_LE_0,...
    (180/pi)*std_beta_fa_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_17.mainLine.LineWidth = 2;
% xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOf_fa_LSF_1_LE_0_Files),')']}); 
xlabel({'Time (sec),';['(n = ',num2str(NumOf_fa_LSF_1_LE_0_Files),')']}); 
% ylabel({'beta';'(degrees)'});
ylabel({'Abdominal flexion';'(degrees)'});
% hold on;
% plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% hold on;
% plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',0:40:80);
axis([0, Tstore(end-100+1), -5, 25])

%beta - fs LSF_1_LE_0
subplot(5,4,18)
Fig2_18 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_fs_LSF_1_LE_0,...
    (180/pi)*std_beta_fs_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_18.mainLine.LineWidth = 2;
% xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOf_fs_LSF_1_LE_0_Files),')']}); 
xlabel({'Time (sec),'; ['(n = ',num2str(NumOf_fs_LSF_1_LE_0_Files),')']}); 
% hold on;
% plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% hold on;
% plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',0:40:80);
axis([0, Tstore(end-100+1), -5, 25])

% %beta - ua LSF_1_LE_0
% subplot(5,4,27)
% Fig2_27 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_ua_LSF_1_LE_0,...
%     (180/pi)*std_beta_ua_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);   
% Fig2_27.mainLine.LineWidth = 2;
% % xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOf_ua_LSF_1_LE_0_Files),')']}); 
% xlabel({'Time (sec),';['(n = ',num2str(NumOf_ua_LSF_1_LE_0_Files),')']}); 
% % hold on;
% % plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% % hold on;
% % plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% % axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% % set(gca,'ytick',0:40:80);
% axis([0, Tstore(end-100+1), -15, 1200])

% %beta - us LSF_1_LE_0
% subplot(5,4,28)
% Fig2_28 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_us_LSF_1_LE_0,...
%     (180/pi)*std_beta_us_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% Fig2_28.mainLine.LineWidth = 2;
% % xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOf_us_LSF_1_LE_0_Files),')']}); 
% xlabel({'Time (sec),'; ['(n = ',num2str(NumOf_us_LSF_1_LE_0_Files),')']}); 
% % hold on;
% % plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% % hold on;
% % plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% % axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% % set(gca,'ytick',-120:40:180);
% axis([0, Tstore(end-100+1), -15, 1200])

%beta - uaW LSF_1_LE_0
subplot(5,4,19)
Fig2_19 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_uaW_LSF_1_LE_0,...
    (180/pi)*std_beta_uaW_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);   
Fig2_19.mainLine.LineWidth = 2;
% xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOf_ua_LSF_1_LE_0_Files),')']}); 
xlabel({'Time (sec),';['(n = ',num2str(NumOf_uaW_LSF_1_LE_0_Files),')']}); 
% hold on;
% plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% hold on;
% plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',0:40:80);
axis([0, Tstore(end-100+1), 0, 500])

%beta - usW LSF_1_LE_0
subplot(5,4,20)
Fig2_20 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_usW_LSF_1_LE_0,...
    (180/pi)*std_beta_usW_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
Fig2_20.mainLine.LineWidth = 2;
% xlabel({'Time (sec),'; 'con'; ['(n = ',num2str(NumOf_us_LSF_1_LE_0_Files),')']}); 
xlabel({'Time (sec),'; ['(n = ',num2str(NumOf_usW_LSF_1_LE_0_Files),')']}); 
% hold on;
% plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% hold on;
% plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% set(gca,'ytick',-120:40:180);
axis([0, Tstore(end-100+1), 0, 500])

set(fig2, 'Position', [0,0, 1250, 2500])

saveas(gcf,'../Figures_MPC/compTmt/LSF_1/Fig2_compTmt_stateVarsMean_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/compTmt/LSF_1/Fig2_compTmt_stateVarsMean_LSF_1_v12h.jpg');

disp('Figure 2 done')

%% Figure for SACNAS 2019
% 
% figSACNAS = figure(3);
% 
% %x motion - expected
% subplot(4,5,1)
% plot(Tstore(1:(end-100+1)), x_g,'r-','LineWidth',2);
% ylabel({'x';'(cm)'});
% axis([0, Tstore(end-100+1), -2, 2])
% 
% %x - fs LSF_1_LE_0
% subplot(4,5,2)
% figSACNAS_2 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_fs_LSF_1_LE_0,...
%     std_x_fs_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% figSACNAS_2.mainLine.LineWidth = 2;
% % title({'fs';['LSF: ',num2str(LSF_1)];'LE: 0'})
% title('fs')
% % axis([0, Tstore(end-100+1), x_min, x_max])
% axis([0, Tstore(end-100+1), -2, 2])
% 
% %x - fa LSF_1_LE_0
% subplot(4,5,3)
% figSACNAS_3 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_fa_LSF_1_LE_0,...
%     std_x_fa_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% figSACNAS_3.mainLine.LineWidth = 2;
% % title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0'})
% title('fa')
% % ylabel({'x';'(cm)'});
% % axis([0, Tstore(end-100+1), x_min, x_max])
% axis([0, Tstore(end-100+1), -2, 2])
% 
% %x - us LSF_1_LE_0
% subplot(4,5,4)
% figSACNAS_4 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_us_LSF_1_LE_0,...
%     std_x_us_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% figSACNAS_4.mainLine.LineWidth = 2;
% % title({'us';['LSF: ',num2str(LSF_1)];'LE: 0'})
% title('us')
% % axis([0, Tstore(end-100+1), x_min, x_max])
% axis([0, Tstore(end-100+1), -2, 2])
% 
% %x - ua LSF_1_LE_0
% subplot(4,5,5)
% figSACNAS_5 = shadedErrorBar(Tstore(1:(end-100+1)),mean_x_ua_LSF_1_LE_0,...
%     std_x_ua_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% figSACNAS_5.mainLine.LineWidth = 2;
% % title({'ua';['LSF: ',num2str(LSF_1)];'LE: 0'})
% title('ua')
% % axis([0, Tstore(end-100+1), x_min, x_max])
% % axis([0, Tstore(end-100+1), -2, 2])
% 
% %y-motion expected
% subplot(4,5,6)
% plot(Tstore(1:(end-100+1)), y_g,'r-','LineWidth',2);
% ylabel({'y';'(cm)'});
% axis([0, Tstore(end-100+1), -11, 11])
% 
% %y - fs LSF_1_LE_0
% subplot(4,5,7)
% figSACNAS_7 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_fs_LSF_1_LE_0,...
%     std_y_fs_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% figSACNAS_7.mainLine.LineWidth = 2;  %was 1.2
% % hold on;
% % plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% % axis([0, Tstore(end-100+1), y_min, y_max])
% axis([0, Tstore(end-100+1), -11, 11])
% 
% %y - fa LSF_1_LE_0
% subplot(4,5,8)
% figSACNAS_8 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_fa_LSF_1_LE_0,...
%     std_y_fa_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% figSACNAS_8.mainLine.LineWidth = 2; 
% % hold on;
% % ylabel({'y';'(cm)'});
% % plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% % axis([0, Tstore(end-100+1), y_min, y_max])
% axis([0, Tstore(end-100+1), -11, 11])
% 
% %y - us LSF_1_LE_0
% subplot(4,5,9)
% figSACNAS_9 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_us_LSF_1_LE_0,...
%     std_y_us_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% figSACNAS_9.mainLine.LineWidth = 2; 
% % hold on;
% % plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% % axis([0, Tstore(end-100+1), y_min, y_max])
% axis([0, Tstore(end-100+1), -11, 11])
% 
% %y - ua LSF_1_LE_0
% subplot(4,5,10)
% figSACNAS_10 = shadedErrorBar(Tstore(1:(end-100+1)),mean_y_ua_LSF_1_LE_0,...
%     std_y_ua_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% figSACNAS_10.mainLine.LineWidth = 2;
% % hold on;
% % plot(Tstore(1:(end-100+1)), y_g,'r--','LineWidth',1);
% % axis([0, Tstore(end-100+1), y_min, y_max])
% % axis([0, Tstore(end-100+1), -11, 11])
% 
% %theta expected motion
% subplot(4,5,11)
% plot(Tstore(1:(end-100+1)), (180/pi).*theta_g,'r-','LineWidth',2);
% ylabel({'Head motion';'(degrees)'});
% axis([0, Tstore(end-100+1), 20, 70])
% 
% %theta - fs LSF_1_LE_0
% subplot(4,5,12)
% figSACNAS_12 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_fs_LSF_1_LE_0,...
%     (180/pi)*std_theta_fs_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% figSACNAS_12.mainLine.LineWidth = 2;
% % axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
% axis([0, Tstore(end-100+1), 20, 70])
% 
% %theta - fa LSF_1_LE_0
% subplot(4,5,13)
% figSACNAS_13 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_fa_LSF_1_LE_0,...
%     (180/pi)*std_theta_fa_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% % ylabel({'theta';'(degrees)'});
% figSACNAS_13.mainLine.LineWidth = 2;
% % axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
% axis([0, Tstore(end-100+1), 20, 70])
% 
% %theta - us LSF_1_LE_0
% subplot(4,5,14)
% figSACNAS_14 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_us_LSF_1_LE_0,...
%     (180/pi)*std_theta_us_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% figSACNAS_14.mainLine.LineWidth = 2;
% % axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
% axis([0, Tstore(end-100+1), 20, 70])
% 
% %theta - ua LSF_1_LE_0
% subplot(4,5,15)
% figSACNAS_15 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_theta_ua_LSF_1_LE_0,...
%     (180/pi)*std_theta_ua_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% figSACNAS_15.mainLine.LineWidth = 2;
% % axis([0, Tstore(end-100+1), (180/pi)*theta_min, (180/pi)*theta_max])
% % axis([0, Tstore(end-100+1), 20, 70])
% 
% %beta expected motion
% subplot(4,5,16)
% % plot(Tstore(1:(end-100+1)), y_g,'r-','LineWidth',2);
% % hold on;
% plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% hold on;
% plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% ylabel({'Abdominal flexion';'(degrees)'});
% axis([0, Tstore(end-100+1), -150, 150])
% 
% %beta - fs LSF_1_LE_0
% subplot(4,5,17)
% figSACNAS_17 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_fs_LSF_1_LE_0,...
%     (180/pi)*std_beta_fs_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% figSACNAS_17.mainLine.LineWidth = 2;
% xlabel({'Time (sec),'; ['(n = ',num2str(NumOf_fs_LSF_1_LE_0_Files),')']}); 
% hold on;
% plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% hold on;
% plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% % axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% % set(gca,'ytick',0:40:80);
% axis([0, Tstore(end-100+1), -150, 150])
% 
% %beta - fa LSF_1_LE_0
% subplot(4,5,18)
% figSACNAS_18 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_fa_LSF_1_LE_0,...
%     (180/pi)*std_beta_fa_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% figSACNAS_18.mainLine.LineWidth = 2;
% xlabel({'Time (sec),'; ['(n = ',num2str(NumOf_fa_LSF_1_LE_0_Files),')']}); 
% % ylabel({'beta';'(degrees)'});
% hold on;
% plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% hold on;
% plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% % axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% % set(gca,'ytick',0:40:80);
% axis([0, Tstore(end-100+1), -150, 150])
% 
% %beta - us LSF_1_LE_0
% subplot(4,5,19)
% figSACNAS_19 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_us_LSF_1_LE_0,...
%     (180/pi)*std_beta_us_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);
% figSACNAS_19.mainLine.LineWidth = 2;
% xlabel({'Time (sec),'; ['(n = ',num2str(NumOf_us_LSF_1_LE_0_Files),')']}); 
% hold on;
% plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% hold on;
% plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% % axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% % set(gca,'ytick',-120:40:180);
% axis([0, Tstore(end-100+1), -150, 150])
% 
% %beta - ua LSF_1_LE_0
% subplot(4,5,20)
% figSACNAS_20 = shadedErrorBar(Tstore(1:(end-100+1)),(180/pi)*mean_beta_ua_LSF_1_LE_0,...
%     (180/pi)*std_beta_ua_LSF_1_LE_0,...
%     'lineprops','-k','transparent',false,'patchSaturation',0.075);   
% figSACNAS_20.mainLine.LineWidth = 2;
% xlabel({'Time (sec),'; ['(n = ',num2str(NumOf_ua_LSF_1_LE_0_Files),')']}); 
% hold on;
% plot(Tstore, (180/pi)*livebetarange,'r--','LineWidth',1)
% hold on;
% plot(Tstore, -(180/pi)*livebetarange,'r--','LineWidth',1)
% % axis([0, Tstore(end-100+1), (180/pi)*beta_min, (180/pi)*beta_max])
% % set(gca,'ytick',0:40:80);
% % axis([0, Tstore(end-100+1), -150, 150])
% 
% set(figSACNAS, 'Position', [0,0, 1250, 2500])
% 
% saveas(gcf,'../Figures_MPC/compTmt/LSF_1/SACNAS2019_compTmt_stateVarsMean_LSF_1_v12h.fig');
% saveas(gcf,'../Figures_MPC/compTmt/LSF_1/SACNAS2019_compTmt_stateVarsMean_LSF_1_v12h.jpg');
% 
% disp('SACNAS 2019 dynamics figure done')

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