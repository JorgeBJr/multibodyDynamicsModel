%2019/7/19 Script developed to plot the work done by the applied efforts
%(i.e. abdominal torque, wing torque, aerodynamic force)
%2020/5/19 Script modified for inertial modification models (used for paper
    %2).
%2021/5/20 Script modified for legnth scale factor (LSF), and length
%extension (LE) modifications. Also used for paper 2. 

%% We the People, in Order to form a perfect Union,...
clc;
clearvars;
close all;
disp(['Time at start: ',datestr(now)])

%% Optional save prompt
prompt1 = 'Do you want to save the stats CSVs? y/n: ';
str1 = input(prompt1, 's');

%% List the directory stuff

%ALL SIMS FOR THIS PAPER ARE ONLY FULLY-ACTUATED (fa) TREATMENT!
listOf_LSF_1_LE_0_Files = dir('../SimData_MPC/fa_Sims/LSF_1/Winstore_*_fa_LSF_1_LE_0.mat');
listOf_LSF_1_LE_p2_Files = dir('../SimData_MPC/fa_Sims/LSF_1/Winstore_*_fa_LSF_1_LE_p2.mat'); 
listOf_LSF_1_LE_p4_Files = dir('../SimData_MPC/fa_Sims/LSF_1/Winstore_*_fa_LSF_1_LE_p4.mat'); 
listOf_LSF_p5_LE_0_Files = dir('../SimData_MPC/fa_Sims/LSF_p5/Winstore_*_fa_LSF_p5_LE_0.mat'); 
listOf_LSF_p5_LE_p2_Files = dir('../SimData_MPC/fa_Sims/LSF_p5/Winstore_*_fa_LSF_p5_LE_p2.mat'); 
listOf_LSF_p5_LE_p4_Files = dir('../SimData_MPC/fa_Sims/LSF_p5/Winstore_*_fa_LSF_p5_LE_p4.mat'); 
listOf_LSF_2_LE_0_Files = dir('../SimData_MPC/fa_Sims/LSF_2/Winstore_*_fa_LSF_2_LE_0.mat');  
listOf_LSF_2_LE_p2_Files = dir('../SimData_MPC/fa_Sims/LSF_2/Winstore_*_fa_LSF_2_LE_p2.mat');
listOf_LSF_2_LE_p4_Files = dir('../SimData_MPC/fa_Sims/LSF_2/Winstore_*_fa_LSF_2_LE_p4.mat');

%Number of files for each treatment
NumOf_LSF_p5_LE_0_Files = numel(listOf_LSF_p5_LE_0_Files);
NumOf_LSF_p5_LE_p2_Files = numel(listOf_LSF_p5_LE_p2_Files);
NumOf_LSF_p5_LE_p4_Files = numel(listOf_LSF_p5_LE_p4_Files);
NumOf_LSF_1_LE_0_Files = numel(listOf_LSF_1_LE_0_Files);
NumOf_LSF_1_LE_p2_Files = numel(listOf_LSF_1_LE_p2_Files);
NumOf_LSF_1_LE_p4_Files = numel(listOf_LSF_1_LE_p4_Files);
NumOf_LSF_2_LE_0_Files = numel(listOf_LSF_2_LE_0_Files);
NumOf_LSF_2_LE_p2_Files = numel(listOf_LSF_2_LE_p2_Files);
NumOf_LSF_2_LE_p4_Files = numel(listOf_LSF_2_LE_p4_Files);

%the asterisk is a wildcard
%The dir function returns a "listing" of an M x 1 "structure." The 
%structure has five fields in this case listing: name, date, byte, isdir, 
%datenum.
%I used the wildcard because I know the number of text files will
%definitely increase as we gather more data.
%For more information enter   help dir   into MATLAB mainframe
disp('Load time vector')
load('../SimData_MPC/Tstore_MPC_hws_sp.mat')

%% Moth morphometric parameters

%Length scale factors (LSF)
LSF_p5 = 0.5; 
LSF_1 = 1; 
LSF_2 = 2; 

%Length extension values (LE)
LE_0 = 0; 
LE_p2 = 0.2; 
LE_p4 = 0.4; 

%The densities of the head and butt do not change regardless of size
rho_head = 0.9; %The density of the head-thorax in g/(cm^3)
rho_butt = 0.4; %The density of the abdomen in g/(cm^3)
g = 980; %g is the acceleration due to gravity in cm/(s^2)

%The following morphometrics will vary by Length Scale Factor (LSF)
%LSF_p5
ahead_LSF_p5 = LSF_p5*0.9; %Major axis of the head-thorax ellipsoid in cm.
abutt_LSF_p5 = LSF_p5*1.9; %Major axis of the abdomen ellipsoid in cm.
bhead_LSF_p5 = LSF_p5*0.5; %Minor axis of the head-thorax ellipsoid in cm.
bbutt_LSF_p5 = LSF_p5*0.75; %Minor axis of the abdomen ellipsoid in cm.
%LSF_1
ahead_LSF_1 = LSF_1*0.9; %Major axis of the head-thorax ellipsoid in cm.
abutt_LSF_1 = LSF_1*1.9; %Major axis of the abdomen ellipsoid in cm.
bhead_LSF_1 = LSF_1*0.5; %Minor axis of the head-thorax ellipsoid in cm.
bbutt_LSF_1 = LSF_1*0.75; %Minor axis of the abdomen ellipsoid in cm.
%LSF_2
ahead_LSF_2 = LSF_2*0.9; %Major axis of the head-thorax ellipsoid in cm.
abutt_LSF_2 = LSF_2*1.9; %Major axis of the abdomen ellipsoid in cm.
bhead_LSF_2 = LSF_2*0.5; %Minor axis of the head-thorax ellipsoid in cm.
bbutt_LSF_2 = LSF_2*0.75; %Minor axis of the abdomen ellipsoid in cm.

%The various length of petiole extension(s) in cm
%LSF_p5
L_pet_LSF_p5_LE_0 = LE_0*(2*(ahead_LSF_p5 + abutt_LSF_p5)); 
L_pet_LSF_p5_LE_p2 = LE_p2*(2*(ahead_LSF_p5 + abutt_LSF_p5)); 
L_pet_LSF_p5_LE_p4 = LE_p4*(2*(ahead_LSF_p5 + abutt_LSF_p5)); 
%LSF_1
L_pet_LSF_1_LE_0 = LE_0*(2*(ahead_LSF_1 + abutt_LSF_1)); %
L_pet_LSF_1_LE_p2 = LE_p2*(2*(ahead_LSF_1 + abutt_LSF_1));
L_pet_LSF_1_LE_p4 = LE_p4*(2*(ahead_LSF_1 + abutt_LSF_1)); 
%LSF_2
L_pet_LSF_2_LE_0 = LE_0*(2*(ahead_LSF_2 + abutt_LSF_2)); 
L_pet_LSF_2_LE_p2 = LE_p2*(2*(ahead_LSF_2 + abutt_LSF_2)); 
L_pet_LSF_2_LE_p4 = LE_p4*(2*(ahead_LSF_2 + abutt_LSF_2)); 

%The various body lengths *including petiole lengths* in cm
%LSF_p5
bl_LSF_p5_LE_0 = (2*(ahead_LSF_p5 + abutt_LSF_p5)) + L_pet_LSF_p5_LE_0;
bl_LSF_p5_LE_p2 = (2*(ahead_LSF_p5 + abutt_LSF_p5)) + L_pet_LSF_p5_LE_p2;
bl_LSF_p5_LE_p4 = (2*(ahead_LSF_p5 + abutt_LSF_p5)) + L_pet_LSF_p5_LE_p4;
%LSF_1
bl_LSF_1_LE_0 = (2*(ahead_LSF_1 + abutt_LSF_1)) + L_pet_LSF_1_LE_0;
bl_LSF_1_LE_p2 = (2*(ahead_LSF_1 + abutt_LSF_1)) + L_pet_LSF_1_LE_p2;
bl_LSF_1_LE_p4 = (2*(ahead_LSF_1 + abutt_LSF_1)) + L_pet_LSF_1_LE_p4;
%LSF_2
bl_LSF_2_LE_0 = (2*(ahead_LSF_2 + abutt_LSF_2)) + L_pet_LSF_2_LE_0;
bl_LSF_2_LE_p2 = (2*(ahead_LSF_2 + abutt_LSF_2)) + L_pet_LSF_2_LE_p2;
bl_LSF_2_LE_p4 = (2*(ahead_LSF_2 + abutt_LSF_2)) + L_pet_LSF_2_LE_p4;

%Mass of the head-thorax ellipsoid (m1), and its variants in grams.
m1_LSF_p5 = rho_head*(4/3)*pi*(bhead_LSF_p5^2)*ahead_LSF_p5; 
m1_LSF_1 = rho_head*(4/3)*pi*(bhead_LSF_1^2)*ahead_LSF_1; 
m1_LSF_2 = rho_head*(4/3)*pi*(bhead_LSF_2^2)*ahead_LSF_2; 

%Mass of the abdomen ellipsoid (m2), and its variants in grams.
m2_LSF_p5 = rho_butt*(4/3)*pi*(bbutt_LSF_p5^2)*abutt_LSF_p5;
m2_LSF_1 = rho_butt*(4/3)*pi*(bbutt_LSF_1^2)*abutt_LSF_1;
m2_LSF_2 = rho_butt*(4/3)*pi*(bbutt_LSF_2^2)*abutt_LSF_2;

%% Assign the imported numbers to internal variables
disp('Assign internal variables')

hws = 500; 
PartPath = 2000;
timestep = numel(Tstore(1:(end-100)))/hws;
tderiv = Tstore(2);
t_end = Tstore(end-100);

%Create the sum of primes signal
signal_amp = 5; %in cm 
prime_f = [0.2, 0.3, 0.5, 0.7, 1.1, 1.7, 2.9, 4.3, 7.9, 13.7, 19.9]; %in Hz
prime_a = (signal_amp./(2*pi.*prime_f)).*(2*pi.*prime_f(1)); %in cm
prime_ph = zeros(1,numel(prime_f)); %Not sure if we'll require a phase

%Goal criteria (for the cost function)
y_g = prime_a(1)*sin(2*pi*prime_f(1)*Tstore(1:(end-100)) + prime_ph(1)) +...
    prime_a(2)*sin(2*pi*prime_f(2)*Tstore(1:(end-100)) + prime_ph(2)) +...
    prime_a(3)*sin(2*pi*prime_f(3)*Tstore(1:(end-100)) + prime_ph(3)) +...
    prime_a(4)*sin(2*pi*prime_f(4)*Tstore(1:(end-100)) + prime_ph(4)) +...
    prime_a(5)*sin(2*pi*prime_f(5)*Tstore(1:(end-100)) + prime_ph(5)) +...
    prime_a(6)*sin(2*pi*prime_f(6)*Tstore(1:(end-100)) + prime_ph(6)) +...
    prime_a(7)*sin(2*pi*prime_f(7)*Tstore(1:(end-100)) + prime_ph(7)) +...
    prime_a(8)*sin(2*pi*prime_f(8)*Tstore(1:(end-100)) + prime_ph(8)) +...
    prime_a(9)*sin(2*pi*prime_f(9)*Tstore(1:(end-100)) + prime_ph(9)) +...
    prime_a(10)*sin(2*pi*prime_f(10)*Tstore(1:(end-100)) + prime_ph(10)) +...
    prime_a(11)*sin(2*pi*prime_f(11)*Tstore(1:(end-100)) + prime_ph(11));

ydot_g = 2*pi*prime_a(1)*prime_f(1)*cos(2*pi*prime_f(1)*Tstore(1:(end-100)) + prime_ph(1)) +...
    2*pi*prime_a(2)*prime_f(2)*cos(2*pi*prime_f(2)*Tstore(1:(end-100)) + prime_ph(2)) +...
    2*pi*prime_a(3)*prime_f(3)*cos(2*pi*prime_f(3)*Tstore(1:(end-100)) + prime_ph(3)) +...
    2*pi*prime_a(4)*prime_f(4)*cos(2*pi*prime_f(4)*Tstore(1:(end-100)) + prime_ph(4)) +...
    2*pi*prime_a(5)*prime_f(5)*cos(2*pi*prime_f(5)*Tstore(1:(end-100)) + prime_ph(5)) +...
    2*pi*prime_a(6)*prime_f(6)*cos(2*pi*prime_f(6)*Tstore(1:(end-100)) + prime_ph(6)) +...
    2*pi*prime_a(7)*prime_f(7)*cos(2*pi*prime_f(7)*Tstore(1:(end-100)) + prime_ph(7)) +...
    2*pi*prime_a(8)*prime_f(8)*cos(2*pi*prime_f(8)*Tstore(1:(end-100)) + prime_ph(8)) +...
    2*pi*prime_a(9)*prime_f(9)*cos(2*pi*prime_f(9)*Tstore(1:(end-100)) + prime_ph(9)) +...
    2*pi*prime_a(10)*prime_f(10)*cos(2*pi*prime_f(10)*Tstore(1:(end-100)) + prime_ph(10)) +...
    2*pi*prime_a(11)*prime_f(11)*cos(2*pi*prime_f(11)*Tstore(1:(end-100)) + prime_ph(11));

y_g = y_g';
x_g = zeros(1,numel(y_g));
theta_g = zeros(1,numel(y_g));
theta_g(1,:) = pi/4;

ydot_g = ydot_g';
xdot_g = zeros(1,numel(y_g));
thetadot_g = zeros(1,numel(y_g));

time_hws(1:PartPath) = Tstore(1:(timestep*hws/PartPath):(end-100));

%% Conversion to mks

%Units of force WERE in g*cm/(s^2), now converted to 10^(-5) Newtons
conv_F_mks = 1/(1000*100); 

%Units of torque WERE in XX, now converted to N*m
conv_torque_mks = 1/((100^2)*1000);

%% Import all the variables
disp('Importing the relevant data')

%Magnitude of applied forces
load('../CuratedData_MPC/fa/LSF_p5/F_fa_LSF_p5_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_p5/F_fa_LSF_p5_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_p5/F_fa_LSF_p5_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/F_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/F_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/F_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_2/F_fa_LSF_2_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_2/F_fa_LSF_2_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_2/F_fa_LSF_2_LE_p4.mat')

%Direction of applied forces
load('../CuratedData_MPC/fa/LSF_p5/alpha_fa_LSF_p5_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_p5/alpha_fa_LSF_p5_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_p5/alpha_fa_LSF_p5_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/alpha_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/alpha_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/alpha_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_2/alpha_fa_LSF_2_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_2/alpha_fa_LSF_2_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_2/alpha_fa_LSF_2_LE_p4.mat')

%Abdominal torque
load('../CuratedData_MPC/fa/LSF_p5/tauAbdo_fa_LSF_p5_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_p5/tauAbdo_fa_LSF_p5_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_p5/tauAbdo_fa_LSF_p5_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/tauAbdo_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/tauAbdo_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/tauAbdo_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_2/tauAbdo_fa_LSF_2_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_2/tauAbdo_fa_LSF_2_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_2/tauAbdo_fa_LSF_2_LE_p4.mat')

%Wing torque
load('../CuratedData_MPC/fa/LSF_p5/tauWing_fa_LSF_p5_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_p5/tauWing_fa_LSF_p5_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_p5/tauWing_fa_LSF_p5_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/tauWing_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/tauWing_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/tauWing_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_2/tauWing_fa_LSF_2_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_2/tauWing_fa_LSF_2_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_2/tauWing_fa_LSF_2_LE_p4.mat')

%x position of pin joint
load('../CuratedData_MPC/fa/LSF_p5/Winstore_x_fa_LSF_p5_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_p5/Winstore_x_fa_LSF_p5_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_p5/Winstore_x_fa_LSF_p5_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_x_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_x_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_x_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_2/Winstore_x_fa_LSF_2_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_2/Winstore_x_fa_LSF_2_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_2/Winstore_x_fa_LSF_2_LE_p4.mat')

%y position of pin joint
load('../CuratedData_MPC/fa/LSF_p5/Winstore_y_fa_LSF_p5_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_p5/Winstore_y_fa_LSF_p5_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_p5/Winstore_y_fa_LSF_p5_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_y_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_y_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_y_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_2/Winstore_y_fa_LSF_2_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_2/Winstore_y_fa_LSF_2_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_2/Winstore_y_fa_LSF_2_LE_p4.mat')

%theta (motion of head-thorax mass with respect to the horizontal)
load('../CuratedData_MPC/fa/LSF_p5/Winstore_theta_fa_LSF_p5_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_p5/Winstore_theta_fa_LSF_p5_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_p5/Winstore_theta_fa_LSF_p5_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_theta_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_theta_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_theta_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_2/Winstore_theta_fa_LSF_2_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_2/Winstore_theta_fa_LSF_2_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_2/Winstore_theta_fa_LSF_2_LE_p4.mat')

%beta (flexion--difference between abdominal & head-thorax mass angles)
load('../CuratedData_MPC/fa/LSF_p5/Winstore_beta_fa_LSF_p5_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_p5/Winstore_beta_fa_LSF_p5_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_p5/Winstore_beta_fa_LSF_p5_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_beta_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_beta_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_beta_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_2/Winstore_beta_fa_LSF_2_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_2/Winstore_beta_fa_LSF_2_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_2/Winstore_beta_fa_LSF_2_LE_p4.mat')

%tracking error
load('../CuratedData_MPC/fa/LSF_p5/trackingerror_fa_LSF_p5_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_p5/trackingerror_fa_LSF_p5_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_p5/trackingerror_fa_LSF_p5_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/trackingerror_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/trackingerror_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/trackingerror_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_2/trackingerror_fa_LSF_2_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_2/trackingerror_fa_LSF_2_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_2/trackingerror_fa_LSF_2_LE_p4.mat')

%% Positional vectors
%IT ALSO LOOKS LIKE WE WON'T NEED L1!
% %Pack the parameter struct (NOT to be altered)
% L1 = LSF_1*0.9; %Length from the thorax-petiole joint to 
%     %the center of the head-thorax in cm

%The following lines were muted out because L3 is ONLY used in 
% %If we *want* the L3 off of the m1 center of mass:
% L3_shifted = LSF_1*0.75; %Length from the thorax-petiole joint 
%     %to the aerodynamic force vector in cm 

%% Preamble to calculating works
%Note (2021/05/23), this preamble also includes other treatments used in
%paper 1. THIS PAPER WILL ONLY USE THE FULLY ACTUATED (fa) TREATMENT

%The work is divided into two main components: work by applied torques 
%(i.e. rotational work), and work by applied force (i.e. rectilinear work).

%Rotational work
%The rotational work done by applied torques are further decomposed into 
%three components:
    %1. Torque by the applied force (IF the applied force IS NOT on the 
    %center of mass for the head-thorax). This is multiplied by the
    %rotation of the head-thorax mass (theta). 
    %**This will only be applicable to treatments "fs" and "us"**
    
    %2. Wing torque, which is by definition, applied about the center of
    %mass of the head-thorax. This is also multiplied by the rotation of 
    %the head-thorax mass (theta). 
    %**This will only be applicable to treatments "fa" and "fs"**
    
    %3. Abdominal torque, which is by definition applied about the pin
    %joint joining the head-thorax mass and the abdomen mass. This is
    %multiplied by the flexion angle (beta).

%Rectilinear work
%The rectilinear work is done only by the applied force, and is a product
%of the applied force with the Euclidean distance the center of mass of the
%head-thorax travels in the given time period.

%% Prescribe state variable matrices

interval = floor((numel(y_g)+1)/PartPath);

%LSF_p5, LE_0
x_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,PartPath);
y_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,PartPath);
theta_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,PartPath);
beta_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,PartPath);
r_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,PartPath);
work_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,PartPath);
work_nd_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,PartPath);

wingWork_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,PartPath);
abdoWork_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,PartPath);
appFWork_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,PartPath);

wingWork_nd_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,PartPath);
abdoWork_nd_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,PartPath);
appFWork_nd_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,PartPath);

%LSF_p5, LE_p2
x_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,PartPath);
y_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,PartPath);
theta_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,PartPath);
beta_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,PartPath);
r_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,PartPath); 
work_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,PartPath);
work_nd_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,PartPath);

wingWork_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,PartPath);
abdoWork_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,PartPath);
appFWork_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,PartPath);

wingWork_nd_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,PartPath);
abdoWork_nd_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,PartPath);
appFWork_nd_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,PartPath);

%LSF_p5, LE_p4
x_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,PartPath);
y_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,PartPath);
theta_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,PartPath);
beta_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,PartPath);
r_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,PartPath);
work_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,PartPath);
work_nd_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,PartPath);

wingWork_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,PartPath);
abdoWork_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,PartPath);
appFWork_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,PartPath);

wingWork_nd_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,PartPath);
abdoWork_nd_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,PartPath);
appFWork_nd_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,PartPath);

%LSF_1, LE_0
x_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,PartPath);
y_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,PartPath);
theta_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,PartPath);
beta_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,PartPath);
r_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,PartPath);
work_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,PartPath);
work_nd_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,PartPath);

wingWork_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,PartPath);
abdoWork_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,PartPath);
appFWork_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,PartPath);

wingWork_nd_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,PartPath);
abdoWork_nd_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,PartPath);
appFWork_nd_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,PartPath);

%LSF_1, LE_p2
x_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,PartPath);
y_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,PartPath);
theta_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,PartPath);
beta_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,PartPath);
r_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,PartPath);
work_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,PartPath);
work_nd_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,PartPath);

wingWork_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,PartPath);
abdoWork_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,PartPath);
appFWork_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,PartPath);

wingWork_nd_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,PartPath);
abdoWork_nd_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,PartPath);
appFWork_nd_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,PartPath);

%LSF_1, LE_p4
x_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,PartPath);
y_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,PartPath);
theta_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,PartPath);
beta_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,PartPath);
r_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,PartPath);
work_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,PartPath);
work_nd_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,PartPath);

wingWork_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,PartPath);
abdoWork_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,PartPath);
appFWork_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,PartPath);

wingWork_nd_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,PartPath);
abdoWork_nd_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,PartPath);
appFWork_nd_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,PartPath);

%LSF_2, LE_0
x_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,PartPath);
y_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,PartPath);
theta_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,PartPath);
beta_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,PartPath);
r_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,PartPath);
work_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,PartPath);
work_nd_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,PartPath);

wingWork_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,PartPath);
abdoWork_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,PartPath);
appFWork_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,PartPath);

wingWork_nd_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,PartPath);
abdoWork_nd_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,PartPath);
appFWork_nd_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,PartPath);

%LSF_2, LE_p2
x_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,PartPath);
y_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,PartPath);
theta_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,PartPath);
beta_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,PartPath);
r_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,PartPath);
work_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,PartPath);
work_nd_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,PartPath);

wingWork_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,PartPath);
abdoWork_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,PartPath);
appFWork_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,PartPath);

wingWork_nd_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,PartPath);
abdoWork_nd_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,PartPath);
appFWork_nd_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,PartPath);

%LSF_2, LE_p4
x_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,PartPath);
y_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,PartPath);
theta_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,PartPath);
beta_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,PartPath);
r_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,PartPath);
work_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,PartPath);
work_nd_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,PartPath);

wingWork_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,PartPath);
abdoWork_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,PartPath);
appFWork_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,PartPath);

wingWork_nd_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,PartPath);
abdoWork_nd_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,PartPath);
appFWork_nd_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,PartPath);

%% Downsize the matrices as appropriate and calculate displacement

%This downsizing only incorporates the x, y, head motion (theta), and 
%flexion (beta) values at the end of the 20 ms time interval.

%LSF_p5, LE_0
x_LSF_p5_LE_0_pre = Winstore_x_fa_LSF_p5_LE_0(:,1:interval:end);
y_LSF_p5_LE_0_pre = Winstore_y_fa_LSF_p5_LE_0(:,1:interval:end);
theta_LSF_p5_LE_0_pre = Winstore_theta_fa_LSF_p5_LE_0(:,1:interval:end);
beta_LSF_p5_LE_0_pre = Winstore_beta_fa_LSF_p5_LE_0(:,1:interval:end);

%LSF_p5, LE_p2
x_LSF_p5_LE_p2_pre = Winstore_x_fa_LSF_p5_LE_p2(:,1:interval:end);
y_LSF_p5_LE_p2_pre = Winstore_y_fa_LSF_p5_LE_p2(:,1:interval:end);
theta_LSF_p5_LE_p2_pre = Winstore_theta_fa_LSF_p5_LE_p2(:,1:interval:end);
beta_LSF_p5_LE_p2_pre = Winstore_beta_fa_LSF_p5_LE_p2(:,1:interval:end);

%LSF_p5, LE_p4
x_LSF_p5_LE_p4_pre = Winstore_x_fa_LSF_p5_LE_p4(:,1:interval:end);
y_LSF_p5_LE_p4_pre = Winstore_y_fa_LSF_p5_LE_p4(:,1:interval:end);
theta_LSF_p5_LE_p4_pre = Winstore_theta_fa_LSF_p5_LE_p4(:,1:interval:end);
beta_LSF_p5_LE_p4_pre = Winstore_beta_fa_LSF_p5_LE_p4(:,1:interval:end);

%LSF_1, LE_0
x_LSF_1_LE_0_pre = Winstore_x_fa_LSF_1_LE_0(:,1:interval:end);
y_LSF_1_LE_0_pre = Winstore_y_fa_LSF_1_LE_0(:,1:interval:end);
theta_LSF_1_LE_0_pre = Winstore_theta_fa_LSF_1_LE_0(:,1:interval:end);
beta_LSF_1_LE_0_pre = Winstore_beta_fa_LSF_1_LE_0(:,1:interval:end);

%LSF_1_LE_p2
x_LSF_1_LE_p2_pre = Winstore_x_fa_LSF_1_LE_p2(:,1:interval:end);
y_LSF_1_LE_p2_pre = Winstore_y_fa_LSF_1_LE_p2(:,1:interval:end);
theta_LSF_1_LE_p2_pre = Winstore_theta_fa_LSF_1_LE_p2(:,1:interval:end);
beta_LSF_1_LE_p2_pre = Winstore_beta_fa_LSF_1_LE_p2(:,1:interval:end);

%LSF_1, LE_p4
x_LSF_1_LE_p4_pre = Winstore_x_fa_LSF_1_LE_p4(:,1:interval:end);
y_LSF_1_LE_p4_pre = Winstore_y_fa_LSF_1_LE_p4(:,1:interval:end);
theta_LSF_1_LE_p4_pre = Winstore_theta_fa_LSF_1_LE_p4(:,1:interval:end);
beta_LSF_1_LE_p4_pre = Winstore_beta_fa_LSF_1_LE_p4(:,1:interval:end);

%LSF_2, LE_0
x_LSF_2_LE_0_pre = Winstore_x_fa_LSF_2_LE_0(:,1:interval:end);
y_LSF_2_LE_0_pre = Winstore_y_fa_LSF_2_LE_0(:,1:interval:end);
theta_LSF_2_LE_0_pre = Winstore_theta_fa_LSF_2_LE_0(:,1:interval:end);
beta_LSF_2_LE_0_pre = Winstore_beta_fa_LSF_2_LE_0(:,1:interval:end);

%LSF_2, LE_p2
x_LSF_2_LE_p2_pre = Winstore_x_fa_LSF_2_LE_p2(:,1:interval:end);
y_LSF_2_LE_p2_pre = Winstore_y_fa_LSF_2_LE_p2(:,1:interval:end);
theta_LSF_2_LE_p2_pre = Winstore_theta_fa_LSF_2_LE_p2(:,1:interval:end);
beta_LSF_2_LE_p2_pre = Winstore_beta_fa_LSF_2_LE_p2(:,1:interval:end);

%LSF_2, LE_p4
x_LSF_2_LE_p4_pre = Winstore_x_fa_LSF_2_LE_p4(:,1:interval:end);
y_LSF_2_LE_p4_pre = Winstore_y_fa_LSF_2_LE_p4(:,1:interval:end);
theta_LSF_2_LE_p4_pre = Winstore_theta_fa_LSF_2_LE_p4(:,1:interval:end);
beta_LSF_2_LE_p4_pre = Winstore_beta_fa_LSF_2_LE_p4(:,1:interval:end);

%LSF_p5_LE_0
if NumOf_LSF_p5_LE_0_Files > 0
    for i = 1:NumOf_LSF_p5_LE_0_Files
        x_LSF_p5_LE_0(i,:) = diff(x_LSF_p5_LE_0_pre(i,:)); %in cm
        y_LSF_p5_LE_0(i,:) = diff(y_LSF_p5_LE_0_pre(i,:)); %in cm
        theta_LSF_p5_LE_0(i,:) = diff(theta_LSF_p5_LE_0_pre(i,:)); %in radians
        beta_LSF_p5_LE_0(i,:) = diff(beta_LSF_p5_LE_0_pre(i,:)); %in radians
    end
    
    r_LSF_p5_LE_0 = sqrt(x_LSF_p5_LE_0.^2 + y_LSF_p5_LE_0.^2); %in cm
end

%LSF_p5_LE_p2
if NumOf_LSF_p5_LE_p2_Files > 0
    for i = 1:NumOf_LSF_p5_LE_p2_Files
        x_LSF_p5_LE_p2(i,:) = diff(x_LSF_p5_LE_p2_pre(i,:)); %in cm
        y_LSF_p5_LE_p2(i,:) = diff(y_LSF_p5_LE_p2_pre(i,:)); %in cm
        theta_LSF_p5_LE_p2(i,:) = diff(theta_LSF_p5_LE_p2_pre(i,:)); %in radians
        beta_LSF_p5_LE_p2(i,:) = diff(beta_LSF_p5_LE_p2_pre(i,:)); %in radians
    end
    
    r_LSF_p5_LE_p2 = sqrt(x_LSF_p5_LE_p2.^2 + y_LSF_p5_LE_p2.^2); %in cm
end

%LSF_p5_LE_p4
if NumOf_LSF_p5_LE_p4_Files > 0
    for i = 1:NumOf_LSF_p5_LE_p4_Files
        x_LSF_p5_LE_p4(i,:) = diff(x_LSF_p5_LE_p4_pre(i,:)); %in cm
        y_LSF_p5_LE_p4(i,:) = diff(y_LSF_p5_LE_p4_pre(i,:)); %in cm
        theta_LSF_p5_LE_p4(i,:) = diff(theta_LSF_p5_LE_p4_pre(i,:)); %in radians
        beta_LSF_p5_LE_p4(i,:) = diff(beta_LSF_p5_LE_p4_pre(i,:)); %in radians
    end
    
    r_LSF_p5_LE_p4 = sqrt(x_LSF_p5_LE_p4.^2 + y_LSF_p5_LE_p4.^2); %in cm
end

%LSF_1_LE_0
if NumOf_LSF_1_LE_0_Files > 0
    for i = 1:NumOf_LSF_1_LE_0_Files
        x_LSF_1_LE_0(i,:) = diff(x_LSF_1_LE_0_pre(i,:)); %in cm
        y_LSF_1_LE_0(i,:) = diff(y_LSF_1_LE_0_pre(i,:)); %in cm
        theta_LSF_1_LE_0(i,:) = diff(theta_LSF_1_LE_0_pre(i,:)); %in radians
        beta_LSF_1_LE_0(i,:) = diff(beta_LSF_1_LE_0_pre(i,:)); %in radians
    end
    
    r_LSF_1_LE_0 = sqrt(x_LSF_1_LE_0.^2 + y_LSF_1_LE_0.^2); %in cm
end

%LSF_1_LE_p2
if NumOf_LSF_1_LE_p2_Files > 0
    for i = 1:NumOf_LSF_1_LE_p2_Files
        x_LSF_1_LE_p2(i,:) = diff(x_LSF_1_LE_p2_pre(i,:)); %in cm
        y_LSF_1_LE_p2(i,:) = diff(y_LSF_1_LE_p2_pre(i,:)); %in cm
        theta_LSF_1_LE_p2(i,:) = diff(theta_LSF_1_LE_p2_pre(i,:)); %in radians
        beta_LSF_1_LE_p2(i,:) = diff(beta_LSF_1_LE_p2_pre(i,:)); %in radians
    end
    
    r_LSF_1_LE_p2 = sqrt(x_LSF_1_LE_p2.^2 + y_LSF_1_LE_p2.^2); %in cm
end

%LSF_1_LE_p4
if NumOf_LSF_1_LE_p4_Files > 0
    for i = 1:NumOf_LSF_1_LE_p4_Files
        x_LSF_1_LE_p4(i,:) = diff(x_LSF_1_LE_p4_pre(i,:)); %in cm
        y_LSF_1_LE_p4(i,:) = diff(y_LSF_1_LE_p4_pre(i,:)); %in cm
        theta_LSF_1_LE_p4(i,:) = diff(theta_LSF_1_LE_p4_pre(i,:)); %in radians
        beta_LSF_1_LE_p4(i,:) = diff(beta_LSF_1_LE_p4_pre(i,:)); %in radians
    end
    
    r_LSF_1_LE_p4 = sqrt(x_LSF_1_LE_p4.^2 + y_LSF_1_LE_p4.^2); %in cm
end

%LSF_2_LE_0
if NumOf_LSF_2_LE_0_Files > 0
    for i = 1:NumOf_LSF_2_LE_0_Files
        x_LSF_2_LE_0(i,:) = diff(x_LSF_2_LE_0_pre(i,:)); %in cm
        y_LSF_2_LE_0(i,:) = diff(y_LSF_2_LE_0_pre(i,:)); %in cm
        theta_LSF_2_LE_0(i,:) = diff(theta_LSF_2_LE_0_pre(i,:)); %in radians
        beta_LSF_2_LE_0(i,:) = diff(beta_LSF_2_LE_0_pre(i,:)); %in radians
    end
    
    r_LSF_2_LE_0 = sqrt(x_LSF_2_LE_0.^2 + y_LSF_2_LE_0.^2); %in cm
end

%LSF_2_LE_p2
if NumOf_LSF_2_LE_p2_Files > 0
    for i = 1:NumOf_LSF_2_LE_p2_Files
        x_LSF_2_LE_p2(i,:) = diff(x_LSF_2_LE_p2_pre(i,:)); %in cm
        y_LSF_2_LE_p2(i,:) = diff(y_LSF_2_LE_p2_pre(i,:)); %in cm
        theta_LSF_2_LE_p2(i,:) = diff(theta_LSF_2_LE_p2_pre(i,:)); %in radians
        beta_LSF_2_LE_p2(i,:) = diff(beta_LSF_2_LE_p2_pre(i,:)); %in radians
    end
    
    r_LSF_2_LE_p2 = sqrt(x_LSF_2_LE_p2.^2 + y_LSF_2_LE_p2.^2); %in cm
end

%LSF_2_LE_p4
if NumOf_LSF_2_LE_p4_Files > 0
    for i = 1:NumOf_LSF_2_LE_p4_Files
        x_LSF_2_LE_p4(i,:) = diff(x_LSF_2_LE_p4_pre(i,:)); %in cm
        y_LSF_2_LE_p4(i,:) = diff(y_LSF_2_LE_p4_pre(i,:)); %in cm
        theta_LSF_2_LE_p4(i,:) = diff(theta_LSF_2_LE_p4_pre(i,:)); %in radians
        beta_LSF_2_LE_p4(i,:) = diff(beta_LSF_2_LE_p4_pre(i,:)); %in radians
    end
    
    r_LSF_2_LE_p4 = sqrt(x_LSF_2_LE_p4.^2 + y_LSF_2_LE_p4.^2); %in cm
end

%% Calculate the mechanical work (raw data)

%From original script, this was the fully actuated ("fa") section
    %WHICH IS THE ONLY TREATMENT WE WILL NEED FOR THIS SCRIPT (and PAPER 2)
%LSF_p5_LE_0
if NumOf_LSF_p5_LE_0_Files > 0
    work_LSF_p5_LE_0 = 0 ...
    + abs((conv_torque_mks.*tauWing_fa_LSF_p5_LE_0).*theta_LSF_p5_LE_0)...
        + abs((conv_torque_mks.*tauAbdo_fa_LSF_p5_LE_0).*beta_LSF_p5_LE_0)...
        + abs((conv_F_mks.*F_fa_LSF_p5_LE_0).*(r_LSF_p5_LE_0/100)); %in Joules
    
    wingWork_LSF_p5_LE_0 =...
        abs((conv_torque_mks.*tauWing_fa_LSF_p5_LE_0).*theta_LSF_p5_LE_0);
    abdoWork_LSF_p5_LE_0 =...
        abs((conv_torque_mks.*tauAbdo_fa_LSF_p5_LE_0).*beta_LSF_p5_LE_0);
    appFWork_LSF_p5_LE_0 =...
        abs((conv_F_mks.*F_fa_LSF_p5_LE_0).*(r_LSF_p5_LE_0/100));
end

%LSF_p5_LE_p2
if NumOf_LSF_p5_LE_p2_Files > 0
    work_LSF_p5_LE_p2 = 0 ...
    + abs((conv_torque_mks.*tauWing_fa_LSF_p5_LE_p2).*theta_LSF_p5_LE_p2)...
        + abs((conv_torque_mks.*tauAbdo_fa_LSF_p5_LE_p2).*beta_LSF_p5_LE_p2)...
        + abs((conv_F_mks.*F_fa_LSF_p5_LE_p2).*(r_LSF_p5_LE_p2/100)); %in Joules
    
    wingWork_LSF_p5_LE_p2 =...
        abs((conv_torque_mks.*tauWing_fa_LSF_p5_LE_p2).*theta_LSF_p5_LE_p2);
    abdoWork_LSF_p5_LE_p2 =...
        abs((conv_torque_mks.*tauAbdo_fa_LSF_p5_LE_p2).*beta_LSF_p5_LE_p2);
    appFWork_LSF_p5_LE_p2 =...
        abs((conv_F_mks.*F_fa_LSF_p5_LE_p2).*(r_LSF_p5_LE_p2/100));
end

%LSF_p5_LE_p4
if NumOf_LSF_p5_LE_p4_Files > 0
    work_LSF_p5_LE_p4 = 0 ...
    + abs((conv_torque_mks.*tauWing_fa_LSF_p5_LE_p4).*theta_LSF_p5_LE_p4)...
        + abs((conv_torque_mks.*tauAbdo_fa_LSF_p5_LE_p4).*beta_LSF_p5_LE_p4)...
        + abs((conv_F_mks.*F_fa_LSF_p5_LE_p4).*(r_LSF_p5_LE_p4/100)); %in Joules
    
    wingWork_LSF_p5_LE_p4 =...
        abs((conv_torque_mks.*tauWing_fa_LSF_p5_LE_p4).*theta_LSF_p5_LE_p4);
    abdoWork_LSF_p5_LE_p4 =...
        abs((conv_torque_mks.*tauAbdo_fa_LSF_p5_LE_p4).*beta_LSF_p5_LE_p4);
    appFWork_LSF_p5_LE_p4 =...
        abs((conv_F_mks.*F_fa_LSF_p5_LE_p4).*(r_LSF_p5_LE_p4/100));
end

%LSF_1_LE_0
if NumOf_LSF_1_LE_0_Files > 0
    work_LSF_1_LE_0 = 0 ...
    + abs((conv_torque_mks.*tauWing_fa_LSF_1_LE_0).*theta_LSF_1_LE_0)...
        + abs((conv_torque_mks.*tauAbdo_fa_LSF_1_LE_0).*beta_LSF_1_LE_0)...
        + abs((conv_F_mks.*F_fa_LSF_1_LE_0).*(r_LSF_1_LE_0/100)); %in Joules
    
    wingWork_LSF_1_LE_0 =...
        abs((conv_torque_mks.*tauWing_fa_LSF_1_LE_0).*theta_LSF_1_LE_0);
    abdoWork_LSF_1_LE_0 =...
        abs((conv_torque_mks.*tauAbdo_fa_LSF_1_LE_0).*beta_LSF_1_LE_0);
    appFWork_LSF_1_LE_0 =...
        abs((conv_F_mks.*F_fa_LSF_1_LE_0).*(r_LSF_1_LE_0/100));
end

%LSF_1_LE_p2
if NumOf_LSF_1_LE_p2_Files > 0
    work_LSF_1_LE_p2 = 0 ...
    + abs((conv_torque_mks.*tauWing_fa_LSF_1_LE_p2).*theta_LSF_1_LE_p2)...
        + abs((conv_torque_mks.*tauAbdo_fa_LSF_1_LE_p2).*beta_LSF_1_LE_p2)...
        + abs((conv_F_mks.*F_fa_LSF_1_LE_p2).*(r_LSF_1_LE_p2/100)); %in Joules
    
    wingWork_LSF_1_LE_p2 =...
        abs((conv_torque_mks.*tauWing_fa_LSF_1_LE_p2).*theta_LSF_1_LE_p2);
    abdoWork_LSF_1_LE_p2 =...
        abs((conv_torque_mks.*tauAbdo_fa_LSF_1_LE_p2).*beta_LSF_1_LE_p2);
    appFWork_LSF_1_LE_p2 =...
        abs((conv_F_mks.*F_fa_LSF_1_LE_p2).*(r_LSF_1_LE_p2/100));
end

%LSF_1_LE_p4
if NumOf_LSF_1_LE_p4_Files > 0
    work_LSF_1_LE_p4 = 0 ...
    + abs((conv_torque_mks.*tauWing_fa_LSF_1_LE_p4).*theta_LSF_1_LE_p4)...
        + abs((conv_torque_mks.*tauAbdo_fa_LSF_1_LE_p4).*beta_LSF_1_LE_p4)...
        + abs((conv_F_mks.*F_fa_LSF_1_LE_p4).*(r_LSF_1_LE_p4/100)); %in Joules
    
    wingWork_LSF_1_LE_p4 =...
        abs((conv_torque_mks.*tauWing_fa_LSF_1_LE_p4).*theta_LSF_1_LE_p4);
    abdoWork_LSF_1_LE_p4 =...
        abs((conv_torque_mks.*tauAbdo_fa_LSF_1_LE_p4).*beta_LSF_1_LE_p4);
    appFWork_LSF_1_LE_p4 =...
        abs((conv_F_mks.*F_fa_LSF_1_LE_p4).*(r_LSF_1_LE_p4/100));
end

%LSF_2_LE_0
if NumOf_LSF_2_LE_0_Files > 0
    work_LSF_2_LE_0 = 0 ...
    + abs((conv_torque_mks.*tauWing_fa_LSF_2_LE_0).*theta_LSF_2_LE_0)...
        + abs((conv_torque_mks.*tauAbdo_fa_LSF_2_LE_0).*beta_LSF_2_LE_0)...
        + abs((conv_F_mks.*F_fa_LSF_2_LE_0).*(r_LSF_2_LE_0/100)); %in Joules
    
    wingWork_LSF_2_LE_0 =...
        abs((conv_torque_mks.*tauWing_fa_LSF_2_LE_0).*theta_LSF_2_LE_0);
    abdoWork_LSF_2_LE_0 =...
        abs((conv_torque_mks.*tauAbdo_fa_LSF_2_LE_0).*beta_LSF_2_LE_0);
    appFWork_LSF_2_LE_0 =...
        abs((conv_F_mks.*F_fa_LSF_2_LE_0).*(r_LSF_2_LE_0/100));
end

%LSF_2_LE_p2
if NumOf_LSF_2_LE_p2_Files > 0
    work_LSF_2_LE_p2 = 0 ...
    + abs((conv_torque_mks.*tauWing_fa_LSF_2_LE_p2).*theta_LSF_2_LE_p2)...
        + abs((conv_torque_mks.*tauAbdo_fa_LSF_2_LE_p2).*beta_LSF_2_LE_p2)...
        + abs((conv_F_mks.*F_fa_LSF_2_LE_p2).*(r_LSF_2_LE_p2/100)); %in Joules
    
    wingWork_LSF_2_LE_p2 =...
        abs((conv_torque_mks.*tauWing_fa_LSF_2_LE_p2).*theta_LSF_2_LE_p2);
    abdoWork_LSF_2_LE_p2 =...
        abs((conv_torque_mks.*tauAbdo_fa_LSF_2_LE_p2).*beta_LSF_2_LE_p2);
    appFWork_LSF_2_LE_p2 =...
        abs((conv_F_mks.*F_fa_LSF_2_LE_p2).*(r_LSF_2_LE_p2/100));
end

%LSF_2_LE_p4
if NumOf_LSF_2_LE_p4_Files > 0
    work_LSF_2_LE_p4 = 0 ...
    + abs((conv_torque_mks.*tauWing_fa_LSF_2_LE_p4).*theta_LSF_2_LE_p4)...
        + abs((conv_torque_mks.*tauAbdo_fa_LSF_2_LE_p4).*beta_LSF_2_LE_p4)...
        + abs((conv_F_mks.*F_fa_LSF_2_LE_p4).*(r_LSF_2_LE_p4/100)); %in Joules
    
    wingWork_LSF_2_LE_p4 =...
        abs((conv_torque_mks.*tauWing_fa_LSF_2_LE_p4).*theta_LSF_2_LE_p4);
    abdoWork_LSF_2_LE_p4 =...
        abs((conv_torque_mks.*tauAbdo_fa_LSF_2_LE_p4).*beta_LSF_2_LE_p4);
    appFWork_LSF_2_LE_p4 =...
        abs((conv_F_mks.*F_fa_LSF_2_LE_p4).*(r_LSF_2_LE_p4/100));
end

%% Calculate the non-dimensionalized mechanical work

%From original script, this was the fully actuated ("fa") section
    %WHICH IS THE ONLY TREATMENT WE WILL NEED FOR THIS SCRIPT (and PAPER 2)
%LSF_p5_LE_0
if NumOf_LSF_p5_LE_0_Files > 0
    work_nd_LSF_p5_LE_0 = (0 ...
    + abs((conv_torque_mks.*tauWing_fa_LSF_p5_LE_0).*theta_LSF_p5_LE_0)...
        + abs((conv_torque_mks.*tauAbdo_fa_LSF_p5_LE_0).*beta_LSF_p5_LE_0)...
        + abs((conv_F_mks.*F_fa_LSF_p5_LE_0).*(r_LSF_p5_LE_0/100)))./...
        (((m1_LSF_p5/2)/1000)*(g/100)*(r_LSF_p5_LE_0/100)); %unitless
    
    wingWork_nd_LSF_p5_LE_0 =...
        abs((conv_torque_mks.*tauWing_fa_LSF_p5_LE_0).*theta_LSF_p5_LE_0)./...
        (((m1_LSF_p5/2)/1000)*(g/100)*(r_LSF_p5_LE_0/100));
    abdoWork_nd_LSF_p5_LE_0 =...
        abs((conv_torque_mks.*tauAbdo_fa_LSF_p5_LE_0).*beta_LSF_p5_LE_0)./...
        (((m1_LSF_p5/2)/1000)*(g/100)*(r_LSF_p5_LE_0/100));
    appFWork_nd_LSF_p5_LE_0 =...
        abs((conv_F_mks.*F_fa_LSF_p5_LE_0).*(r_LSF_p5_LE_0/100))./...
        (((m1_LSF_p5/2)/1000)*(g/100)*(r_LSF_p5_LE_0/100));
end

%LSF_p5_LE_p2
if NumOf_LSF_p5_LE_p2_Files > 0
    work_nd_LSF_p5_LE_p2 = (0 ...
    + abs((conv_torque_mks.*tauWing_fa_LSF_p5_LE_p2).*theta_LSF_p5_LE_p2)...
        + abs((conv_torque_mks.*tauAbdo_fa_LSF_p5_LE_p2).*beta_LSF_p5_LE_p2)...
        + abs((conv_F_mks.*F_fa_LSF_p5_LE_p2).*(r_LSF_p5_LE_p2/100)))./...
        (((m1_LSF_p5/2)/1000)*(g/100)*(r_LSF_p5_LE_p2/100)); %unitless
    
    wingWork_nd_LSF_p5_LE_p2 =...
        abs((conv_torque_mks.*tauWing_fa_LSF_p5_LE_p2).*theta_LSF_p5_LE_p2)./...
        (((m1_LSF_p5/2)/1000)*(g/100)*(r_LSF_p5_LE_p2/100));
    abdoWork_nd_LSF_p5_LE_p2 =...
        abs((conv_torque_mks.*tauAbdo_fa_LSF_p5_LE_p2).*beta_LSF_p5_LE_p2)./...
        (((m1_LSF_p5/2)/1000)*(g/100)*(r_LSF_p5_LE_p2/100));
    appFWork_nd_LSF_p5_LE_p2 =...
        abs((conv_F_mks.*F_fa_LSF_p5_LE_p2).*(r_LSF_p5_LE_p2/100))./...
        (((m1_LSF_p5/2)/1000)*(g/100)*(r_LSF_p5_LE_p2/100));
end

%LSF_p5_LE_p4
if NumOf_LSF_p5_LE_p4_Files > 0
    work_nd_LSF_p5_LE_p4 = (0 ...
    + abs((conv_torque_mks.*tauWing_fa_LSF_p5_LE_p4).*theta_LSF_p5_LE_p4)...
        + abs((conv_torque_mks.*tauAbdo_fa_LSF_p5_LE_p4).*beta_LSF_p5_LE_p4)...
        + abs((conv_F_mks.*F_fa_LSF_p5_LE_p4).*(r_LSF_p5_LE_p4/100)))./...
        (((m1_LSF_p5/2)/1000)*(g/100)*(r_LSF_p5_LE_p4/100)); %unitless
    
    wingWork_nd_LSF_p5_LE_p4 =...
        abs((conv_torque_mks.*tauWing_fa_LSF_p5_LE_p4).*theta_LSF_p5_LE_p4)./...
        (((m1_LSF_p5/2)/1000)*(g/100)*(r_LSF_p5_LE_p4/100));
    abdoWork_nd_LSF_p5_LE_p4 =...
        abs((conv_torque_mks.*tauAbdo_fa_LSF_p5_LE_p4).*beta_LSF_p5_LE_p4)./...
        (((m1_LSF_p5/2)/1000)*(g/100)*(r_LSF_p5_LE_p4/100));
    appFWork_nd_LSF_p5_LE_p4 =...
        abs((conv_F_mks.*F_fa_LSF_p5_LE_p4).*(r_LSF_p5_LE_p4/100))./...
        (((m1_LSF_p5/2)/1000)*(g/100)*(r_LSF_p5_LE_p4/100));
end

%LSF_1_LE_0
if NumOf_LSF_1_LE_0_Files > 0
    work_nd_LSF_1_LE_0 = (0 ...
    + abs((conv_torque_mks.*tauWing_fa_LSF_1_LE_0).*theta_LSF_1_LE_0)...
        + abs((conv_torque_mks.*tauAbdo_fa_LSF_1_LE_0).*beta_LSF_1_LE_0)...
        + abs((conv_F_mks.*F_fa_LSF_1_LE_0).*(r_LSF_1_LE_0/100)))./...
        (((m1_LSF_1/2)/1000)*(g/100)*(r_LSF_1_LE_0/100)); %unitless
    
    wingWork_nd_LSF_1_LE_0 =...
        abs((conv_torque_mks.*tauWing_fa_LSF_1_LE_0).*theta_LSF_1_LE_0)./...
        (((m1_LSF_1/2)/1000)*(g/100)*(r_LSF_1_LE_0/100));
    abdoWork_nd_LSF_1_LE_0 =...
        abs((conv_torque_mks.*tauAbdo_fa_LSF_1_LE_0).*beta_LSF_1_LE_0)./...
        (((m1_LSF_1/2)/1000)*(g/100)*(r_LSF_1_LE_0/100));
    appFWork_nd_LSF_1_LE_0 =...
        abs((conv_F_mks.*F_fa_LSF_1_LE_0).*(r_LSF_1_LE_0/100))./...
        (((m1_LSF_1/2)/1000)*(g/100)*(r_LSF_1_LE_0/100));
end

%LSF_1_LE_p2
if NumOf_LSF_1_LE_p2_Files > 0
    work_nd_LSF_1_LE_p2 = (0 ...
    + abs((conv_torque_mks.*tauWing_fa_LSF_1_LE_p2).*theta_LSF_1_LE_p2)...
        + abs((conv_torque_mks.*tauAbdo_fa_LSF_1_LE_p2).*beta_LSF_1_LE_p2)...
        + abs((conv_F_mks.*F_fa_LSF_1_LE_p2).*(r_LSF_1_LE_p2/100)))./...
        (((m1_LSF_1/2)/1000)*(g/100)*(r_LSF_1_LE_p2/100)); %unitless
    
    wingWork_nd_LSF_1_LE_p2 =...
        abs((conv_torque_mks.*tauWing_fa_LSF_1_LE_p2).*theta_LSF_1_LE_p2)./...
        (((m1_LSF_1/2)/1000)*(g/100)*(r_LSF_1_LE_p2/100));
    abdoWork_nd_LSF_1_LE_p2 =...
        abs((conv_torque_mks.*tauAbdo_fa_LSF_1_LE_p2).*beta_LSF_1_LE_p2)./...
        (((m1_LSF_1/2)/1000)*(g/100)*(r_LSF_1_LE_p2/100));
    appFWork_nd_LSF_1_LE_p2 =...
        abs((conv_F_mks.*F_fa_LSF_1_LE_p2).*(r_LSF_1_LE_p2/100))./...
        (((m1_LSF_1/2)/1000)*(g/100)*(r_LSF_1_LE_p2/100));
end

%LSF_1_LE_p4
if NumOf_LSF_1_LE_p4_Files > 0
    work_nd_LSF_1_LE_p4 = (0 ...
    + abs((conv_torque_mks.*tauWing_fa_LSF_1_LE_p4).*theta_LSF_1_LE_p4)...
        + abs((conv_torque_mks.*tauAbdo_fa_LSF_1_LE_p4).*beta_LSF_1_LE_p4)...
        + abs((conv_F_mks.*F_fa_LSF_1_LE_p4).*(r_LSF_1_LE_p4/100)))./...
        (((m1_LSF_1/2)/1000)*(g/100)*(r_LSF_1_LE_p4/100)); %unitless
    
    wingWork_nd_LSF_1_LE_p4 =...
        abs((conv_torque_mks.*tauWing_fa_LSF_1_LE_p4).*theta_LSF_1_LE_p4)./...
        (((m1_LSF_1/2)/1000)*(g/100)*(r_LSF_1_LE_p4/100));
    abdoWork_nd_LSF_1_LE_p4 =...
        abs((conv_torque_mks.*tauAbdo_fa_LSF_1_LE_p4).*beta_LSF_1_LE_p4)./...
        (((m1_LSF_1/2)/1000)*(g/100)*(r_LSF_1_LE_p4/100));
    appFWork_nd_LSF_1_LE_p4 =...
        abs((conv_F_mks.*F_fa_LSF_1_LE_p4).*(r_LSF_1_LE_p4/100))./...
        (((m1_LSF_1/2)/1000)*(g/100)*(r_LSF_1_LE_p4/100));
end

%LSF_2_LE_0
if NumOf_LSF_2_LE_0_Files > 0
    work_nd_LSF_2_LE_0 = (0 ...
    + abs((conv_torque_mks.*tauWing_fa_LSF_2_LE_0).*theta_LSF_2_LE_0)...
        + abs((conv_torque_mks.*tauAbdo_fa_LSF_2_LE_0).*beta_LSF_2_LE_0)...
        + abs((conv_F_mks.*F_fa_LSF_2_LE_0).*(r_LSF_2_LE_0/100)))./...
        (((m1_LSF_2/2)/1000)*(g/100)*(r_LSF_2_LE_0/100)); %unitless
    
    wingWork_nd_LSF_2_LE_0 =...
        abs((conv_torque_mks.*tauWing_fa_LSF_2_LE_0).*theta_LSF_2_LE_0)./...
        (((m1_LSF_2/2)/1000)*(g/100)*(r_LSF_2_LE_0/100));
    abdoWork_nd_LSF_2_LE_0 =...
        abs((conv_torque_mks.*tauAbdo_fa_LSF_2_LE_0).*beta_LSF_2_LE_0)./...
        (((m1_LSF_2/2)/1000)*(g/100)*(r_LSF_2_LE_0/100));
    appFWork_nd_LSF_2_LE_0 =...
        abs((conv_F_mks.*F_fa_LSF_2_LE_0).*(r_LSF_2_LE_0/100))./...
        (((m1_LSF_2/2)/1000)*(g/100)*(r_LSF_2_LE_0/100));
end

%LSF_2_LE_p2
if NumOf_LSF_2_LE_p2_Files > 0
    work_nd_LSF_2_LE_p2 = (0 ...
    + abs((conv_torque_mks.*tauWing_fa_LSF_2_LE_p2).*theta_LSF_2_LE_p2)...
        + abs((conv_torque_mks.*tauAbdo_fa_LSF_2_LE_p2).*beta_LSF_2_LE_p2)...
        + abs((conv_F_mks.*F_fa_LSF_2_LE_p2).*(r_LSF_2_LE_p2/100)))./...
        (((m1_LSF_2/2)/1000)*(g/100)*(r_LSF_2_LE_p2/100)); %unitless
    
    wingWork_nd_LSF_2_LE_p2 =...
        abs((conv_torque_mks.*tauWing_fa_LSF_2_LE_p2).*theta_LSF_2_LE_p2)./...
        (((m1_LSF_2/2)/1000)*(g/100)*(r_LSF_2_LE_p2/100));
    abdoWork_nd_LSF_2_LE_p2 =...
        abs((conv_torque_mks.*tauAbdo_fa_LSF_2_LE_p2).*beta_LSF_2_LE_p2)./...
        (((m1_LSF_2/2)/1000)*(g/100)*(r_LSF_2_LE_p2/100));
    appFWork_nd_LSF_2_LE_p2 =...
        abs((conv_F_mks.*F_fa_LSF_2_LE_p2).*(r_LSF_2_LE_p2/100))./...
        (((m1_LSF_2/2)/1000)*(g/100)*(r_LSF_2_LE_p2/100));
end

%LSF_2_LE_p4
if NumOf_LSF_2_LE_p4_Files > 0
    work_nd_LSF_2_LE_p4 = (0 ...
    + abs((conv_torque_mks.*tauWing_fa_LSF_2_LE_p4).*theta_LSF_2_LE_p4)...
        + abs((conv_torque_mks.*tauAbdo_fa_LSF_2_LE_p4).*beta_LSF_2_LE_p4)...
        + abs((conv_F_mks.*F_fa_LSF_2_LE_p4).*(r_LSF_2_LE_p4/100)))./...
        (((m1_LSF_2/2)/1000)*(g/100)*(r_LSF_2_LE_p4/100)); %unitless
    
    wingWork_nd_LSF_2_LE_p4 =...
        abs((conv_torque_mks.*tauWing_fa_LSF_2_LE_p4).*theta_LSF_2_LE_p4)./...
        (((m1_LSF_2/2)/1000)*(g/100)*(r_LSF_2_LE_p4/100));
    abdoWork_nd_LSF_2_LE_p4 =...
        abs((conv_torque_mks.*tauAbdo_fa_LSF_2_LE_p4).*beta_LSF_2_LE_p4)./...
        (((m1_LSF_2/2)/1000)*(g/100)*(r_LSF_2_LE_p4/100));
    appFWork_nd_LSF_2_LE_p4 =...
        abs((conv_F_mks.*F_fa_LSF_2_LE_p4).*(r_LSF_2_LE_p4/100))./...
        (((m1_LSF_2/2)/1000)*(g/100)*(r_LSF_2_LE_p4/100));
end

%% Calculate normalized tracking error

%Import tracking error and normalize it to body length.

%Normalized tracking error -- All values are unitless
%LSF_p5
nd_te_LSF_p5_LE_0 = (trackingerror_fa_LSF_p5_LE_0(:,1:interval:end))./bl_LSF_p5_LE_0; 
nd_te_LSF_p5_LE_p2 = (trackingerror_fa_LSF_p5_LE_p2(:,1:interval:end))./bl_LSF_p5_LE_p2; 
nd_te_LSF_p5_LE_p4 = (trackingerror_fa_LSF_p5_LE_p4(:,1:interval:end))./bl_LSF_p5_LE_p4; 
%LSF_1
nd_te_LSF_1_LE_0 = (trackingerror_fa_LSF_1_LE_0(:,1:interval:end))./bl_LSF_1_LE_0; 
nd_te_LSF_1_LE_p2 = (trackingerror_fa_LSF_1_LE_p2(:,1:interval:end))./bl_LSF_1_LE_p2; 
nd_te_LSF_1_LE_p4 = (trackingerror_fa_LSF_1_LE_p4(:,1:interval:end))./bl_LSF_1_LE_p4; 
%LSF_2
nd_te_LSF_2_LE_0 = (trackingerror_fa_LSF_2_LE_0(:,1:interval:end))./bl_LSF_2_LE_0; 
nd_te_LSF_2_LE_p2 = (trackingerror_fa_LSF_2_LE_p2(:,1:interval:end))./bl_LSF_2_LE_p2; 
nd_te_LSF_2_LE_p4 = (trackingerror_fa_LSF_2_LE_p4(:,1:interval:end))./bl_LSF_2_LE_p4; 

%% Calculate the averages 

%Mean of work terms (and its component terms)
mean_work_LSF_p5_LE_0 = NaN(1,PartPath);
mean_work_LSF_p5_LE_p2 = NaN(1,PartPath);
mean_work_LSF_p5_LE_p4 = NaN(1,PartPath);
mean_work_LSF_1_LE_0 = NaN(1,PartPath);
mean_work_LSF_1_LE_p2 = NaN(1,PartPath);
mean_work_LSF_1_LE_p4 = NaN(1,PartPath);
mean_work_LSF_2_LE_0 = NaN(1,PartPath);
mean_work_LSF_2_LE_p2 = NaN(1,PartPath);
mean_work_LSF_2_LE_p4 = NaN(1,PartPath);

%Standard deviation of work terms (and its component terms)
std_work_LSF_p5_LE_0 = NaN(1,PartPath);
std_work_LSF_p5_LE_p2 = NaN(1,PartPath);
std_work_LSF_p5_LE_p4 = NaN(1,PartPath);
std_work_LSF_1_LE_0 = NaN(1,PartPath);
std_work_LSF_1_LE_p2 = NaN(1,PartPath);
std_work_LSF_1_LE_p4 = NaN(1,PartPath);
std_work_LSF_2_LE_0 = NaN(1,PartPath);
std_work_LSF_2_LE_p2 = NaN(1,PartPath);
std_work_LSF_2_LE_p4 = NaN(1,PartPath);

%Mean of work per full run (wpfr)
mean_wpfr_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,1);
mean_wpfr_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,1);
mean_wpfr_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,1);
mean_wpfr_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,1);
mean_wpfr_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,1);
mean_wpfr_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,1);
mean_wpfr_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,1);
mean_wpfr_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,1);
mean_wpfr_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,1);

mean_nd_wpfr_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,1);
mean_nd_wpfr_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,1);
mean_nd_wpfr_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,1);
mean_nd_wpfr_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,1);
mean_nd_wpfr_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,1);
mean_nd_wpfr_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,1);
mean_nd_wpfr_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,1);
mean_nd_wpfr_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,1);
mean_nd_wpfr_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,1);

%Mean of work per full run (wpfr) per term
mean_wingWork_pfr_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,1);
mean_wingWork_pfr_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,1);
mean_wingWork_pfr_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,1);
mean_wingWork_pfr_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,1);
mean_wingWork_pfr_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,1);
mean_wingWork_pfr_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,1);
mean_wingWork_pfr_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,1);
mean_wingWork_pfr_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,1);
mean_wingWork_pfr_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,1);

mean_abdoWork_pfr_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,1);
mean_abdoWork_pfr_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,1);
mean_abdoWork_pfr_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,1);
mean_abdoWork_pfr_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,1);
mean_abdoWork_pfr_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,1);
mean_abdoWork_pfr_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,1);
mean_abdoWork_pfr_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,1);
mean_abdoWork_pfr_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,1);
mean_abdoWork_pfr_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,1);

mean_appFWork_pfr_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,1);
mean_appFWork_pfr_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,1);
mean_appFWork_pfr_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,1);
mean_appFWork_pfr_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,1);
mean_appFWork_pfr_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,1);
mean_appFWork_pfr_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,1);
mean_appFWork_pfr_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,1);
mean_appFWork_pfr_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,1);
mean_appFWork_pfr_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,1);

mean_nd_wingWork_pfr_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,1);
mean_nd_wingWork_pfr_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,1);
mean_nd_wingWork_pfr_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,1);
mean_nd_wingWork_pfr_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,1);
mean_nd_wingWork_pfr_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,1);
mean_nd_wingWork_pfr_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,1);
mean_nd_wingWork_pfr_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,1);
mean_nd_wingWork_pfr_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,1);
mean_nd_wingWork_pfr_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,1);

mean_nd_abdoWork_pfr_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,1);
mean_nd_abdoWork_pfr_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,1);
mean_nd_abdoWork_pfr_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,1);
mean_nd_abdoWork_pfr_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,1);
mean_nd_abdoWork_pfr_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,1);
mean_nd_abdoWork_pfr_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,1);
mean_nd_abdoWork_pfr_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,1);
mean_nd_abdoWork_pfr_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,1);
mean_nd_abdoWork_pfr_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,1);

mean_nd_appFWork_pfr_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,1);
mean_nd_appFWork_pfr_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,1);
mean_nd_appFWork_pfr_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,1);
mean_nd_appFWork_pfr_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,1);
mean_nd_appFWork_pfr_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,1);
mean_nd_appFWork_pfr_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,1);
mean_nd_appFWork_pfr_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,1);
mean_nd_appFWork_pfr_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,1);
mean_nd_appFWork_pfr_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,1);

%Mean of tracking error per full run (tepfr)
mean_nd_tepfr_LSF_p5_LE_0 = zeros(NumOf_LSF_p5_LE_0_Files,1);
mean_nd_tepfr_LSF_p5_LE_p2 = zeros(NumOf_LSF_p5_LE_p2_Files,1);
mean_nd_tepfr_LSF_p5_LE_p4 = zeros(NumOf_LSF_p5_LE_p4_Files,1);
mean_nd_tepfr_LSF_1_LE_0 = zeros(NumOf_LSF_1_LE_0_Files,1);
mean_nd_tepfr_LSF_1_LE_p2 = zeros(NumOf_LSF_1_LE_p2_Files,1);
mean_nd_tepfr_LSF_1_LE_p4 = zeros(NumOf_LSF_1_LE_p4_Files,1);
mean_nd_tepfr_LSF_2_LE_0 = zeros(NumOf_LSF_2_LE_0_Files,1);
mean_nd_tepfr_LSF_2_LE_p2 = zeros(NumOf_LSF_2_LE_p2_Files,1);
mean_nd_tepfr_LSF_2_LE_p4 = zeros(NumOf_LSF_2_LE_p4_Files,1);

%Calculate the means and standard deviation of each work treatment
for i = 1:PartPath
    mean_work_LSF_p5_LE_0(1,i) = nanmean(work_LSF_p5_LE_0(:,i));
    mean_work_LSF_p5_LE_p2(1,i) = nanmean(work_LSF_p5_LE_p2(:,i));
    mean_work_LSF_p5_LE_p4(1,i) = nanmean(work_LSF_p5_LE_p4(:,i));
    mean_work_LSF_1_LE_0(1,i) = nanmean(work_LSF_1_LE_0(:,i));
    mean_work_LSF_1_LE_p2(1,i) = nanmean(work_LSF_1_LE_p2(:,i));
    mean_work_LSF_1_LE_p4(1,i) = nanmean(work_LSF_1_LE_p4(:,i));
    mean_work_LSF_2_LE_0(1,i) = nanmean(work_LSF_2_LE_0(:,i));
    mean_work_LSF_2_LE_p2(1,i) = nanmean(work_LSF_2_LE_p2(:,i));
    mean_work_LSF_2_LE_p4(1,i) = nanmean(work_LSF_2_LE_p4(:,i));
    
    std_work_LSF_p5_LE_0(1,i) = nanstd(work_LSF_p5_LE_0(:,i));
    std_work_LSF_p5_LE_p2(1,i) = nanstd(work_LSF_p5_LE_p2(:,i));
    std_work_LSF_p5_LE_p4(1,i) = nanstd(work_LSF_p5_LE_p4(:,i));
    std_work_LSF_1_LE_0(1,i) = nanstd(work_LSF_1_LE_0(:,i));
    std_work_LSF_1_LE_p2(1,i) = nanstd(work_LSF_1_LE_p2(:,i));
    std_work_LSF_1_LE_p4(1,i) = nanstd(work_LSF_1_LE_p4(:,i));
    std_work_LSF_2_LE_0(1,i) = nanstd(work_LSF_2_LE_0(:,i));
    std_work_LSF_2_LE_p2(1,i) = nanstd(work_LSF_2_LE_p2(:,i));
    std_work_LSF_2_LE_p4(1,i) = nanstd(work_LSF_2_LE_p4(:,i));
    
end

%wpfr means "work per full run"
for i = 1:NumOf_LSF_p5_LE_0_Files
    mean_wpfr_LSF_p5_LE_0(i,1) = nanmean(work_LSF_p5_LE_0(i,:));
    mean_nd_wpfr_LSF_p5_LE_0(i,1) = nanmean(work_nd_LSF_p5_LE_0(i,:));
    mean_nd_tepfr_LSF_p5_LE_0(i,1) = nanmean(nd_te_LSF_p5_LE_0(i,2:end));
    
    mean_wingWork_pfr_LSF_p5_LE_0(i,1) = nanmean(wingWork_LSF_p5_LE_0(i,:));
    mean_abdoWork_pfr_LSF_p5_LE_0(i,1) = nanmean(abdoWork_LSF_p5_LE_0(i,:));
    mean_appFWork_pfr_LSF_p5_LE_0(i,1) = nanmean(appFWork_LSF_p5_LE_0(i,:));
    
    mean_nd_wingWork_pfr_LSF_p5_LE_0(i,1) = nanmean(wingWork_nd_LSF_p5_LE_0(i,:));
    mean_nd_abdoWork_pfr_LSF_p5_LE_0(i,1) = nanmean(abdoWork_nd_LSF_p5_LE_0(i,:));
    mean_nd_appFWork_pfr_LSF_p5_LE_0(i,1) = nanmean(appFWork_nd_LSF_p5_LE_0(i,:));
end

for i = 1:NumOf_LSF_p5_LE_p2_Files
    mean_wpfr_LSF_p5_LE_p2(i,1) = nanmean(work_LSF_p5_LE_p2(i,:));
    mean_nd_wpfr_LSF_p5_LE_p2(i,1) = nanmean(work_nd_LSF_p5_LE_p2(i,:));
    mean_nd_tepfr_LSF_p5_LE_p2(i,1) = nanmean(nd_te_LSF_p5_LE_p2(i,2:end));
    
    mean_wingWork_pfr_LSF_p5_LE_p2(i,1) = nanmean(wingWork_LSF_p5_LE_p2(i,:));
    mean_abdoWork_pfr_LSF_p5_LE_p2(i,1) = nanmean(abdoWork_LSF_p5_LE_p2(i,:));
    mean_appFWork_pfr_LSF_p5_LE_p2(i,1) = nanmean(appFWork_LSF_p5_LE_p2(i,:));
    
    mean_nd_wingWork_pfr_LSF_p5_LE_p2(i,1) = nanmean(wingWork_nd_LSF_p5_LE_p2(i,:));
    mean_nd_abdoWork_pfr_LSF_p5_LE_p2(i,1) = nanmean(abdoWork_nd_LSF_p5_LE_p2(i,:));
    mean_nd_appFWork_pfr_LSF_p5_LE_p2(i,1) = nanmean(appFWork_nd_LSF_p5_LE_p2(i,:));
end

for i = 1:NumOf_LSF_p5_LE_p4_Files
    mean_wpfr_LSF_p5_LE_p4(i,1) = nanmean(work_LSF_p5_LE_p4(i,:));
    mean_nd_wpfr_LSF_p5_LE_p4(i,1) = nanmean(work_nd_LSF_p5_LE_p4(i,:));
    mean_nd_tepfr_LSF_p5_LE_p4(i,1) = nanmean(nd_te_LSF_p5_LE_p4(i,2:end));
    
    mean_wingWork_pfr_LSF_p5_LE_p4(i,1) = nanmean(wingWork_LSF_p5_LE_p4(i,:));
    mean_abdoWork_pfr_LSF_p5_LE_p4(i,1) = nanmean(abdoWork_LSF_p5_LE_p4(i,:));
    mean_appFWork_pfr_LSF_p5_LE_p4(i,1) = nanmean(appFWork_LSF_p5_LE_p4(i,:));
    
    mean_nd_wingWork_pfr_LSF_p5_LE_p4(i,1) = nanmean(wingWork_nd_LSF_p5_LE_p4(i,:));
    mean_nd_abdoWork_pfr_LSF_p5_LE_p4(i,1) = nanmean(abdoWork_nd_LSF_p5_LE_p4(i,:));
    mean_nd_appFWork_pfr_LSF_p5_LE_p4(i,1) = nanmean(appFWork_nd_LSF_p5_LE_p4(i,:));
end

for i = 1:NumOf_LSF_1_LE_0_Files
    mean_wpfr_LSF_1_LE_0(i,1) = nanmean(work_LSF_1_LE_0(i,:));
    mean_nd_wpfr_LSF_1_LE_0(i,1) = nanmean(work_nd_LSF_1_LE_0(i,:));
    mean_nd_tepfr_LSF_1_LE_0(i,1) = nanmean(nd_te_LSF_1_LE_0(i,2:end));
    
    mean_wingWork_pfr_LSF_1_LE_0(i,1) = nanmean(wingWork_LSF_1_LE_0(i,:));
    mean_abdoWork_pfr_LSF_1_LE_0(i,1) = nanmean(abdoWork_LSF_1_LE_0(i,:));
    mean_appFWork_pfr_LSF_1_LE_0(i,1) = nanmean(appFWork_LSF_1_LE_0(i,:));
    
    mean_nd_wingWork_pfr_LSF_1_LE_0(i,1) = nanmean(wingWork_nd_LSF_1_LE_0(i,:));
    mean_nd_abdoWork_pfr_LSF_1_LE_0(i,1) = nanmean(abdoWork_nd_LSF_1_LE_0(i,:));
    mean_nd_appFWork_pfr_LSF_1_LE_0(i,1) = nanmean(appFWork_nd_LSF_1_LE_0(i,:));
end

for i = 1:NumOf_LSF_1_LE_p2_Files
    mean_wpfr_LSF_1_LE_p2(i,1) = nanmean(work_LSF_1_LE_p2(i,:));
    mean_nd_wpfr_LSF_1_LE_p2(i,1) = nanmean(work_nd_LSF_1_LE_p2(i,:));
    mean_nd_tepfr_LSF_1_LE_p2(i,1) = nanmean(nd_te_LSF_1_LE_p2(i,2:end));
    
    mean_wingWork_pfr_LSF_1_LE_p2(i,1) = nanmean(wingWork_LSF_1_LE_p2(i,:));
    mean_abdoWork_pfr_LSF_1_LE_p2(i,1) = nanmean(abdoWork_LSF_1_LE_p2(i,:));
    mean_appFWork_pfr_LSF_1_LE_p2(i,1) = nanmean(appFWork_LSF_1_LE_p2(i,:));
    
    mean_nd_wingWork_pfr_LSF_1_LE_p2(i,1) = nanmean(wingWork_nd_LSF_1_LE_p2(i,:));
    mean_nd_abdoWork_pfr_LSF_1_LE_p2(i,1) = nanmean(abdoWork_nd_LSF_1_LE_p2(i,:));
    mean_nd_appFWork_pfr_LSF_1_LE_p2(i,1) = nanmean(appFWork_nd_LSF_1_LE_p2(i,:));
end

for i = 1:NumOf_LSF_1_LE_p4_Files
    mean_wpfr_LSF_1_LE_p4(i,1) = nanmean(work_LSF_1_LE_p4(i,:));
    mean_nd_wpfr_LSF_1_LE_p4(i,1) = nanmean(work_nd_LSF_1_LE_p4(i,:));
    mean_nd_tepfr_LSF_1_LE_p4(i,1) = nanmean(nd_te_LSF_1_LE_p4(i,2:end));
    
    mean_wingWork_pfr_LSF_1_LE_p4(i,1) = nanmean(wingWork_LSF_1_LE_p4(i,:));
    mean_abdoWork_pfr_LSF_1_LE_p4(i,1) = nanmean(abdoWork_LSF_1_LE_p4(i,:));
    mean_appFWork_pfr_LSF_1_LE_p4(i,1) = nanmean(appFWork_LSF_1_LE_p4(i,:));
    
    mean_nd_wingWork_pfr_LSF_1_LE_p4(i,1) = nanmean(wingWork_nd_LSF_1_LE_p4(i,:));
    mean_nd_abdoWork_pfr_LSF_1_LE_p4(i,1) = nanmean(abdoWork_nd_LSF_1_LE_p4(i,:));
    mean_nd_appFWork_pfr_LSF_1_LE_p4(i,1) = nanmean(appFWork_nd_LSF_1_LE_p4(i,:));
end

for i = 1:NumOf_LSF_2_LE_0_Files
    mean_wpfr_LSF_2_LE_0(i,1) = nanmean(work_LSF_2_LE_0(i,:));
    mean_nd_wpfr_LSF_2_LE_0(i,1) = nanmean(work_nd_LSF_2_LE_0(i,:));
    mean_nd_tepfr_LSF_2_LE_0(i,1) = nanmean(nd_te_LSF_2_LE_0(i,2:end));
    
    mean_wingWork_pfr_LSF_2_LE_0(i,1) = nanmean(wingWork_LSF_2_LE_0(i,:));
    mean_abdoWork_pfr_LSF_2_LE_0(i,1) = nanmean(abdoWork_LSF_2_LE_0(i,:));
    mean_appFWork_pfr_LSF_2_LE_0(i,1) = nanmean(appFWork_LSF_2_LE_0(i,:));
    
    mean_nd_wingWork_pfr_LSF_2_LE_0(i,1) = nanmean(wingWork_nd_LSF_2_LE_0(i,:));
    mean_nd_abdoWork_pfr_LSF_2_LE_0(i,1) = nanmean(abdoWork_nd_LSF_2_LE_0(i,:));
    mean_nd_appFWork_pfr_LSF_2_LE_0(i,1) = nanmean(appFWork_nd_LSF_2_LE_0(i,:));
end

for i = 1:NumOf_LSF_2_LE_p2_Files
    mean_wpfr_LSF_2_LE_p2(i,1) = nanmean(work_LSF_2_LE_p2(i,:));
    mean_nd_wpfr_LSF_2_LE_p2(i,1) = nanmean(work_nd_LSF_2_LE_p2(i,:));
    mean_nd_tepfr_LSF_2_LE_p2(i,1) = nanmean(nd_te_LSF_2_LE_p2(i,2:end));
    
    mean_wingWork_pfr_LSF_2_LE_p2(i,1) = nanmean(wingWork_LSF_2_LE_p2(i,:));
    mean_abdoWork_pfr_LSF_2_LE_p2(i,1) = nanmean(abdoWork_LSF_2_LE_p2(i,:));
    mean_appFWork_pfr_LSF_2_LE_p2(i,1) = nanmean(appFWork_LSF_2_LE_p2(i,:));
    
    mean_nd_wingWork_pfr_LSF_2_LE_p2(i,1) = nanmean(wingWork_nd_LSF_2_LE_p2(i,:));
    mean_nd_abdoWork_pfr_LSF_2_LE_p2(i,1) = nanmean(abdoWork_nd_LSF_2_LE_p2(i,:));
    mean_nd_appFWork_pfr_LSF_2_LE_p2(i,1) = nanmean(appFWork_nd_LSF_2_LE_p2(i,:));
end

for i = 1:NumOf_LSF_2_LE_p4_Files
    mean_wpfr_LSF_2_LE_p4(i,1) = nanmean(work_LSF_2_LE_p4(i,:));
    mean_nd_wpfr_LSF_2_LE_p4(i,1) = nanmean(work_nd_LSF_2_LE_p4(i,:));
    mean_nd_tepfr_LSF_2_LE_p4(i,1) = nanmean(nd_te_LSF_2_LE_p4(i,2:end));
    
    mean_wingWork_pfr_LSF_2_LE_p4(i,1) = nanmean(wingWork_LSF_2_LE_p4(i,:));
    mean_abdoWork_pfr_LSF_2_LE_p4(i,1) = nanmean(abdoWork_LSF_2_LE_p4(i,:));
    mean_appFWork_pfr_LSF_2_LE_p4(i,1) = nanmean(appFWork_LSF_2_LE_p4(i,:));
    
    mean_nd_wingWork_pfr_LSF_2_LE_p4(i,1) = nanmean(wingWork_nd_LSF_2_LE_p4(i,:));
    mean_nd_abdoWork_pfr_LSF_2_LE_p4(i,1) = nanmean(abdoWork_nd_LSF_2_LE_p4(i,:));
    mean_nd_appFWork_pfr_LSF_2_LE_p4(i,1) = nanmean(appFWork_nd_LSF_2_LE_p4(i,:));
end

%% Generate the stats files

%Curate for Stats csv file
size_export_LSF_p5_LE_0(1:numel(mean_wpfr_LSF_p5_LE_0),1) = {'LSF_p5'};
shape_export_LSF_p5_LE_0(1:numel(mean_wpfr_LSF_p5_LE_0),1) = {'LE_0'};
task_export_LSF_p5_LE_0(1:numel(mean_wpfr_LSF_p5_LE_0),1) = {'vertical'};
size_export_LSF_p5_LE_p2(1:numel(mean_wpfr_LSF_p5_LE_p2),1) = {'LSF_p5'};
shape_export_LSF_p5_LE_p2(1:numel(mean_wpfr_LSF_p5_LE_p2),1) = {'LE_p2'};
task_export_LSF_p5_LE_p2(1:numel(mean_wpfr_LSF_p5_LE_p2),1) = {'vertical'};
size_export_LSF_p5_LE_p4(1:numel(mean_wpfr_LSF_p5_LE_p4),1) = {'LSF_p5'};
shape_export_LSF_p5_LE_p4(1:numel(mean_wpfr_LSF_p5_LE_p4),1) = {'LE_p4'};
task_export_LSF_p5_LE_p4(1:numel(mean_wpfr_LSF_p5_LE_p4),1) = {'vertical'};
size_export_LSF_1_LE_0(1:numel(mean_wpfr_LSF_1_LE_0),1) = {'LSF_1'};
shape_export_LSF_1_LE_0(1:numel(mean_wpfr_LSF_1_LE_0),1) = {'LE_0'};
task_export_LSF_1_LE_0(1:numel(mean_wpfr_LSF_1_LE_0),1) = {'vertical'};
size_export_LSF_1_LE_p2(1:numel(mean_wpfr_LSF_1_LE_p2),1) = {'LSF_1'};
shape_export_LSF_1_LE_p2(1:numel(mean_wpfr_LSF_1_LE_p2),1) = {'LE_p2'};
task_export_LSF_1_LE_p2(1:numel(mean_wpfr_LSF_1_LE_p2),1) = {'vertical'};
size_export_LSF_1_LE_p4(1:numel(mean_wpfr_LSF_1_LE_p4),1) = {'LSF_1'};
shape_export_LSF_1_LE_p4(1:numel(mean_wpfr_LSF_1_LE_p4),1) = {'LE_p4'};
task_export_LSF_1_LE_p4(1:numel(mean_wpfr_LSF_1_LE_p4),1) = {'vertical'};
size_export_LSF_2_LE_0(1:numel(mean_wpfr_LSF_2_LE_0),1) = {'LSF_2'};
shape_export_LSF_2_LE_0(1:numel(mean_wpfr_LSF_2_LE_0),1) = {'LE_0'};
task_export_LSF_2_LE_0(1:numel(mean_wpfr_LSF_2_LE_0),1) = {'vertical'};
size_export_LSF_2_LE_p2(1:numel(mean_wpfr_LSF_2_LE_p2),1) = {'LSF_2'};
shape_export_LSF_2_LE_p2(1:numel(mean_wpfr_LSF_2_LE_p2),1) = {'LE_p2'};
task_export_LSF_2_LE_p2(1:numel(mean_wpfr_LSF_2_LE_p2),1) = {'vertical'};
size_export_LSF_2_LE_p4(1:numel(mean_wpfr_LSF_2_LE_p4),1) = {'LSF_2'};
shape_export_LSF_2_LE_p4(1:numel(mean_wpfr_LSF_2_LE_p4),1) = {'LE_p4'};
task_export_LSF_2_LE_p4(1:numel(mean_wpfr_LSF_2_LE_p4),1) = {'vertical'};

%sish below means "SIze" and "SHape"
sishManip_export_LSF_p5_LE_0(1:numel(mean_wpfr_LSF_p5_LE_0),1) = {'L05_E0'}; %L05_E0
sishManip_export_LSF_p5_LE_p2(1:numel(mean_wpfr_LSF_p5_LE_p2),1) = {'L05_Ep2'}; %L05_Ep2
sishManip_export_LSF_p5_LE_p4(1:numel(mean_wpfr_LSF_p5_LE_p4),1) = {'L05_Ep4'}; %L05_Ep4
sishManip_export_LSF_1_LE_0(1:numel(mean_wpfr_LSF_1_LE_0),1) = {'L10_E0'}; %L10_E0
sishManip_export_LSF_1_LE_p2(1:numel(mean_wpfr_LSF_1_LE_p2),1) = {'L10_Ep2'}; %L10_Ep2
sishManip_export_LSF_1_LE_p4(1:numel(mean_wpfr_LSF_1_LE_p4),1) = {'L10_Ep4'}; %L10_Ep4
sishManip_export_LSF_2_LE_0(1:numel(mean_wpfr_LSF_2_LE_0),1) = {'L20_E0'}; %L20_E0
sishManip_export_LSF_2_LE_p2(1:numel(mean_wpfr_LSF_2_LE_p2),1) = {'L20_Ep2'}; %L20_Ep2
sishManip_export_LSF_2_LE_p4(1:numel(mean_wpfr_LSF_2_LE_p4),1) = {'L20_Ep4'}; %L20_Ep4

export_size = [size_export_LSF_p5_LE_0; size_export_LSF_p5_LE_p2; ...
    size_export_LSF_p5_LE_p4; size_export_LSF_1_LE_0; ...
    size_export_LSF_1_LE_p2; size_export_LSF_1_LE_p4;...
    size_export_LSF_2_LE_0; size_export_LSF_2_LE_p2; ...
    size_export_LSF_2_LE_p4];

export_shape = [shape_export_LSF_p5_LE_0; shape_export_LSF_p5_LE_p2; ...
    shape_export_LSF_p5_LE_p4; shape_export_LSF_1_LE_0; ...
    shape_export_LSF_1_LE_p2; shape_export_LSF_1_LE_p4;...
    shape_export_LSF_2_LE_0; shape_export_LSF_2_LE_p2; ...
    shape_export_LSF_2_LE_p4];

export_task = [task_export_LSF_p5_LE_0; task_export_LSF_p5_LE_p2; ...
    task_export_LSF_p5_LE_p4; task_export_LSF_1_LE_0; ...
    task_export_LSF_1_LE_p2; task_export_LSF_1_LE_p4;...
    task_export_LSF_2_LE_0; task_export_LSF_2_LE_p2;...
    task_export_LSF_2_LE_p4];

export_sishManip = [sishManip_export_LSF_p5_LE_0;...
    sishManip_export_LSF_p5_LE_p2; sishManip_export_LSF_p5_LE_p4;...
    sishManip_export_LSF_1_LE_0; sishManip_export_LSF_1_LE_p2;...
    sishManip_export_LSF_1_LE_p4; sishManip_export_LSF_2_LE_0;...
    sishManip_export_LSF_2_LE_p2; sishManip_export_LSF_2_LE_p4];

%Work per full run (in Joules)
wpfr_export_numbers = [mean_wpfr_LSF_p5_LE_0; mean_wpfr_LSF_p5_LE_p2;...
    mean_wpfr_LSF_p5_LE_p4; mean_wpfr_LSF_1_LE_0; mean_wpfr_LSF_1_LE_p2;...
    mean_wpfr_LSF_1_LE_p4; mean_wpfr_LSF_2_LE_0; mean_wpfr_LSF_2_LE_p2;...
    mean_wpfr_LSF_2_LE_p4];

wpfr_export_numbers = 1000.*wpfr_export_numbers; %convert to milliJoules

%Non-dimensional work per full run
nd_wpfr_export_numbers = [mean_nd_wpfr_LSF_p5_LE_0;...
    mean_nd_wpfr_LSF_p5_LE_p2; mean_nd_wpfr_LSF_p5_LE_p4;...
    mean_nd_wpfr_LSF_1_LE_0; mean_nd_wpfr_LSF_1_LE_p2;...
    mean_nd_wpfr_LSF_1_LE_p4; mean_nd_wpfr_LSF_2_LE_0;...
    mean_nd_wpfr_LSF_2_LE_p2; mean_nd_wpfr_LSF_2_LE_p4];

%Non-dimensional tracking error per full run
nd_tepfr_export_numbers = [mean_nd_tepfr_LSF_p5_LE_0;...
    mean_nd_tepfr_LSF_p5_LE_p2; mean_nd_tepfr_LSF_p5_LE_p4;...
    mean_nd_tepfr_LSF_1_LE_0; mean_nd_tepfr_LSF_1_LE_p2;...
    mean_nd_tepfr_LSF_1_LE_p4; mean_nd_tepfr_LSF_2_LE_0;...
    mean_nd_tepfr_LSF_2_LE_p2; mean_nd_tepfr_LSF_2_LE_p4;];

%Cost of transport * non-dimensional tracking error per full run
nd_product_export_numbers = [mean_nd_wpfr_LSF_p5_LE_0.*mean_nd_tepfr_LSF_p5_LE_0;...
    mean_nd_wpfr_LSF_p5_LE_p2.*mean_nd_tepfr_LSF_p5_LE_p2;...
    mean_nd_wpfr_LSF_p5_LE_p4.*mean_nd_tepfr_LSF_p5_LE_p4;...
    mean_nd_wpfr_LSF_1_LE_0.*mean_nd_tepfr_LSF_1_LE_0;...
    mean_nd_wpfr_LSF_1_LE_p2.*mean_nd_tepfr_LSF_1_LE_p2;...
    mean_nd_wpfr_LSF_1_LE_p4.*mean_nd_tepfr_LSF_1_LE_p4;...
    mean_nd_wpfr_LSF_2_LE_0.*mean_nd_tepfr_LSF_2_LE_0;...
    mean_nd_wpfr_LSF_2_LE_p2.*mean_nd_tepfr_LSF_2_LE_p2;...
    mean_nd_wpfr_LSF_2_LE_p4.*mean_nd_tepfr_LSF_2_LE_p4];

if strcmp(str1,'y')
    %Work per full run
    wpfr_write = fopen('../StatsStuff_MPC/sishManip_wpfr_vert.csv','w');
    %Write the names to a CSV file
    for k = 1:numel(export_size)
        fprintf(wpfr_write,'%f,',wpfr_export_numbers(k,1)); 
                        %The %f is floating point precision number,
                        %the \n indicates a new line
        fprintf(wpfr_write,'%s,', export_size{k,1}); 
                        %The %s is for string formatting,
                        %Make sure the comma is where it is!
        fprintf(wpfr_write,'%s,', export_shape{k,1}); 
                        %The %s is for string formatting,
                        %Make sure the comma is where it is!
        fprintf(wpfr_write,'%s,', export_task{k,1}); 
                        %The %s is for string formatting,
                        %Make sure the comma is where it is!
        fprintf(wpfr_write,'%s\n', export_sishManip{k,1}); 
                        %The %s is for string formatting,
                        %Make sure the comma is where it is!
    end
    fclose(wpfr_write);

    %Non-dimensionalized work per full run
    nd_wpfr_write = fopen('../StatsStuff_MPC/sishManip_nd_wpfr_vert.csv','w');
    %Write the names to a CSV file
    for k = 1:numel(export_size)
        fprintf(nd_wpfr_write,'%f,',nd_wpfr_export_numbers(k,1)); 
                        %The %f is floating point precision number,
                        %the \n indicates a new line
        fprintf(nd_wpfr_write,'%s,', export_size{k,1}); 
                        %The %s is for string formatting,
                        %Make sure the comma is where it is!
        fprintf(nd_wpfr_write,'%s,', export_shape{k,1}); 
                        %The %s is for string formatting,
                        %Make sure the comma is where it is!
        fprintf(nd_wpfr_write,'%s,', export_task{k,1}); 
                        %The %s is for string formatting,
                        %Make sure the comma is where it is!
        fprintf(nd_wpfr_write,'%s\n', export_sishManip{k,1}); 
                        %The %s is for string formatting,
                        %Make sure the comma is where it is!
    end
    fclose(nd_wpfr_write);

    %Non-dimensionalized tracking error per full run
    nd_tepfr_write = fopen('../StatsStuff_MPC/sishManip_nd_tepfr_vert.csv','w');
    %Write the names to a CSV file
    for k = 1:numel(export_size)
        fprintf(nd_tepfr_write,'%f,',nd_tepfr_export_numbers(k,1)); 
                        %The %f is floating point precision number,
                        %the \n indicates a new line
        fprintf(nd_tepfr_write,'%s,', export_size{k,1}); 
                        %The %s is for string formatting,
                        %Make sure the comma is where it is!
        fprintf(nd_tepfr_write,'%s,', export_shape{k,1}); 
                        %The %s is for string formatting,
                        %Make sure the comma is where it is!
        fprintf(nd_tepfr_write,'%s,', export_task{k,1}); 
                        %The %s is for string formatting,
                        %Make sure the comma is where it is!
        fprintf(nd_tepfr_write,'%s\n', export_sishManip{k,1}); 
                        %The %s is for string formatting,
                        %Make sure the comma is where it is!
    end
    fclose(nd_tepfr_write);

    %Product of non-dimensionalize work and non-dimensionalized tracking error 
        %per full run
    nd_product_write = fopen('../StatsStuff_MPC/sishManip_nd_product_vert.csv','w');
    %Write the names to a CSV file
    for k = 1:numel(export_size)
        fprintf(nd_product_write,'%f,',nd_product_export_numbers(k,1)); 
                        %The %f is floating point precision number,
                        %the \n indicates a new line
        fprintf(nd_product_write,'%s,', export_size{k,1}); 
                        %The %s is for string formatting,
                        %Make sure the comma is where it is!
        fprintf(nd_product_write,'%s,', export_shape{k,1}); 
                        %The %s is for string formatting,
                        %Make sure the comma is where it is!
        fprintf(nd_product_write,'%s,', export_task{k,1}); 
                        %The %s is for string formatting,
                        %Make sure the comma is where it is!
        fprintf(nd_product_write,'%s\n', export_sishManip{k,1}); 
                        %The %s is for string formatting,
                        %Make sure the comma is where it is!
    end
    fclose(nd_product_write);
end
%% Find the min & max values for work

workMax_vec = [max(max(mean_work_LSF_p5_LE_0 + std_work_LSF_p5_LE_0)),...
    max(max(mean_work_LSF_p5_LE_p2 + std_work_LSF_p5_LE_p2)),...
    max(max(mean_work_LSF_p5_LE_p4 + std_work_LSF_p5_LE_p4)),...
    max(max(mean_work_LSF_1_LE_0 + std_work_LSF_1_LE_0)),...
    max(max(mean_work_LSF_1_LE_p2 + std_work_LSF_1_LE_p2)),...
    max(max(mean_work_LSF_1_LE_p4 + std_work_LSF_1_LE_p4)),...
    max(max(mean_work_LSF_2_LE_0 + std_work_LSF_2_LE_0)),...
    max(max(mean_work_LSF_2_LE_p2 + std_work_LSF_2_LE_p2)),...
    max(max(mean_work_LSF_2_LE_p4 + std_work_LSF_2_LE_p4))];
work_max = max(workMax_vec);

workMin_vec = [min(min(mean_work_LSF_p5_LE_0 - std_work_LSF_p5_LE_0)),...
    min(min(mean_work_LSF_p5_LE_p2 - std_work_LSF_p5_LE_p2)),...
    min(min(mean_work_LSF_p5_LE_p4 - std_work_LSF_p5_LE_p4)),...
    min(min(mean_work_LSF_1_LE_0 - std_work_LSF_1_LE_0)),...
    min(min(mean_work_LSF_1_LE_p2 - std_work_LSF_1_LE_p2)),...
    min(min(mean_work_LSF_1_LE_p4 - std_work_LSF_1_LE_p4)),...
    min(min(mean_work_LSF_2_LE_0 - std_work_LSF_2_LE_0)),...
    min(min(mean_work_LSF_2_LE_p2 - std_work_LSF_2_LE_p2)),...
    min(min(mean_work_LSF_2_LE_p4 - std_work_LSF_2_LE_p4))];
work_min = min(workMin_vec);

%wpfr max and min
meanworkMax_vec = [max(max(mean_wpfr_LSF_p5_LE_0)),...
    max(max(mean_wpfr_LSF_p5_LE_p2)),...
    max(max(mean_wpfr_LSF_p5_LE_p4)),...
    max(max(mean_wpfr_LSF_1_LE_0)),...
    max(max(mean_wpfr_LSF_1_LE_p2)),...
    max(max(mean_wpfr_LSF_1_LE_p4)),...
    max(max(mean_wpfr_LSF_2_LE_0)),...
    max(max(mean_wpfr_LSF_2_LE_p2)),...
    max(max(mean_wpfr_LSF_2_LE_p4))];
meanwork_max = max(meanworkMax_vec);

meanworkMin_vec = [min(min(mean_wpfr_LSF_p5_LE_0)),...
    min(min(mean_wpfr_LSF_p5_LE_p2)),...
    min(min(mean_wpfr_LSF_p5_LE_p4)),...
    min(min(mean_wpfr_LSF_1_LE_0)),...
    min(min(mean_wpfr_LSF_1_LE_p2)),...
    min(min(mean_wpfr_LSF_1_LE_p4)),...
    min(min(mean_wpfr_LSF_2_LE_0)),...
    min(min(mean_wpfr_LSF_2_LE_p2)),...
    min(min(mean_wpfr_LSF_2_LE_p4))];
meanwork_min = min(meanworkMin_vec);

%non-dimensionalized wpfr max and min
meanworkMax_nd_vec = [max(max(mean_nd_wpfr_LSF_p5_LE_0)),...
    max(max(mean_nd_wpfr_LSF_p5_LE_p2)),...
    max(max(mean_nd_wpfr_LSF_p5_LE_p4)),...
    max(max(mean_nd_wpfr_LSF_1_LE_0)),...
    max(max(mean_nd_wpfr_LSF_1_LE_p2)),...
    max(max(mean_nd_wpfr_LSF_1_LE_p4)),...
    max(max(mean_nd_wpfr_LSF_2_LE_0)),...
    max(max(mean_nd_wpfr_LSF_2_LE_p2)),...
    max(max(mean_nd_wpfr_LSF_2_LE_p4))];
meanwork_nd_max = max(meanworkMax_nd_vec);

meanworkMin_nd_vec = [min(min(mean_nd_wpfr_LSF_p5_LE_0)),...
    min(min(mean_nd_wpfr_LSF_p5_LE_p2)),...
    min(min(mean_nd_wpfr_LSF_p5_LE_p4)),...
    min(min(mean_nd_wpfr_LSF_1_LE_0)),...
    min(min(mean_nd_wpfr_LSF_1_LE_p2)),...
    min(min(mean_nd_wpfr_LSF_1_LE_p4)),...
    min(min(mean_nd_wpfr_LSF_2_LE_0)),...
    min(min(mean_nd_wpfr_LSF_2_LE_p2)),...
    min(min(mean_nd_wpfr_LSF_2_LE_p4))];
meanwork_nd_min = min(meanworkMin_nd_vec);

%non-dimensionalized tepfr max and min
meantepfrMax_nd_vec = [max(max(mean_nd_tepfr_LSF_p5_LE_0)),...
    max(max(mean_nd_tepfr_LSF_p5_LE_p2)),...
    max(max(mean_nd_tepfr_LSF_p5_LE_p4)),...
    max(max(mean_nd_tepfr_LSF_1_LE_0)),...
    max(max(mean_nd_tepfr_LSF_1_LE_p2)),...
    max(max(mean_nd_tepfr_LSF_1_LE_p4)),...
    max(max(mean_nd_tepfr_LSF_2_LE_0)),...
    max(max(mean_nd_tepfr_LSF_2_LE_p2)),...
    max(max(mean_nd_tepfr_LSF_2_LE_p4))];
meantepfr_nd_max = max(meantepfrMax_nd_vec);

meantepfrMin_nd_vec = [min(min(mean_nd_tepfr_LSF_p5_LE_0)),...
    min(min(mean_nd_tepfr_LSF_p5_LE_p2)),...
    min(min(mean_nd_tepfr_LSF_p5_LE_p4)),...
    min(min(mean_nd_tepfr_LSF_1_LE_0)),...
    min(min(mean_nd_tepfr_LSF_1_LE_p2)),...
    min(min(mean_nd_tepfr_LSF_1_LE_p4)),...
    min(min(mean_nd_tepfr_LSF_2_LE_0)),...
    min(min(mean_nd_tepfr_LSF_2_LE_p2)),...
    min(min(mean_nd_tepfr_LSF_2_LE_p4))];
meantepfr_nd_min = min(meantepfrMin_nd_vec);

%non-dimensionalized product of work*tracking error max and min
meanProductMax_vec = [max(max(mean_nd_wpfr_LSF_p5_LE_0.*mean_nd_tepfr_LSF_p5_LE_0)),...
    max(max(mean_nd_wpfr_LSF_p5_LE_p2.*mean_nd_tepfr_LSF_p5_LE_p2)),...
    max(max(mean_nd_wpfr_LSF_p5_LE_p4.*mean_nd_tepfr_LSF_p5_LE_p4)),...
    max(max(mean_nd_wpfr_LSF_1_LE_0.*mean_nd_tepfr_LSF_1_LE_0)),...
    max(max(mean_nd_wpfr_LSF_1_LE_p2.*mean_nd_tepfr_LSF_1_LE_p2)),...
    max(max(mean_nd_wpfr_LSF_1_LE_p4.*mean_nd_tepfr_LSF_1_LE_p4)),...
    max(max(mean_nd_wpfr_LSF_2_LE_0.*mean_nd_tepfr_LSF_2_LE_0)),...
    max(max(mean_nd_wpfr_LSF_2_LE_p2.*mean_nd_tepfr_LSF_2_LE_p2)),...
    max(max(mean_nd_wpfr_LSF_2_LE_p4.*mean_nd_tepfr_LSF_2_LE_p4))];
meanProduct_max = max(meanProductMax_vec);

meanProductMin_vec = [min(min(mean_nd_wpfr_LSF_p5_LE_0.*mean_nd_tepfr_LSF_p5_LE_0)),...
    min(min(mean_nd_wpfr_LSF_p5_LE_p2.*mean_nd_tepfr_LSF_p5_LE_p2)),...
    min(min(mean_nd_wpfr_LSF_p5_LE_p4.*mean_nd_tepfr_LSF_p5_LE_p4)),...
    min(min(mean_nd_wpfr_LSF_1_LE_0.*mean_nd_tepfr_LSF_1_LE_0)),...
    min(min(mean_nd_wpfr_LSF_1_LE_p2.*mean_nd_tepfr_LSF_1_LE_p2)),...
    min(min(mean_nd_wpfr_LSF_1_LE_p4.*mean_nd_tepfr_LSF_1_LE_p4)),...
    min(min(mean_nd_wpfr_LSF_2_LE_0.*mean_nd_tepfr_LSF_2_LE_0)),...
    min(min(mean_nd_wpfr_LSF_2_LE_p2.*mean_nd_tepfr_LSF_2_LE_p2)),...
    min(min(mean_nd_wpfr_LSF_2_LE_p4.*mean_nd_tepfr_LSF_2_LE_p4))];
meanProduct_min = min(meanProductMin_vec);
stop
%% Plot the work individualized

fig1 = figure(1);
subplot(3,3,1)
if NumOf_LSF_p5_LE_0_Files > 0
    plot(time_hws,work_LSF_p5_LE_0,'LineWidth',2)
end
title('LSF: 0.5, LE: 0')
% axis([0, time_hws(end), 0, 0.06])
ylabel('Work (Joules)')

subplot(3,3,2)
if NumOf_LSF_p5_LE_p2_Files > 0
    plot(time_hws,work_LSF_p5_LE_p2,'LineWidth',2)
end
title('LSF: 0.5, LE: 0.2')
% axis([0, time_hws(end), 0, 0.06])
% ylabel({'fs';'Work (Joules)'})

subplot(3,3,3)
if NumOf_LSF_p5_LE_p4_Files > 0
    plot(time_hws,work_LSF_p5_LE_p4,'LineWidth',2)
end
title('LSF: 0.5, LE: 0.4')
% axis([0, time_hws(end), 0, 0.06])
% ylabel({'ua';'Work (Joules)'})

subplot(3,3,4)
if NumOf_LSF_1_LE_0_Files > 0
    plot(time_hws,work_LSF_1_LE_0,'LineWidth',2)
end
title('LSF: 1, LE: 0')
% axis([0, time_hws(end), 0, 0.06])
ylabel('Work (Joules)')

subplot(3,3,5)
if NumOf_LSF_1_LE_p2_Files > 0
    plot(time_hws,work_LSF_1_LE_p2,'LineWidth',2)
end
title('LSF: 1, LE: 0.2')
% axis([0, time_hws(end), 0, 0.06])
% ylabel('Work (Joules)')

subplot(3,3,6)
if NumOf_LSF_1_LE_p4_Files > 0
    plot(time_hws,work_LSF_1_LE_p4,'LineWidth',2)
end
title('LSF: 1, LE: 0.4')
% axis([0, time_hws(end), 0, 0.06])
% ylabel('Work (Joules)')

subplot(3,3,7)
if NumOf_LSF_2_LE_0_Files > 0
    plot(time_hws,work_LSF_2_LE_0,'LineWidth',2)
end
title('LSF: 2, LE: 0')
% axis([0, time_hws(end), 0, 0.06])
ylabel('Work (Joules)')
xlabel('Time (seconds)')

subplot(3,3,8)
if NumOf_LSF_2_LE_p2_Files > 0
    plot(time_hws,work_LSF_2_LE_p2,'LineWidth',2)
end
title('LSF: 1, LE: 0.2')
% axis([0, time_hws(end), 0, 0.06])
% ylabel('Work (Joules)')
xlabel('Time (seconds)')

subplot(3,3,9)
if NumOf_LSF_2_LE_p4_Files > 0
    plot(time_hws,work_LSF_2_LE_p4,'LineWidth',2)
end
title('LSF: 1, LE: 0.4')
% axis([0, time_hws(end), 0, 0.06])
% ylabel('Work (Joules)')
xlabel('Time (seconds)')

set(fig1, 'Position', [120,120, 1000, 600])

saveas(gcf,'../Figures_MPC/compSize/Fig26_sishManip_workIndiv.fig');
saveas(gcf,'../Figures_MPC/compSize/Fig26_sishManip_workIndiv.jpg');

disp('Figure 26 done')

%% Plot the work individualized by work term

fig1a = figure(1);
subplot(9,3,1)
if NumOf_LSF_p5_LE_0_Files > 0
    plot(time_hws,wingWork_LSF_p5_LE_0,'LineWidth',2)
end
% title('LSF: 0.5, LE: 0')
title('tauWing*theta')
axis([0, time_hws(end), 0, 0.00013])
ylabel('Work (Joules)')

subplot(9,3,2)
if NumOf_LSF_p5_LE_0_Files > 0
    plot(time_hws,abdoWork_LSF_p5_LE_0,'LineWidth',2)
end
title({'LSF: 0.5, LE: 0';'tauAbdo*beta'})
% title('abdomen work')
axis([0, time_hws(end), 0, 0.00013])

subplot(9,3,3)
if NumOf_LSF_p5_LE_0_Files > 0
    plot(time_hws,appFWork_LSF_p5_LE_0,'LineWidth',2)
end
title('F*distanceTraveled')
axis([0, time_hws(end), 0, 0.00013])

subplot(9,3,4)
if NumOf_LSF_p5_LE_p2_Files > 0
    plot(time_hws,wingWork_LSF_p5_LE_p2,'LineWidth',2)
end
title('tauWing*theta')
axis([0, time_hws(end), 0, 0.00013])

subplot(9,3,5)
if NumOf_LSF_p5_LE_p2_Files > 0
    plot(time_hws,abdoWork_LSF_p5_LE_p2,'LineWidth',2)
end
title({'LSF: 0.5, LE: 0.2';'tauAbdo*beta'})
axis([0, time_hws(end), 0, 0.00013])

subplot(9,3,6)
if NumOf_LSF_p5_LE_p2_Files > 0
    plot(time_hws,appFWork_LSF_p5_LE_p2,'LineWidth',2)
end
title('F*distanceTraveled')
axis([0, time_hws(end), 0, 0.00013])

subplot(9,3,7)
if NumOf_LSF_p5_LE_p4_Files > 0
    plot(time_hws,wingWork_LSF_p5_LE_p4,'LineWidth',2)
end
title('tauWing*theta')
axis([0, time_hws(end), 0, 0.00013])

subplot(9,3,8)
if NumOf_LSF_p5_LE_p4_Files > 0
    plot(time_hws,abdoWork_LSF_p5_LE_p4,'LineWidth',2)
end
title({'LSF: 0.5, LE: 0.4';'tauAbdo*beta'})
axis([0, time_hws(end), 0, 0.00013])

subplot(9,3,9)
if NumOf_LSF_p5_LE_p4_Files > 0
    plot(time_hws,appFWork_LSF_p5_LE_p4,'LineWidth',2)
end
title('F*distanceTraveled')
axis([0, time_hws(end), 0, 0.00013])

subplot(9,3,10)
if NumOf_LSF_1_LE_0_Files > 0
    plot(time_hws,wingWork_LSF_1_LE_0,'LineWidth',2)
end
title('tauWing*theta')
axis([0, time_hws(end), 0, 0.0008])
ylabel('Work (Joules)')

subplot(9,3,11)
if NumOf_LSF_1_LE_0_Files > 0
    plot(time_hws,abdoWork_LSF_1_LE_0,'LineWidth',2)
end
title({'LSF: 1, LE: 0';'tauAbdo*beta'})
axis([0, time_hws(end), 0, 0.0008])

subplot(9,3,12)
if NumOf_LSF_1_LE_0_Files > 0
    plot(time_hws,appFWork_LSF_1_LE_0,'LineWidth',2)
end
title('F*distanceTraveled')
axis([0, time_hws(end), 0, 0.0008])

subplot(9,3,13)
if NumOf_LSF_1_LE_p2_Files > 0
    plot(time_hws,wingWork_LSF_1_LE_p2,'LineWidth',2)
end
title('tauWing*theta')
axis([0, time_hws(end), 0, 0.0008])
% ylabel('Work (Joules)')

subplot(9,3,14)
if NumOf_LSF_1_LE_p2_Files > 0
    plot(time_hws,abdoWork_LSF_1_LE_p2,'LineWidth',2)
end
title({'LSF: 1, LE: 0.2';'tauAbdo*beta'})
axis([0, time_hws(end), 0, 0.0008])

subplot(9,3,15)
if NumOf_LSF_1_LE_p2_Files > 0
    plot(time_hws,appFWork_LSF_1_LE_p2,'LineWidth',2)
end
title('F*distanceTraveled')
axis([0, time_hws(end), 0, 0.0008])

subplot(9,3,16)
if NumOf_LSF_1_LE_p4_Files > 0
    plot(time_hws,wingWork_LSF_1_LE_p4,'LineWidth',2)
end
title('tauWing*theta')
axis([0, time_hws(end), 0, 0.0008])
% ylabel('Work (Joules)')

subplot(9,3,17)
if NumOf_LSF_1_LE_p4_Files > 0
    plot(time_hws,abdoWork_LSF_1_LE_p4,'LineWidth',2)
end
title({'LSF: 1, LE: 0.4';'tauAbdo*beta'})
axis([0, time_hws(end), 0, 0.0008])

subplot(9,3,18)
if NumOf_LSF_1_LE_p4_Files > 0
    plot(time_hws,appFWork_LSF_1_LE_p4,'LineWidth',2)
end
title('F*distanceTraveled')
axis([0, time_hws(end), 0, 0.0008])

subplot(9,3,19)
if NumOf_LSF_2_LE_0_Files > 0
    plot(time_hws,wingWork_LSF_2_LE_0,'LineWidth',2)
end
title('tauWing*theta')
axis([0, time_hws(end), 0, 0.04])
ylabel('Work (Joules)')

subplot(9,3,20)
if NumOf_LSF_2_LE_0_Files > 0
    plot(time_hws,abdoWork_LSF_2_LE_0,'LineWidth',2)
end
title({'LSF: 2, LE: 0';'tauAbdo*beta'})
xlabel('Time (seconds)')
axis([0, time_hws(end), 0, 0.04])

subplot(9,3,21)
if NumOf_LSF_2_LE_0_Files > 0
    plot(time_hws,appFWork_LSF_2_LE_0,'LineWidth',2)
end
title('F*distanceTraveled')
axis([0, time_hws(end), 0, 0.04])

subplot(9,3,22)
if NumOf_LSF_2_LE_p2_Files > 0
    plot(time_hws,wingWork_LSF_2_LE_p2,'LineWidth',2)
end
title('tauWing*theta')
axis([0, time_hws(end), 0, 0.04])

subplot(9,3,23)
if NumOf_LSF_2_LE_p2_Files > 0
    plot(time_hws,abdoWork_LSF_2_LE_p2,'LineWidth',2)
end
title({'LSF: 2, LE: 0.2';'tauAbdo*beta'})
xlabel('Time (seconds)')
axis([0, time_hws(end), 0, 0.04])

subplot(9,3,24)
if NumOf_LSF_2_LE_p2_Files > 0
    plot(time_hws,appFWork_LSF_2_LE_p2,'LineWidth',2)
end
title('F*distanceTraveled')
axis([0, time_hws(end), 0, 0.04])

subplot(9,3,25)
if NumOf_LSF_2_LE_p4_Files > 0
    plot(time_hws,wingWork_LSF_2_LE_p4,'LineWidth',2)
end
title('tauWing*theta')
axis([0, time_hws(end), 0, 0.04])

subplot(9,3,26)
if NumOf_LSF_2_LE_p4_Files > 0
    plot(time_hws,abdoWork_LSF_2_LE_p4,'LineWidth',2)
end
title({'LSF: 2, LE: 0.4';'tauAbdo*beta'})
xlabel('Time (seconds)')
axis([0, time_hws(end), 0, 0.04])

subplot(9,3,27)
if NumOf_LSF_2_LE_p4_Files > 0
    plot(time_hws,appFWork_LSF_2_LE_p4,'LineWidth',2)
end
title('F*distanceTraveled')
axis([0, time_hws(end), 0, 0.04])

set(fig1a, 'Position', [120,120, 1000, 600])

saveas(gcf,'../Figures_MPC/compSize/Fig26a_sishManip_workIndivByTerm.fig');
saveas(gcf,'../Figures_MPC/compSize/Fig26a_sishManip_workIndivByTerm.jpg');

disp('Figure 26a done')

%% Plot the average of work w.r.t. time

fig2 = figure(2);
subplot(3,3,1)
if NumOf_LSF_p5_LE_0_Files > 0
    Fig2_1 = shadedErrorBar(time_hws,mean_work_LSF_p5_LE_0,std_work_LSF_p5_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
    Fig2_1.mainLine.LineWidth = 2;
end
axis([0, time_hws(end), work_min, work_max])
title('LSF: 0.5, LE: 0')
% axis([0, time_hws(end), 0, 0.06])
ylabel('Work (Joules)')

subplot(3,3,2)
if NumOf_LSF_p5_LE_p2_Files > 0
    Fig2_2 = shadedErrorBar(time_hws,mean_work_LSF_p5_LE_p2,std_work_LSF_p5_LE_p2,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
    Fig2_2.mainLine.LineWidth = 2;
end
axis([0, time_hws(end), work_min, work_max])
% axis([0, time_hws(end), 0, 0.00051])
title('LSF: 0.5, LE: 0.2')

subplot(3,3,3)
if NumOf_LSF_p5_LE_p4_Files > 0
    Fig2_3 = shadedErrorBar(time_hws,mean_work_LSF_p5_LE_p4,std_work_LSF_p5_LE_p4,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
    Fig2_3.mainLine.LineWidth = 2;
end
axis([0, time_hws(end), work_min, work_max])
% axis([0, time_hws(end), 0, 0.00051])
title('LSF: 0.5, LE: 0.4')

subplot(3,3,4)
if NumOf_LSF_1_LE_0_Files > 0
    Fig2_4 = shadedErrorBar(time_hws,mean_work_LSF_1_LE_0,std_work_LSF_1_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
    Fig2_4.mainLine.LineWidth = 2;
end
axis([0, time_hws(end), work_min, work_max])
% axis([0, time_hws(end), 0, 0.00051])
ylabel('Work (Joules)')
title('LSF: 1, LE: 0')

subplot(3,3,5)
if NumOf_LSF_1_LE_p2_Files > 0
    Fig2_5 = shadedErrorBar(time_hws,mean_work_LSF_1_LE_p2,std_work_LSF_1_LE_p2,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
    Fig2_5.mainLine.LineWidth = 2;
end
axis([0, time_hws(end), work_min, work_max])
% axis([0, time_hws(end), 0, 0.00051])
title('LSF: 1, LE: 0.2')

subplot(3,3,6)
if NumOf_LSF_1_LE_p4_Files > 0
    Fig2_6 = shadedErrorBar(time_hws,mean_work_LSF_1_LE_p4,std_work_LSF_1_LE_p4,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
    Fig2_6.mainLine.LineWidth = 2;
end
axis([0, time_hws(end), work_min, work_max])
% axis([0, time_hws(end), 0, 0.00051])
title('LSF: 1, LE: 0.4')

subplot(3,3,7)
if NumOf_LSF_2_LE_0_Files > 0
    Fig2_7 = shadedErrorBar(time_hws,mean_work_LSF_2_LE_0,std_work_LSF_2_LE_0,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
    Fig2_7.mainLine.LineWidth = 2;
end
axis([0, time_hws(end), work_min, work_max])
% axis([0, time_hws(end), 0, 0.00051])
ylabel('Work (Joules)')
title('LSF: 2, LE: 0')
xlabel('Time (seconds)')

subplot(3,3,8)
if NumOf_LSF_2_LE_p2_Files > 0
    Fig2_8 = shadedErrorBar(time_hws,mean_work_LSF_2_LE_p2,std_work_LSF_2_LE_p2,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
    Fig2_8.mainLine.LineWidth = 2;
end
axis([0, time_hws(end), work_min, work_max])
% axis([0, time_hws(end), 0, 0.00051])
title('LSF: 2, LE: 0.2')
xlabel('Time (seconds)')

subplot(3,3,9)
if NumOf_LSF_2_LE_p4_Files > 0
    Fig2_9 = shadedErrorBar(time_hws,mean_work_LSF_2_LE_p4,std_work_LSF_2_LE_p4,...
    'lineprops','-k','transparent',false,'patchSaturation',0.075);
    Fig2_9.mainLine.LineWidth = 2;
end
axis([0, time_hws(end), work_min, work_max])
% axis([0, time_hws(end), 0, 0.00051])
title('LSF: 2, LE: 0.4')
xlabel('Time (seconds)')

set(fig2, 'Position', [120,120, 1000, 600])

saveas(gcf,'../Figures_MPC/compSize/Fig27_sishManip_workAvg.fig');
saveas(gcf,'../Figures_MPC/compSize/Fig27_sishManip_workAvg.jpg');

disp('Figure 27 done')

%% Work per full run distribution plot

fig3 = figure(3);
%Work - LSF_p5_LE_0
subplot(1,9,1)
title('LSF: 0.5, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_0_Files),')'])
if NumOf_LSF_p5_LE_0_Files > 0
    distributionPlot(mean_wpfr_LSF_p5_LE_0,'addSpread',true)
    axis([0, 2, meanwork_min, meanwork_max])
    % set(gca,'ytick',0:0.5:1);
%     axis([0, 2, 0.00015, 0.000177])
end
ylabel('Mean of work for each full run (Joules)')

%Work - LSF_p5_LE_p2
subplot(1,9,2)
title('LSF: 0.5, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p2_Files),')'])
if NumOf_LSF_p5_LE_p2_Files > 0
    distributionPlot(mean_wpfr_LSF_p5_LE_p2,'addSpread',true)
    axis([0, 2, meanwork_min, meanwork_max])
%     axis([0, 2, 0.00015, 0.000177])
end

%Work - LSF_p5_LE_p4
subplot(1,9,3)
title('LSF: 0.5, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p4_Files),')'])
if NumOf_LSF_p5_LE_p4_Files > 0
    distributionPlot(mean_wpfr_LSF_p5_LE_p4,'addSpread',true)
    axis([0, 2, meanwork_min, meanwork_max])
%     axis([0, 2, 0.125, 0.2])
end

%Work - LSF_1_LE_0
subplot(1,9,4)
title('LSF: 1, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_0_Files),')'])
if NumOf_LSF_1_LE_0_Files > 0
    distributionPlot(mean_wpfr_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, meanwork_min, meanwork_max])
%     axis([0, 2, 0.125, 0.2])
end

%Work - LSF_1_LE_p2
subplot(1,9,5)
title('LSF: 1, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p2_Files),')'])
if NumOf_LSF_1_LE_p2_Files > 0
    distributionPlot(mean_wpfr_LSF_1_LE_p2,'addSpread',true)
    axis([0, 2, meanwork_min, meanwork_max])
%     axis([0, 2, 0.125, 0.2])
end

%Work - LSF_1_LE_p4
subplot(1,9,6)
title('LSF: 1, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p4_Files),')'])
if NumOf_LSF_1_LE_p4_Files > 0
    distributionPlot(mean_wpfr_LSF_1_LE_p4,'addSpread',true)
    axis([0, 2, meanwork_min, meanwork_max])
%     axis([0, 2, 0.125, 0.2])
end

%Work - LSF_2_LE_0
subplot(1,9,7)
title('LSF: 2, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_0_Files),')'])
if NumOf_LSF_2_LE_0_Files > 0
    distributionPlot(mean_wpfr_LSF_2_LE_0,'addSpread',true)
    axis([0, 2, meanwork_min, meanwork_max])
%     axis([0, 2, 0.125, 0.2])
end

%Work - LSF_2_LE_p2
subplot(1,9,8)
title('LSF: 2, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p2_Files),')'])
if NumOf_LSF_2_LE_p2_Files > 0
    distributionPlot(mean_wpfr_LSF_2_LE_p2,'addSpread',true)
    axis([0, 2, meanwork_min, meanwork_max])
%     axis([0, 2, 0.125, 0.2])
end

%Work - LSF_2_LE_p4
subplot(1,9,9)
title('LSF: 2, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p4_Files),')'])
if NumOf_LSF_2_LE_p4_Files > 0
    distributionPlot(mean_wpfr_LSF_2_LE_p4,'addSpread',true)
    axis([0, 2, meanwork_min, meanwork_max])
%     axis([0, 2, 0.125, 0.2])
end

set(fig3, 'Position', [120,120, 1000, 600])

saveas(gcf,'../Figures_MPC/compSize/Fig28_sishManip_workViolin.fig');
saveas(gcf,'../Figures_MPC/compSize/Fig28_sishManip_workViolin.jpg');

disp('Figure 28 done')

%% Just LSF_p5's work 
fig3a = figure(3);
%Work - LSF_p5_LE_0
subplot(1,3,1)
title('LSF: 0.5, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_0_Files),')'])
if NumOf_LSF_p5_LE_0_Files > 0
    distributionPlot(mean_wpfr_LSF_p5_LE_0,'addSpread',true)
%     axis([0, 2, meanwork_min, meanwork_max])
    % set(gca,'ytick',0:0.5:1);
    axis([0, 2, 4.8e-6, 10e-6])
end
ylabel('Mean of work for each full run (Joules)')

%Work - LSF_p5_LE_p2
subplot(1,3,2)
title('LSF: 0.5, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p2_Files),')'])
if NumOf_LSF_p5_LE_p2_Files > 0
    distributionPlot(mean_wpfr_LSF_p5_LE_p2,'addSpread',true)
%     axis([0, 2, meanwork_min, meanwork_max])
    axis([0, 2, 4.8e-6, 10e-6])
end

%Work - LSF_p5_LE_p4
subplot(1,3,3)
title('LSF: 0.5, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p4_Files),')'])
if NumOf_LSF_p5_LE_p4_Files > 0
    distributionPlot(mean_wpfr_LSF_p5_LE_p4,'addSpread',true)
%     axis([0, 2, meanwork_min, meanwork_max])
    axis([0, 2, 4.8e-6, 10e-6])
end

set(fig3a, 'Position', [120,120, 1000, 600])

saveas(gcf,'../Figures_MPC/compSize/Fig28a_sishManip_workViolin_LSF_p5.fig');
saveas(gcf,'../Figures_MPC/compSize/Fig28a_sishManip_workViolin_LSF_p5.jpg');

disp('Figure 28a done')

%% Just LSF_1's work 
fig3b = figure(3);

%Work - LSF_1_LE_0
subplot(1,3,1)
title('LSF: 1, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_0_Files),')'])
if NumOf_LSF_1_LE_0_Files > 0
    distributionPlot(mean_wpfr_LSF_1_LE_0,'addSpread',true)
%     axis([0, 2, meanwork_min, meanwork_max])
    axis([0, 2, 0.00012, 0.00018])
end
ylabel('Mean of work for each full run (Joules)')


%Work - LSF_1_LE_p2
subplot(1,3,2)
title('LSF: 1, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p2_Files),')'])
if NumOf_LSF_1_LE_p2_Files > 0
    distributionPlot(mean_wpfr_LSF_1_LE_p2,'addSpread',true)
%     axis([0, 2, meanwork_min, meanwork_max])
    axis([0, 2, 0.00012, 0.00018])
end

%Work - LSF_1_LE_p4
subplot(1,3,3)
title('LSF: 1, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p4_Files),')'])
if NumOf_LSF_1_LE_p4_Files > 0
    distributionPlot(mean_wpfr_LSF_1_LE_p4,'addSpread',true)
%     axis([0, 2, meanwork_min, meanwork_max])
    axis([0, 2, 0.00012, 0.00018])
end

set(fig3b, 'Position', [120,120, 1000, 600])

saveas(gcf,'../Figures_MPC/compSize/Fig28b_sishManip_workViolin_LSF_1.fig');
saveas(gcf,'../Figures_MPC/compSize/Fig28b_sishManip_workViolin_LSF_1.jpg');

disp('Figure 28b done')

%% Non-dimensionalized tracking error per full run distribution plot

fig4 = figure(4);
%non-dim tepfr - LSF_p5_LE_0
subplot(1,9,1)
title('LSF: 0.5, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_0_Files),')'])
if NumOf_LSF_p5_LE_0_Files > 0
    distributionPlot(mean_nd_tepfr_LSF_p5_LE_0,'addSpread',true)
    axis([0, 2, meantepfr_nd_min, meantepfr_nd_max])
%     axis([0, 2, 0.00015, 0.000177])
end
ylabel('Mean of non-dimensionalized tracking error for each full run (unitless)')

%non-dim tepfr - LSF_p5_LE_p2
subplot(1,9,2)
title('LSF: 0.5, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p2_Files),')'])
if NumOf_LSF_p5_LE_p2_Files > 0
    distributionPlot(mean_nd_tepfr_LSF_p5_LE_p2,'addSpread',true)
    axis([0, 2, meantepfr_nd_min, meantepfr_nd_max])
%     axis([0, 2, 0.00015, 0.000177])
end

%non-dim tepfr - LSF_p5_LE_p4
subplot(1,9,3)
title('LSF: 0.5, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p4_Files),')'])
if NumOf_LSF_p5_LE_p4_Files > 0
    distributionPlot(mean_nd_tepfr_LSF_p5_LE_p4,'addSpread',true)
    axis([0, 2, meantepfr_nd_min, meantepfr_nd_max])
%     axis([0, 2, 0.125, 0.2])
end

%non-dim tepfr - LSF_1_LE_0
subplot(1,9,4)
title('LSF: 1, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_0_Files),')'])
if NumOf_LSF_1_LE_0_Files > 0
    distributionPlot(mean_nd_tepfr_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, meantepfr_nd_min, meantepfr_nd_max])
%     axis([0, 2, 0.125, 0.2])
end

%non-dim tepfr - LSF_1_LE_p2
subplot(1,9,5)
title('LSF: 1, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p2_Files),')'])
if NumOf_LSF_1_LE_p2_Files > 0
    distributionPlot(mean_nd_tepfr_LSF_1_LE_p2,'addSpread',true)
    axis([0, 2, meantepfr_nd_min, meantepfr_nd_max])
%     axis([0, 2, 0.125, 0.2])
end

%non-dim tepfr - LSF_1_LE_p4
subplot(1,9,6)
title('LSF: 1, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p4_Files),')'])
if NumOf_LSF_1_LE_p4_Files > 0
    distributionPlot(mean_nd_tepfr_LSF_1_LE_p4,'addSpread',true)
    axis([0, 2, meantepfr_nd_min, meantepfr_nd_max])
%     axis([0, 2, 0.125, 0.2])
end

%non-dim tepfr - LSF_2_LE_0
subplot(1,9,7)
title('LSF: 2, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_0_Files),')'])
if NumOf_LSF_2_LE_0_Files > 0
    distributionPlot(mean_nd_tepfr_LSF_2_LE_0,'addSpread',true)
    axis([0, 2, meantepfr_nd_min, meantepfr_nd_max])
%     axis([0, 2, 0.125, 0.2])
end

%non-dim tepfr - LSF_2_LE_p2
subplot(1,9,8)
title('LSF: 2, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p2_Files),')'])
if NumOf_LSF_2_LE_p2_Files > 0
    distributionPlot(mean_nd_tepfr_LSF_2_LE_p2,'addSpread',true)
    axis([0, 2, meantepfr_nd_min, meantepfr_nd_max])
%     axis([0, 2, 0.125, 0.2])
end

%non-dim tepfr - LSF_2_LE_p4
subplot(1,9,9)
title('LSF: 2, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p4_Files),')'])
if NumOf_LSF_2_LE_p4_Files > 0
    distributionPlot(mean_nd_tepfr_LSF_2_LE_p4,'addSpread',true)
    axis([0, 2, meantepfr_nd_min, meantepfr_nd_max])
%     axis([0, 2, 0.125, 0.2])
end

set(fig4, 'Position', [120,120, 1000, 600])

saveas(gcf,'../Figures_MPC/compSize/Fig29a_sishManip_nonDimTEpfr.fig');
saveas(gcf,'../Figures_MPC/compSize/Fig29a_sishManip_nonDimTEpfr.jpg');

disp('Figure 29a done')

%% Non-dimensionalized work per full run distribution plot

fig5 = figure(5);
%non-dim work - LSF_p5_LE_0
subplot(1,9,1)
title('LSF: 0.5, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_0_Files),')'])
if NumOf_LSF_p5_LE_0_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_p5_LE_0,'addSpread',true)
    axis([0, 2, meanwork_nd_min, meanwork_nd_max])
    % set(gca,'ytick',0:0.5:1);
%     axis([0, 2, 0.00015, 0.000177])
end
ylabel('Mean of non-dimensionalized work for each full run (unitless)')

%non-dim work - LSF_p5_LE_p2
subplot(1,9,2)
title('LSF: 0.5, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p2_Files),')'])
if NumOf_LSF_p5_LE_p2_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_p5_LE_p2,'addSpread',true)
    axis([0, 2, meanwork_nd_min, meanwork_nd_max])
%     axis([0, 2, 0.00015, 0.000177])
end

%non-dim work - LSF_p5_LE_p4
subplot(1,9,3)
title('LSF: 0.5, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p4_Files),')'])
if NumOf_LSF_p5_LE_p4_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_p5_LE_p4,'addSpread',true)
    axis([0, 2, meanwork_nd_min, meanwork_nd_max])
%     axis([0, 2, 0.125, 0.2])
end

%non-dim work - LSF_1_LE_0
subplot(1,9,4)
title('LSF: 1, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_0_Files),')'])
if NumOf_LSF_1_LE_0_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, meanwork_nd_min, meanwork_nd_max])
%     axis([0, 2, 0.125, 0.2])
end

%non-dim work - LSF_1_LE_p2
subplot(1,9,5)
title('LSF: 1, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p2_Files),')'])
if NumOf_LSF_1_LE_p2_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_1_LE_p2,'addSpread',true)
    axis([0, 2, meanwork_nd_min, meanwork_nd_max])
%     axis([0, 2, 0.125, 0.2])
end

%non-dim work - LSF_1_LE_p4
subplot(1,9,6)
title('LSF: 1, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p4_Files),')'])
if NumOf_LSF_1_LE_p4_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_1_LE_p4,'addSpread',true)
    axis([0, 2, meanwork_nd_min, meanwork_nd_max])
%     axis([0, 2, 0.125, 0.2])
end

%non-dim work - LSF_2_LE_0
subplot(1,9,7)
title('LSF: 2, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_0_Files),')'])
if NumOf_LSF_2_LE_0_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_2_LE_0,'addSpread',true)
    axis([0, 2, meanwork_nd_min, meanwork_nd_max])
%     axis([0, 2, 0.125, 0.2])
end

%non-dim work - LSF_2_LE_p2
subplot(1,9,8)
title('LSF: 2, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p2_Files),')'])
if NumOf_LSF_2_LE_p2_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_2_LE_p2,'addSpread',true)
    axis([0, 2, meanwork_nd_min, meanwork_nd_max])
%     axis([0, 2, 0.125, 0.2])
end

%non-dim work - LSF_2_LE_p4
subplot(1,9,9)
title('LSF: 2, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p4_Files),')'])
if NumOf_LSF_2_LE_p4_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_2_LE_p4,'addSpread',true)
    axis([0, 2, meanwork_nd_min, meanwork_nd_max])
%     axis([0, 2, 0.125, 0.2])
end

set(fig5, 'Position', [120,120, 1000, 600])

saveas(gcf,'../Figures_MPC/compSize/Fig30_sishManip_nonDimWorkViolin.fig');
saveas(gcf,'../Figures_MPC/compSize/Fig30_sishManip_nonDimWorkViolin.jpg');

disp('Figure 30 done')

%% Non-dimensionalized work per full run distribution plot -- JUST LSF_p5

fig5a = figure(5);
%non-dim work - LSF_p5_LE_0
subplot(1,3,1)
title('LSF: 0.5, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_0_Files),')'])
if NumOf_LSF_p5_LE_0_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_p5_LE_0,'addSpread',true)
%     axis([0, 2, meanwork_nd_min, meanwork_nd_max])
    % set(gca,'ytick',0:0.5:1);
    axis([0, 2, 14, 20])
end
ylabel('Mean of non-dimensionalized work for each full run (unitless)')

%non-dim work - LSF_p5_LE_p2
subplot(1,3,2)
title('LSF: 0.5, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p2_Files),')'])
if NumOf_LSF_p5_LE_p2_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_p5_LE_p2,'addSpread',true)
%     axis([0, 2, meanwork_nd_min, meanwork_nd_max])
    axis([0, 2, 14, 20])
end

%non-dim work - LSF_p5_LE_p4
subplot(1,3,3)
title('LSF: 0.5, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p4_Files),')'])
if NumOf_LSF_p5_LE_p4_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_p5_LE_p4,'addSpread',true)
%     axis([0, 2, meanwork_nd_min, meanwork_nd_max])
    axis([0, 2, 14, 20])
end

set(fig5a, 'Position', [120,120, 1000, 600])

saveas(gcf,'../Figures_MPC/compSize/Fig30a_sishManip_nonDimWorkViolin_LSF_p5.fig');
saveas(gcf,'../Figures_MPC/compSize/Fig30a_sishManip_nonDimWorkViolin_LSF_p5.jpg');

disp('Figure 30a done')

%% Non-dimensionalized work per full run distribution plot -- JUST LSF_1

fig5b = figure(5);
%non-dim work - LSF_1_LE_0
subplot(1,3,1)
title('LSF: 1, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_0_Files),')'])
if NumOf_LSF_1_LE_0_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_1_LE_0,'addSpread',true)
%     axis([0, 2, meanwork_nd_min, meanwork_nd_max])
    axis([0, 2, 60, 105])
end
ylabel('Mean of non-dimensionalized work for each full run (unitless)')

%non-dim work - LSF_1_LE_p2
subplot(1,3,2)
title('LSF: 1, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p2_Files),')'])
if NumOf_LSF_1_LE_p2_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_1_LE_p2,'addSpread',true)
%     axis([0, 2, meanwork_nd_min, meanwork_nd_max])
    axis([0, 2, 60, 105])
end

%non-dim work - LSF_1_LE_p4
subplot(1,3,3)
title('LSF: 1, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p4_Files),')'])
if NumOf_LSF_1_LE_p4_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_1_LE_p4,'addSpread',true)
%     axis([0, 2, meanwork_nd_min, meanwork_nd_max])
    axis([0, 2, 60, 105])
end

set(fig5b, 'Position', [120,120, 1000, 600])

saveas(gcf,'../Figures_MPC/compSize/Fig30b_sishManip_nonDimWorkViolin_LSF_1.fig');
saveas(gcf,'../Figures_MPC/compSize/Fig30b_sishManip_nonDimWorkViolin_LSF_1.jpg');

disp('Figure 30b done')

%% Non-dimensionalized product of work&tracking error per full run distribution plot

fig6 = figure(6);
%non-dim product - LSF_p5_LE_0
subplot(1,9,1)
title('LSF: 0.5, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_0_Files),')'])
if NumOf_LSF_p5_LE_0_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_p5_LE_0.*mean_nd_tepfr_LSF_p5_LE_0,'addSpread',true)
    axis([0, 2, meanProduct_min, meanProduct_max])
    % set(gca,'ytick',0:0.5:1);
%     axis([0, 2, 0.00015, 0.000177])
end
ylabel('Mean of non-dimensionalized work*tracking error for each full run (unitless)')

%non-dim product - LSF_p5_LE_p2
subplot(1,9,2)
title('LSF: 0.5, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p2_Files),')'])
if NumOf_LSF_p5_LE_p2_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_p5_LE_p2.*mean_nd_tepfr_LSF_p5_LE_p2,'addSpread',true)
    axis([0, 2, meanProduct_min, meanProduct_max])
%     axis([0, 2, 0.00015, 0.000177])
end

%non-dim product - LSF_p5_LE_p4
subplot(1,9,3)
title('LSF: 0.5, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p4_Files),')'])
if NumOf_LSF_p5_LE_p4_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_p5_LE_p4.*mean_nd_tepfr_LSF_p5_LE_p4,'addSpread',true)
    axis([0, 2, meanProduct_min, meanProduct_max])
%     axis([0, 2, 0.125, 0.2])
end

%non-dim product - LSF_1_LE_0
subplot(1,9,4)
title('LSF: 1, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_0_Files),')'])
if NumOf_LSF_1_LE_0_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_1_LE_0.*mean_nd_tepfr_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, meanProduct_min, meanProduct_max])
%     axis([0, 2, 0.125, 0.2])
end

%non-dim product - LSF_1_LE_p2
subplot(1,9,5)
title('LSF: 1, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p2_Files),')'])
if NumOf_LSF_1_LE_p2_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_1_LE_p2.*mean_nd_tepfr_LSF_1_LE_p2,'addSpread',true)
    axis([0, 2, meanProduct_min, meanProduct_max])
%     axis([0, 2, 0.125, 0.2])
end

%non-dim product - LSF_1_LE_p4
subplot(1,9,6)
title('LSF: 1, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p4_Files),')'])
if NumOf_LSF_1_LE_p4_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_1_LE_p4.*mean_nd_tepfr_LSF_1_LE_p4,'addSpread',true)
    axis([0, 2, meanProduct_min, meanProduct_max])
%     axis([0, 2, 0.125, 0.2])
end

%non-dim product - LSF_2_LE_0
subplot(1,9,7)
title('LSF: 2, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_0_Files),')'])
if NumOf_LSF_2_LE_0_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_2_LE_0.*mean_nd_tepfr_LSF_2_LE_0,'addSpread',true)
    axis([0, 2, meanProduct_min, meanProduct_max])
%     axis([0, 2, 0.125, 0.2])
end

%non-dim product - LSF_2_LE_p2
subplot(1,9,8)
title('LSF: 2, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p2_Files),')'])
if NumOf_LSF_2_LE_p2_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_2_LE_p2.*mean_nd_tepfr_LSF_2_LE_p2,'addSpread',true)
    axis([0, 2, meanProduct_min, meanProduct_max])
%     axis([0, 2, 0.125, 0.2])
end

%non-dim product - LSF_2_LE_p4
subplot(1,9,9)
title('LSF: 2, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p4_Files),')'])
if NumOf_LSF_2_LE_p4_Files > 0
    distributionPlot(mean_nd_wpfr_LSF_2_LE_p4.*mean_nd_tepfr_LSF_2_LE_p4,'addSpread',true)
    axis([0, 2, meanProduct_min, meanProduct_max])
%     axis([0, 2, 0.125, 0.2])
end

set(fig6, 'Position', [120,120, 1000, 600])

saveas(gcf,'../Figures_MPC/compSize/Fig31_sishManip_nonDimProductViolin.fig');
saveas(gcf,'../Figures_MPC/compSize/Fig31_sishManip_nonDimProductViolin.jpg');

disp('Figure 31 done')

%% Check the distance travelled for each inertial modificiation

figure(7);

subplot(3,3,1)
plot(time_hws,r_LSF_p5_LE_0)
axis([0, 10, 0, 1])
title('LSF: 0.5, LE: 0')

subplot(3,3,2)
plot(time_hws,r_LSF_p5_LE_p2)
axis([0, 10, 0, 1])
title({'Vertically oscillating';'LSF: 0.5, LE: 0.2'})

subplot(3,3,3)
plot(time_hws,r_LSF_p5_LE_p4)
axis([0, 10, 0, 1])
title('LSF: 0.5, LE: 0.4')

subplot(3,3,4)
plot(time_hws,r_LSF_1_LE_0)
axis([0, 10, 0, 1])
title('LSF: 1, LE: 0')
ylabel('Distance traveled (cm)')

subplot(3,3,5)
plot(time_hws,r_LSF_1_LE_p2)
axis([0, 10, 0, 1])
title('LSF: 1, LE: 0.2')

subplot(3,3,6)
plot(time_hws,r_LSF_1_LE_p4)
axis([0, 10, 0, 1])
title('LSF: 1, LE: 0.4')

subplot(3,3,7)
plot(time_hws,r_LSF_2_LE_0)
axis([0, 10, 0, 1])
title('LSF: 2, LE: 0')
xlabel('Time (seconds)')

subplot(3,3,8)
plot(time_hws,r_LSF_2_LE_p2)
axis([0, 10, 0, 1])
title('LSF: 2, LE: 0.2')
xlabel('Time (seconds)')

subplot(3,3,9)
plot(time_hws,r_LSF_2_LE_p4)
axis([0, 10, 0, 1])
title('LSF: 2, LE: 0.4')
xlabel('Time (seconds)')

%% Check the x-distance travelled for each inertial modificiation

figure(8);

subplot(3,3,1)
plot(time_hws,x_LSF_p5_LE_0)
axis([0, 10, -1, 1])
title('LSF: 0.5, LE: 0')

subplot(3,3,2)
plot(time_hws,x_LSF_p5_LE_p2)
axis([0, 10, -1, 1])
title({'Vertically oscillating';'LSF: 0.5, LE: 0.2'})

subplot(3,3,3)
plot(time_hws,x_LSF_p5_LE_p4)
axis([0, 10, -1, 1])
title('LSF: 0.5, LE: 0.4')

subplot(3,3,4)
plot(time_hws,x_LSF_1_LE_0)
axis([0, 10, -0.25, 1])
title('LSF: 1, LE: 0')
ylabel('x-distance traveled (cm)')

subplot(3,3,5)
plot(time_hws,x_LSF_1_LE_p2)
axis([0, 10, -0.25, 1])
title('LSF: 1, LE: 0.2')

subplot(3,3,6)
plot(time_hws,x_LSF_1_LE_p4)
axis([0, 10, -0.25, 1])
title('LSF: 1, LE: 0.4')

subplot(3,3,7)
plot(time_hws,x_LSF_2_LE_0)
axis([0, 10, -0.25, 1])
title('LSF: 2, LE: 0')
xlabel('Time (seconds)')

subplot(3,3,8)
plot(time_hws,x_LSF_2_LE_p2)
axis([0, 10, -0.25, 1])
title('LSF: 2, LE: 0.2')
xlabel('Time (seconds)')

subplot(3,3,9)
plot(time_hws,x_LSF_2_LE_p4)
axis([0, 10, -0.25, 1])
title('LSF: 2, LE: 0.4')
xlabel('Time (seconds)')

%% Check the y-distance travelled for each inertial modificiation

figure(9);

subplot(3,3,1)
plot(time_hws,y_LSF_p5_LE_0)
axis([0, 10, -0.25, 1])
title('LSF: 0.5, LE: 0')

subplot(3,3,2)
plot(time_hws,y_LSF_p5_LE_p2)
axis([0, 10, -0.25, 1])
title({'Vertically oscillating';'LSF: 0.5, LE: 0.2'})

subplot(3,3,3)
plot(time_hws,y_LSF_p5_LE_p4)
axis([0, 10, -0.25, 1])
title('LSF: 0.5, LE: 0.4')

subplot(3,3,4)
plot(time_hws,y_LSF_1_LE_0)
axis([0, 10, -0.25, 1])
title('LSF: 1, LE: 0')
ylabel('y-distance traveled (cm)')

subplot(3,3,5)
plot(time_hws,y_LSF_1_LE_p2)
axis([0, 10, -0.25, 1])
title('LSF: 1, LE: 0.2')

subplot(3,3,6)
plot(time_hws,y_LSF_1_LE_p4)
axis([0, 10, -0.25, 1])
title('LSF: 1, LE: 0.4')

subplot(3,3,7)
plot(time_hws,y_LSF_2_LE_0)
axis([0, 10, -0.25, 1])
title('LSF: 2, LE: 0')
xlabel('Time (seconds)')

subplot(3,3,8)
plot(time_hws,y_LSF_2_LE_p2)
axis([0, 10, -0.25, 1])
title('LSF: 2, LE: 0.2')
xlabel('Time (seconds)')

subplot(3,3,9)
plot(time_hws,y_LSF_2_LE_p4)
axis([0, 10, -0.25, 1])
title('LSF: 2, LE: 0.4')
xlabel('Time (seconds)')

%% Work per full run by term distribution plot -- JUST LSF_p5

fig10 = figure(10);

%wing work -  LSF_p5_LE_0
subplot(1,9,1)
title('tauWing*theta')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_0_Files),')'])
if NumOf_LSF_p5_LE_0_Files > 0
    distributionPlot(mean_wingWork_pfr_LSF_p5_LE_0./mean_wpfr_LSF_p5_LE_0,'addSpread',true)
    axis([0, 2, 0, 1])
    % set(gca,'ytick',0:0.5:1);
%     axis([0, 2, 0.00015, 0.000177])
end
ylabel('Fraction of work for each work term (unitless)')

%abdo work -  LSF_p5_LE_0
subplot(1,9,2)
title({'LSF: 0.5, LE: 0';'tauAbdo*beta'})
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_0_Files),')'])
if NumOf_LSF_p5_LE_0_Files > 0
    distributionPlot(mean_abdoWork_pfr_LSF_p5_LE_0./mean_wpfr_LSF_p5_LE_0,'addSpread',true)
    axis([0, 2, 0, 1])
%     axis([0, 2, 0.00015, 0.000177])
end

%F*distance work -  LSF_p5_LE_0
subplot(1,9,3)
title('F*dist')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_0_Files),')'])
if NumOf_LSF_p5_LE_0_Files > 0
    distributionPlot(mean_appFWork_pfr_LSF_p5_LE_0./mean_wpfr_LSF_p5_LE_0,'addSpread',true)
    axis([0, 2, 0, 1])
%     axis([0, 2, 0.00015, 0.000177])
end

%wing work -  LSF_p5_LE_p2
subplot(1,9,4)
title('tauWing*theta')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p2_Files),')'])
if NumOf_LSF_p5_LE_p2_Files > 0
    distributionPlot(mean_wingWork_pfr_LSF_p5_LE_p2./mean_wpfr_LSF_p5_LE_p2,'addSpread',true)
    axis([0, 2, 0, 1])
    % set(gca,'ytick',0:0.5:1);
%     axis([0, 2, 0.00015, 0.000177])
end

%abdo work -  LSF_p5_LE_p2
subplot(1,9,5)
title({'LSF: 0.5, LE: 0.2';'tauAbdo*beta'})
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p2_Files),')'])
if NumOf_LSF_p5_LE_p2_Files > 0
    distributionPlot(mean_abdoWork_pfr_LSF_p5_LE_p2./mean_wpfr_LSF_p5_LE_p2,'addSpread',true)
    axis([0, 2, 0, 1])
%     axis([0, 2, 0.00015, 0.000177])
end

%F*distance work -  LSF_p5_LE_p2
subplot(1,9,6)
title('F*dist')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p2_Files),')'])
if NumOf_LSF_p5_LE_p2_Files > 0
    distributionPlot(mean_appFWork_pfr_LSF_p5_LE_p2./mean_wpfr_LSF_p5_LE_p2,'addSpread',true)
    axis([0, 2, 0, 1])
%     axis([0, 2, 0.00015, 0.000177])
end

%wing work -  LSF_p5_LE_p4
subplot(1,9,7)
title('tauWing*theta')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p4_Files),')'])
if NumOf_LSF_p5_LE_p4_Files > 0
    distributionPlot(mean_wingWork_pfr_LSF_p5_LE_p4./mean_wpfr_LSF_p5_LE_p4,'addSpread',true)
    axis([0, 2, 0, 1])
    % set(gca,'ytick',0:0.5:1);
%     axis([0, 2, 0.00015, 0.000177])
end

%abdo work -  LSF_p5_LE_p4
subplot(1,9,8)
title({'LSF: 0.5, LE: 0.4';'tauAbdo*beta'})
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p4_Files),')'])
if NumOf_LSF_p5_LE_p4_Files > 0
    distributionPlot(mean_abdoWork_pfr_LSF_p5_LE_p4./mean_wpfr_LSF_p5_LE_p4,'addSpread',true)
    axis([0, 2, 0, 1])
%     axis([0, 2, 0.00015, 0.000177])
end

%F*distance work -  LSF_p5_LE_p4
subplot(1,9,9)
title('F*dist')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p4_Files),')'])
if NumOf_LSF_p5_LE_p4_Files > 0
    distributionPlot(mean_appFWork_pfr_LSF_p5_LE_p4./mean_wpfr_LSF_p5_LE_p4,'addSpread',true)
    axis([0, 2, 0, 1])
%     axis([0, 2, 0.00015, 0.000177])
end

set(fig10, 'Position', [120,120, 1000, 600])

saveas(gcf,'../Figures_MPC/compSize/Fig32a_sishManip_workTermsViolin_LSF_p5.fig');
saveas(gcf,'../Figures_MPC/compSize/Fig32a_sishManip_workTermsViolin_LSF_p5.jpg');

disp('Figure 32a done')

%% Work per full run by term distribution plot -- JUST LSF_1

fig11 = figure(11);

%wing work -  LSF_1_LE_0
subplot(1,9,1)
title('tauWing*theta')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_0_Files),')'])
if NumOf_LSF_1_LE_0_Files > 0
    distributionPlot(mean_wingWork_pfr_LSF_1_LE_0./mean_wpfr_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, 0, 1])
    % set(gca,'ytick',0:0.5:1);
%     axis([0, 2, 0.00015, 0.000177])
end
ylabel('Fraction of work for each work term (unitless)')

%abdo work -  LSF_1_LE_0
subplot(1,9,2)
title({'LSF: 1, LE: 0';'tauAbdo*beta'})
xlabel(['(n=',num2str(NumOf_LSF_1_LE_0_Files),')'])
if NumOf_LSF_1_LE_0_Files > 0
    distributionPlot(mean_abdoWork_pfr_LSF_1_LE_0./mean_wpfr_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, 0, 1])
%     axis([0, 2, 0.00015, 0.000177])
end

%F*distance work -  LSF_1_LE_0
subplot(1,9,3)
title('F*dist')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_0_Files),')'])
if NumOf_LSF_1_LE_0_Files > 0
    distributionPlot(mean_appFWork_pfr_LSF_1_LE_0./mean_wpfr_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, 0, 1])
%     axis([0, 2, 0.00015, 0.000177])
end

%wing work -  LSF_1_LE_p2
subplot(1,9,4)
title('tauWing*theta')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p2_Files),')'])
if NumOf_LSF_1_LE_p2_Files > 0
    distributionPlot(mean_wingWork_pfr_LSF_1_LE_p2./mean_wpfr_LSF_1_LE_p2,'addSpread',true)
    axis([0, 2, 0, 1])
    % set(gca,'ytick',0:0.5:1);
%     axis([0, 2, 0.00015, 0.000177])
end

%abdo work -  LSF_1_LE_p2
subplot(1,9,5)
title({'LSF: 1, LE: 0.2';'tauAbdo*beta'})
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p2_Files),')'])
if NumOf_LSF_1_LE_p2_Files > 0
    distributionPlot(mean_abdoWork_pfr_LSF_1_LE_p2./mean_wpfr_LSF_1_LE_p2,'addSpread',true)
    axis([0, 2, 0, 1])
%     axis([0, 2, 0.00015, 0.000177])
end

%F*distance work -  LSF_1_LE_p2
subplot(1,9,6)
title('F*dist')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p2_Files),')'])
if NumOf_LSF_1_LE_p2_Files > 0
    distributionPlot(mean_appFWork_pfr_LSF_1_LE_p2./mean_wpfr_LSF_1_LE_p2,'addSpread',true)
    axis([0, 2, 0, 1])
%     axis([0, 2, 0.00015, 0.000177])
end

%wing work -  LSF_1_LE_p4
subplot(1,9,7)
title('tauWing*theta')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p4_Files),')'])
if NumOf_LSF_1_LE_p4_Files > 0
    distributionPlot(mean_wingWork_pfr_LSF_1_LE_p4./mean_wpfr_LSF_1_LE_p4,'addSpread',true)
    axis([0, 2, 0, 1])
    % set(gca,'ytick',0:0.5:1);
%     axis([0, 2, 0.00015, 0.000177])
end

%abdo work -  LSF_1_LE_p4
subplot(1,9,8)
title({'LSF: 1, LE: 0.4';'tauAbdo*beta'})
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p4_Files),')'])
if NumOf_LSF_1_LE_p4_Files > 0
    distributionPlot(mean_abdoWork_pfr_LSF_1_LE_p4./mean_wpfr_LSF_1_LE_p4,'addSpread',true)
    axis([0, 2, 0, 1])
%     axis([0, 2, 0.00015, 0.000177])
end

%F*distance work -  LSF_1_LE_p4
subplot(1,9,9)
title('F*dist')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p4_Files),')'])
if NumOf_LSF_1_LE_p4_Files > 0
    distributionPlot(mean_appFWork_pfr_LSF_1_LE_p4./mean_wpfr_LSF_1_LE_p4,'addSpread',true)
    axis([0, 2, 0, 1])
%     axis([0, 2, 0.00015, 0.000177])
end

set(fig11, 'Position', [120,120, 1000, 600])

saveas(gcf,'../Figures_MPC/compSize/Fig32b_sishManip_workTermsViolin_LSF_1.fig');
saveas(gcf,'../Figures_MPC/compSize/Fig32b_sishManip_workTermsViolin_LSF_1.jpg');

disp('Figure 32b done')

%% Work per full run by term distribution plot -- JUST LSF_2

fig12 = figure(12);

%wing work -  LSF_2_LE_0
subplot(1,9,1)
title('tauWing*theta')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_0_Files),')'])
if NumOf_LSF_2_LE_0_Files > 0
    distributionPlot(mean_wingWork_pfr_LSF_2_LE_0./mean_wpfr_LSF_2_LE_0,'addSpread',true)
    axis([0, 2, 0, 1])
    % set(gca,'ytick',0:0.5:1);
end
ylabel('Fraction of work for each work term (unitless)')

%abdo work -  LSF_2_LE_0
subplot(1,9,2)
title({'LSF: 2, LE: 0';'tauAbdo*beta'})
xlabel(['(n=',num2str(NumOf_LSF_2_LE_0_Files),')'])
if NumOf_LSF_2_LE_0_Files > 0
    distributionPlot(mean_abdoWork_pfr_LSF_2_LE_0./mean_wpfr_LSF_2_LE_0,'addSpread',true)
    axis([0, 2, 0, 1])
end

%F*distance work -  LSF_2_LE_0
subplot(1,9,3)
title('F*dist')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_0_Files),')'])
if NumOf_LSF_2_LE_0_Files > 0
    distributionPlot(mean_appFWork_pfr_LSF_2_LE_0./mean_wpfr_LSF_2_LE_0,'addSpread',true)
    axis([0, 2, 0, 1])
end

%wing work -  LSF_2_LE_p2
subplot(1,9,4)
title('tauWing*theta')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p2_Files),')'])
if NumOf_LSF_2_LE_p2_Files > 0
    distributionPlot(mean_wingWork_pfr_LSF_2_LE_p2./mean_wpfr_LSF_2_LE_p2,'addSpread',true)
    axis([0, 2, 0, 1])
    % set(gca,'ytick',0:0.5:1);
end

%abdo work -  LSF_2_LE_p2
subplot(1,9,5)
title({'LSF: 2, LE: 0.2';'tauAbdo*beta'})
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p2_Files),')'])
if NumOf_LSF_2_LE_p2_Files > 0
    distributionPlot(mean_abdoWork_pfr_LSF_2_LE_p2./mean_wpfr_LSF_2_LE_p2,'addSpread',true)
    axis([0, 2, 0, 1])
end

%F*distance work -  LSF_2_LE_p2
subplot(1,9,6)
title('F*dist')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p2_Files),')'])
if NumOf_LSF_2_LE_p2_Files > 0
    distributionPlot(mean_appFWork_pfr_LSF_2_LE_p2./mean_wpfr_LSF_2_LE_p2,'addSpread',true)
    axis([0, 2, 0, 1])
end

%wing work -  LSF_2_LE_p4
subplot(1,9,7)
title('tauWing*theta')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p4_Files),')'])
if NumOf_LSF_2_LE_p4_Files > 0
    distributionPlot(mean_wingWork_pfr_LSF_2_LE_p4./mean_wpfr_LSF_2_LE_p4,'addSpread',true)
    axis([0, 2, 0, 1])
    % set(gca,'ytick',0:0.5:1);
end

%abdo work -  LSF_2_LE_p4
subplot(1,9,8)
title({'LSF: 2, LE: 0.4';'tauAbdo*beta'})
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p4_Files),')'])
if NumOf_LSF_2_LE_p4_Files > 0
    distributionPlot(mean_abdoWork_pfr_LSF_2_LE_p4./mean_wpfr_LSF_2_LE_p4,'addSpread',true)
    axis([0, 2, 0, 1])
end

%F*distance work -  LSF_2_LE_p4
subplot(1,9,9)
title('F*dist')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p4_Files),')'])
if NumOf_LSF_2_LE_p4_Files > 0
    distributionPlot(mean_appFWork_pfr_LSF_2_LE_p4./mean_wpfr_LSF_2_LE_p4,'addSpread',true)
    axis([0, 2, 0, 1])
end

set(fig12, 'Position', [120,120, 1000, 600])

saveas(gcf,'../Figures_MPC/compSize/Fig32c_sishManip_workTermsViolin_LSF_2.fig');
saveas(gcf,'../Figures_MPC/compSize/Fig32c_sishManip_workTermsViolin_LSF_2.jpg');

disp('Figure 32c done')

%% Non-dimensional work per full run by term distribution plot -- JUST LSF_p5

fig13 = figure(13);

%wing work_nd -  LSF_p5_LE_0
subplot(1,9,1)
title('tauWing*theta')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_0_Files),')'])
if NumOf_LSF_p5_LE_0_Files > 0
    distributionPlot(mean_nd_wingWork_pfr_LSF_p5_LE_0./mean_nd_wpfr_LSF_p5_LE_0,'addSpread',true)
    axis([0, 2, 0, 1])
    % set(gca,'ytick',0:0.5:1);
end
ylabel('Fraction of non-dimensional work for each work term (unitless)')

%abdo work_nd -  LSF_p5_LE_0
subplot(1,9,2)
title({'LSF: 0.5, LE: 0';'tauAbdo*beta'})
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_0_Files),')'])
if NumOf_LSF_p5_LE_0_Files > 0
    distributionPlot(mean_nd_abdoWork_pfr_LSF_p5_LE_0./mean_nd_wpfr_LSF_p5_LE_0,'addSpread',true)
    axis([0, 2, 0, 1])
end

%F*distance work_nd -  LSF_p5_LE_0
subplot(1,9,3)
title('F*dist')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_0_Files),')'])
if NumOf_LSF_p5_LE_0_Files > 0
    distributionPlot(mean_nd_appFWork_pfr_LSF_p5_LE_0./mean_nd_wpfr_LSF_p5_LE_0,'addSpread',true)
    axis([0, 2, 0, 1])
end

%wing work_nd -  LSF_p5_LE_p2
subplot(1,9,4)
title('tauWing*theta')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p2_Files),')'])
if NumOf_LSF_p5_LE_p2_Files > 0
    distributionPlot(mean_nd_wingWork_pfr_LSF_p5_LE_p2./mean_nd_wpfr_LSF_p5_LE_p2,'addSpread',true)
    axis([0, 2, 0, 1])
    % set(gca,'ytick',0:0.5:1);
end

%abdo work_nd -  LSF_p5_LE_p2
subplot(1,9,5)
title({'LSF: 0.5, LE: 0.2';'tauAbdo*beta'})
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p2_Files),')'])
if NumOf_LSF_p5_LE_p2_Files > 0
    distributionPlot(mean_nd_abdoWork_pfr_LSF_p5_LE_p2./mean_nd_wpfr_LSF_p5_LE_p2,'addSpread',true)
    axis([0, 2, 0, 1])
end

%F*distance work_nd -  LSF_p5_LE_p2
subplot(1,9,6)
title('F*dist')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p2_Files),')'])
if NumOf_LSF_p5_LE_p2_Files > 0
    distributionPlot(mean_nd_appFWork_pfr_LSF_p5_LE_p2./mean_nd_wpfr_LSF_p5_LE_p2,'addSpread',true)
    axis([0, 2, 0, 1])
end

%wing work_nd -  LSF_p5_LE_p4
subplot(1,9,7)
title('tauWing*theta')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p4_Files),')'])
if NumOf_LSF_p5_LE_p4_Files > 0
    distributionPlot(mean_nd_wingWork_pfr_LSF_p5_LE_p4./mean_nd_wpfr_LSF_p5_LE_p4,'addSpread',true)
    axis([0, 2, 0, 1])
end

%abdo work_nd -  LSF_p5_LE_p4
subplot(1,9,8)
title({'LSF: 0.5, LE: 0.4';'tauAbdo*beta'})
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p4_Files),')'])
if NumOf_LSF_p5_LE_p4_Files > 0
    distributionPlot(mean_nd_abdoWork_pfr_LSF_p5_LE_p4./mean_nd_wpfr_LSF_p5_LE_p4,'addSpread',true)
    axis([0, 2, 0, 1])
end

%F*distance work_nd -  LSF_p5_LE_p4
subplot(1,9,9)
title('F*dist')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p4_Files),')'])
if NumOf_LSF_p5_LE_p4_Files > 0
    distributionPlot(mean_nd_appFWork_pfr_LSF_p5_LE_p4./mean_nd_wpfr_LSF_p5_LE_p4,'addSpread',true)
    axis([0, 2, 0, 1])
end

set(fig13, 'Position', [120,120, 1000, 600])

saveas(gcf,'../Figures_MPC/compSize/Fig33a_sishManip_nd_workTermsViolin_LSF_p5.fig');
saveas(gcf,'../Figures_MPC/compSize/Fig33a_sishManip_nd_workTermsViolin_LSF_p5.jpg');

disp('Figure 33a done')

%% Non-dimensional work per full run by term distribution plot -- JUST LSF_1

fig14 = figure(14);

%wing work_nd -  LSF_1_LE_0
subplot(1,9,1)
title('tauWing*theta')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_0_Files),')'])
if NumOf_LSF_1_LE_0_Files > 0
    distributionPlot(mean_nd_wingWork_pfr_LSF_1_LE_0./mean_nd_wpfr_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, 0, 1])
    % set(gca,'ytick',0:0.5:1);
end
ylabel('Fraction of non-dimensional work for each work term (unitless)')

%abdo work_nd -  LSF_1_LE_0
subplot(1,9,2)
title({'LSF: 1, LE: 0';'tauAbdo*beta'})
xlabel(['(n=',num2str(NumOf_LSF_1_LE_0_Files),')'])
if NumOf_LSF_1_LE_0_Files > 0
    distributionPlot(mean_nd_abdoWork_pfr_LSF_1_LE_0./mean_nd_wpfr_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, 0, 1])
end

%F*distance work_nd -  LSF_1_LE_0
subplot(1,9,3)
title('F*dist')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_0_Files),')'])
if NumOf_LSF_1_LE_0_Files > 0
    distributionPlot(mean_nd_appFWork_pfr_LSF_1_LE_0./mean_nd_wpfr_LSF_1_LE_0,'addSpread',true)
    axis([0, 2, 0, 1])
end

%wing work_nd -  LSF_1_LE_p2
subplot(1,9,4)
title('tauWing*theta')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p2_Files),')'])
if NumOf_LSF_1_LE_p2_Files > 0
    distributionPlot(mean_nd_wingWork_pfr_LSF_1_LE_p2./mean_nd_wpfr_LSF_1_LE_p2,'addSpread',true)
    axis([0, 2, 0, 1])
    % set(gca,'ytick',0:0.5:1);
end

%abdo work_nd -  LSF_1_LE_p2
subplot(1,9,5)
title({'LSF: 1, LE: 0.2';'tauAbdo*beta'})
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p2_Files),')'])
if NumOf_LSF_1_LE_p2_Files > 0
    distributionPlot(mean_nd_abdoWork_pfr_LSF_1_LE_p2./mean_nd_wpfr_LSF_1_LE_p2,'addSpread',true)
    axis([0, 2, 0, 1])
end

%F*distance work_nd -  LSF_1_LE_p2
subplot(1,9,6)
title('F*dist')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p2_Files),')'])
if NumOf_LSF_1_LE_p2_Files > 0
    distributionPlot(mean_nd_appFWork_pfr_LSF_1_LE_p2./mean_nd_wpfr_LSF_1_LE_p2,'addSpread',true)
    axis([0, 2, 0, 1])
end

%wing work_nd -  LSF_1_LE_p4
subplot(1,9,7)
title('tauWing*theta')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p4_Files),')'])
if NumOf_LSF_1_LE_p4_Files > 0
    distributionPlot(mean_nd_wingWork_pfr_LSF_1_LE_p4./mean_nd_wpfr_LSF_1_LE_p4,'addSpread',true)
    axis([0, 2, 0, 1])
end

%abdo work_nd -  LSF_1_LE_p4
subplot(1,9,8)
title({'LSF: 1, LE: 0.4';'tauAbdo*beta'})
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p4_Files),')'])
if NumOf_LSF_1_LE_p4_Files > 0
    distributionPlot(mean_nd_abdoWork_pfr_LSF_1_LE_p4./mean_nd_wpfr_LSF_1_LE_p4,'addSpread',true)
    axis([0, 2, 0, 1])
end

%F*distance work_nd -  LSF_1_LE_p4
subplot(1,9,9)
title('F*dist')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p4_Files),')'])
if NumOf_LSF_1_LE_p4_Files > 0
    distributionPlot(mean_nd_appFWork_pfr_LSF_1_LE_p4./mean_nd_wpfr_LSF_1_LE_p4,'addSpread',true)
    axis([0, 2, 0, 1])
end

set(fig14, 'Position', [120,120, 1000, 600])

saveas(gcf,'../Figures_MPC/compSize/Fig33b_sishManip_nd_workTermsViolin_LSF_1.fig');
saveas(gcf,'../Figures_MPC/compSize/Fig33b_sishManip_nd_workTermsViolin_LSF_1.jpg');

disp('Figure 33b done')

%% Non-dimensional work per full run by term distribution plot -- JUST LSF_2

fig15 = figure(15);

%wing work_nd -  LSF_2_LE_0
subplot(1,9,1)
title('tauWing*theta')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_0_Files),')'])
if NumOf_LSF_2_LE_0_Files > 0
    distributionPlot(mean_nd_wingWork_pfr_LSF_2_LE_0./mean_nd_wpfr_LSF_2_LE_0,'addSpread',true)
    axis([0, 2, 0, 1])
    % set(gca,'ytick',0:0.5:1);
end
ylabel('Fraction of non-dimensional work for each work term (unitless)')

%abdo work_nd -  LSF_2_LE_0
subplot(1,9,2)
title({'LSF: 2, LE: 0';'tauAbdo*beta'})
xlabel(['(n=',num2str(NumOf_LSF_2_LE_0_Files),')'])
if NumOf_LSF_2_LE_0_Files > 0
    distributionPlot(mean_nd_abdoWork_pfr_LSF_2_LE_0./mean_nd_wpfr_LSF_2_LE_0,'addSpread',true)
    axis([0, 2, 0, 1])
end

%F*distance work_nd -  LSF_2_LE_0
subplot(1,9,3)
title('F*dist')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_0_Files),')'])
if NumOf_LSF_2_LE_0_Files > 0
    distributionPlot(mean_nd_appFWork_pfr_LSF_2_LE_0./mean_nd_wpfr_LSF_2_LE_0,'addSpread',true)
    axis([0, 2, 0, 1])
end

%wing work_nd -  LSF_2_LE_p2
subplot(1,9,4)
title('tauWing*theta')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p2_Files),')'])
if NumOf_LSF_2_LE_p2_Files > 0
    distributionPlot(mean_nd_wingWork_pfr_LSF_2_LE_p2./mean_nd_wpfr_LSF_2_LE_p2,'addSpread',true)
    axis([0, 2, 0, 1])
    % set(gca,'ytick',0:0.5:1);
end

%abdo work_nd -  LSF_2_LE_p2
subplot(1,9,5)
title({'LSF: 2, LE: 0.2';'tauAbdo*beta'})
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p2_Files),')'])
if NumOf_LSF_2_LE_p2_Files > 0
    distributionPlot(mean_nd_abdoWork_pfr_LSF_2_LE_p2./mean_nd_wpfr_LSF_2_LE_p2,'addSpread',true)
    axis([0, 2, 0, 1])
end

%F*distance work_nd -  LSF_2_LE_p2
subplot(1,9,6)
title('F*dist')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p2_Files),')'])
if NumOf_LSF_2_LE_p2_Files > 0
    distributionPlot(mean_nd_appFWork_pfr_LSF_2_LE_p2./mean_nd_wpfr_LSF_2_LE_p2,'addSpread',true)
    axis([0, 2, 0, 1])
end

%wing work_nd -  LSF_2_LE_p4
subplot(1,9,7)
title('tauWing*theta')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p4_Files),')'])
if NumOf_LSF_2_LE_p4_Files > 0
    distributionPlot(mean_nd_wingWork_pfr_LSF_2_LE_p4./mean_nd_wpfr_LSF_2_LE_p4,'addSpread',true)
    axis([0, 2, 0, 1])
end

%abdo work_nd -  LSF_2_LE_p4
subplot(1,9,8)
title({'LSF: 2, LE: 0.4';'tauAbdo*beta'})
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p4_Files),')'])
if NumOf_LSF_2_LE_p4_Files > 0
    distributionPlot(mean_nd_abdoWork_pfr_LSF_2_LE_p4./mean_nd_wpfr_LSF_2_LE_p4,'addSpread',true)
    axis([0, 2, 0, 1])
end

%F*distance work_nd -  LSF_2_LE_p4
subplot(1,9,9)
title('F*dist')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p4_Files),')'])
if NumOf_LSF_2_LE_p4_Files > 0
    distributionPlot(mean_nd_appFWork_pfr_LSF_2_LE_p4./mean_nd_wpfr_LSF_2_LE_p4,'addSpread',true)
    axis([0, 2, 0, 1])
end

set(fig15, 'Position', [120,120, 1000, 600])

saveas(gcf,'../Figures_MPC/compSize/Fig33c_sishManip_nd_workTermsViolin_LSF_2.fig');
saveas(gcf,'../Figures_MPC/compSize/Fig33c_sishManip_nd_workTermsViolin_LSF_2.jpg');

disp('Figure 33c done')

%% Tracking error per full run distribution plot (in cm)

fig16 = figure(16);
%tepfr - LSF_p5_LE_0
subplot(1,9,1)
title('LSF: 0.5, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_0_Files),')'])
if NumOf_LSF_p5_LE_0_Files > 0
    distributionPlot(bl_LSF_p5_LE_0.*mean_nd_tepfr_LSF_p5_LE_0,'addSpread',true)
%     axis([0, 2, meantepfr_nd_min, meantepfr_nd_max])
    axis([0, 2, 0.1, 0.32])
    yticks(0.1:0.05:0.3)
end
ylabel({'Mean of tracking error for each full run (cm)';...
    'Vertical flower tracking'})

%tepfr - LSF_p5_LE_p2
subplot(1,9,2)
title('LSF: 0.5, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p2_Files),')'])
if NumOf_LSF_p5_LE_p2_Files > 0
    distributionPlot(bl_LSF_p5_LE_p2.*mean_nd_tepfr_LSF_p5_LE_p2,'addSpread',true)
%     axis([0, 2, meantepfr_nd_min, meantepfr_nd_max])
    axis([0, 2, 0.1, 0.32])
    yticks(0.1:0.05:0.3)
end

%tepfr - LSF_p5_LE_p4
subplot(1,9,3)
title('LSF: 0.5, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_p5_LE_p4_Files),')'])
if NumOf_LSF_p5_LE_p4_Files > 0
    distributionPlot(bl_LSF_p5_LE_p4.*mean_nd_tepfr_LSF_p5_LE_p4,'addSpread',true)
%     axis([0, 2, meantepfr_nd_min, meantepfr_nd_max])
    axis([0, 2, 0.1, 0.32])
    yticks(0.1:0.05:0.3)
end

%tepfr - LSF_1_LE_0
subplot(1,9,4)
title('LSF: 1, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_0_Files),')'])
if NumOf_LSF_1_LE_0_Files > 0
    distributionPlot(bl_LSF_1_LE_0.*mean_nd_tepfr_LSF_1_LE_0,'addSpread',true)
%     axis([0, 2, meantepfr_nd_min, meantepfr_nd_max])
    axis([0, 2, 0.1, 0.32])
    yticks(0.1:0.05:0.3)
end

%tepfr - LSF_1_LE_p2
subplot(1,9,5)
title('LSF: 1, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p2_Files),')'])
if NumOf_LSF_1_LE_p2_Files > 0
    distributionPlot(bl_LSF_1_LE_p2.*mean_nd_tepfr_LSF_1_LE_p2,'addSpread',true)
%     axis([0, 2, meantepfr_nd_min, meantepfr_nd_max])
    axis([0, 2, 0.1, 0.32])
    yticks(0.1:0.05:0.3)
end

%tepfr - LSF_1_LE_p4
subplot(1,9,6)
title('LSF: 1, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_1_LE_p4_Files),')'])
if NumOf_LSF_1_LE_p4_Files > 0
    distributionPlot(bl_LSF_1_LE_p4.*mean_nd_tepfr_LSF_1_LE_p4,'addSpread',true)
%     axis([0, 2, meantepfr_nd_min, meantepfr_nd_max])
    axis([0, 2, 0.1, 0.32])
    yticks(0.1:0.05:0.3)
end

%tepfr - LSF_2_LE_0
subplot(1,9,7)
title('LSF: 2, LE: 0')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_0_Files),')'])
if NumOf_LSF_2_LE_0_Files > 0
    distributionPlot(bl_LSF_2_LE_0.*mean_nd_tepfr_LSF_2_LE_0,'addSpread',true)
%     axis([0, 2, meantepfr_nd_min, meantepfr_nd_max])
    axis([0, 2, 0.1, 0.32])
    yticks(0.1:0.05:0.3)
end

%tepfr - LSF_2_LE_p2
subplot(1,9,8)
title('LSF: 2, LE: 0.2')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p2_Files),')'])
if NumOf_LSF_2_LE_p2_Files > 0
    distributionPlot(bl_LSF_2_LE_p2.*mean_nd_tepfr_LSF_2_LE_p2,'addSpread',true)
%     axis([0, 2, meantepfr_nd_min, meantepfr_nd_max])
    axis([0, 2, 0.1, 0.32])
    yticks(0.1:0.05:0.3)
end

%tepfr - LSF_2_LE_p4
subplot(1,9,9)
title('LSF: 2, LE: 0.4')
xlabel(['(n=',num2str(NumOf_LSF_2_LE_p4_Files),')'])
if NumOf_LSF_2_LE_p4_Files > 0
    distributionPlot(bl_LSF_2_LE_p4.*mean_nd_tepfr_LSF_2_LE_p4,'addSpread',true)
%     axis([0, 2, meantepfr_nd_min, meantepfr_nd_max])
    axis([0, 2, 0.1, 0.32])
    yticks(0.1:0.05:0.3)
end

set(fig16, 'Position', [120,120, 1000, 600])

saveas(gcf,'../Figures_MPC/compSize/Fig29b_sishManip_TEpfr.fig');
saveas(gcf,'../Figures_MPC/compSize/Fig29b_sishManip_TEpfr.jpg');

disp('Figure 29b done')
stop
%% Sample Tracking error figure

figure(17)
plot(trackingerror_fa_LSF_1_LE_0(1,1:(end-1)),'-k','LineWidth', 2)
% hold on;
% plot(numel(trackingerror_fa_LSF_1_LE_0(1,1:(end-1)))/2,...
%     mean(trackingerror_fa_LSF_1_LE_0(1,1:(end-1))),'or','MarkerSize',20,...
%     'MarkerFaceColor','r')
ylabel('Tracking error (in cm)')
axis([0, inf, 0, 1])
yticks(0:0.25:1)

saveas(gcf,'../Figures_MPC/compSize/teSample_1_LSF_1_LE_0.jpg');
%%
figure(18)
plot(trackingerror_fa_LSF_1_LE_0(1,1:(end-1)),'Color', [0.55, 0.55, 0.55],...
    'LineWidth', 2)
hold on;
plot(numel(trackingerror_fa_LSF_1_LE_0(1,1:(end-1)))/2,...
    mean(trackingerror_fa_LSF_1_LE_0(1,1:(end-1))),'or','MarkerSize',16,...
    'MarkerFaceColor','r')
ylabel('Tracking error (in cm)')
axis([0, inf, 0, 1])
yticks(0:0.25:1)

saveas(gcf,'../Figures_MPC/compSize/teSample_2_LSF_1_LE_0.jpg');

%% Stats for morphological modifications, vertical oscillating signal

%General form is...
%[p,tbl,stats] = kruskalwallis(x,group)
stop
%% non-dim. tracking error per full run, vertically oscillating, morphMods
[p_ndTEpfr,tbl_ndTEpfr,stats_ndTEpfr] = ...
    kruskalwallis(nd_tepfr_export_numbers,export_sishManip);
figure;
c_ndTEpfr = multcompare(stats_ndTEpfr,'Alpha',0.05,'CType','bonferroni');

%% mech. work per full run, vertically oscillating, morphMods
[p_wpfr,tbl_wpfr,stats_wpfr] = ...
    kruskalwallis(wpfr_export_numbers,export_sishManip);
figure;
c_wpfr = multcompare(stats_wpfr,'Alpha',0.05,'CType','bonferroni');

%% Cost of Transport per full run, vertically oscillating, morphMods
[p_nd_wpfr,tbl_nd_wpfr,stats_nd_wpfr] = ...
    kruskalwallis(nd_wpfr_export_numbers,export_sishManip);
figure;
c_nd_wpfr = multcompare(stats_nd_wpfr,'Alpha',0.05,'CType','bonferroni');

%% Final window
disp('Completely done, son')
disp(['Time at end: ',datestr(now)])