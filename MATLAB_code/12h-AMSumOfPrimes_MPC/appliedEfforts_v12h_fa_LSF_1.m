%7/9/2018 Script developed to plot the F, alpha, tau distributions
%8/23/2018 Script modified for LengthScaleFactor (LSF)
%10/23/2018 Script modified for size extension (LE)
%1/11/19 Script modified for Model Predictive Control (MPC)
    %1/16/19 Script corrected previous indexing error.
%6/18/19 Script modified to compare between and within treatments 

%% We the People, in Order to form a perfect Union,...
clc;
clearvars;
close all;
disp(['Time at start: ',datestr(now)])

%% List the directory stuff
listOf_fa_LSF_1_LE_0_Files = dir('../SimData_MPC/fa_Sims/LSF_1/Winstore_*_fa_LSF_1_LE_0.mat'); 
listOf_fa_LSF_1_LE_p2_Files = dir('../SimData_MPC/fa_Sims/LSF_1/Winstore_*_fa_LSF_1_LE_p2.mat'); 
listOf_fa_LSF_1_LE_p4_Files = dir('../SimData_MPC/fa_Sims/LSF_1/Winstore_*_fa_LSF_1_LE_p4.mat'); 

LSF_1 = 1;

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

NumOf_fa_LSF_1_LE_0_Files = numel(listOf_fa_LSF_1_LE_0_Files);
NumOf_fa_LSF_1_LE_p2_Files = numel(listOf_fa_LSF_1_LE_p2_Files);
NumOf_fa_LSF_1_LE_p4_Files = numel(listOf_fa_LSF_1_LE_p4_Files);

%% Import all the variables
disp('Importing the relevant data')

%Magnitude of applied forces
load('../CuratedData_MPC/fa/LSF_1/F_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/F_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/F_fa_LSF_1_LE_p4.mat')

%Direction of applied forces
load('../CuratedData_MPC/fa/LSF_1/alpha_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/alpha_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/alpha_fa_LSF_1_LE_p4.mat')

%Abdominal torque
load('../CuratedData_MPC/fa/LSF_1/tauAbdo_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/tauAbdo_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/tauAbdo_fa_LSF_1_LE_p4.mat')

%Wing torque
load('../CuratedData_MPC/fa/LSF_1/tauWing_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/tauWing_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/tauWing_fa_LSF_1_LE_p4.mat')

%% Find the min & max values for...
%Forces
Fmax_vec = [max(max(F_fa_LSF_1_LE_0)),...
    max(max(F_fa_LSF_1_LE_p2)),...
    max(max(F_fa_LSF_1_LE_p4))];
F_max = max(Fmax_vec);

Fmin_vec = [min(min(F_fa_LSF_1_LE_0)),...
    min(min(F_fa_LSF_1_LE_p2)),...
    min(min(F_fa_LSF_1_LE_p4))];
F_min = min(Fmin_vec);

%Abdominal torque
tauAbdoMax_vec = [max(max(tauAbdo_fa_LSF_1_LE_0)),...
    max(max(tauAbdo_fa_LSF_1_LE_p2)),...
    max(max(tauAbdo_fa_LSF_1_LE_p4))];
tauAbdo_max = max(tauAbdoMax_vec);

tauAbdoMin_vec = [min(min(tauAbdo_fa_LSF_1_LE_0)),...
    min(min(tauAbdo_fa_LSF_1_LE_p2)),...
    min(min(tauAbdo_fa_LSF_1_LE_p4))];
tauAbdo_min = min(tauAbdoMin_vec);

%Wing torque
tauWingMax_vec = [max(max(tauWing_fa_LSF_1_LE_0)),...
    max(max(tauWing_fa_LSF_1_LE_p2)),...
    max(max(tauWing_fa_LSF_1_LE_p4))];
tauWing_max = max(tauWingMax_vec);

tauWingMin_vec = [min(min(tauWing_fa_LSF_1_LE_0)),...
    min(min(tauWing_fa_LSF_1_LE_p2)),...
    min(min(tauWing_fa_LSF_1_LE_p4))];
tauWing_min = min(tauWingMin_vec);

%% Conversion to mks

%Units of force WERE in g*cm/(s^2), now converted to 10^(-5) Newtons
conv_F_mks = 1/(1000*100); 

%Units of torque WERE in XX, now converted to N*m
conv_torque_mks = 1/((100^2)*1000);

zeroline(1:numel(time_hws)) = 0;

%% Magnitude of applied force plot w.r.t. time
fig11 = figure(11);

%y-motion
subplot(4,1,1)
plot(Tstore(1:(end-100)),y_g,'-k','LineWidth',2)
title(['fa, LSF: ',num2str(LSF_1)]);
ylabel('y goal (cm)'); 
xlabel('Time (s)');
axis([0, Tstore(end), -11, 11])

%fa - LSF_1_LE_0
subplot(4,1,2)
if NumOf_fa_LSF_1_LE_0_Files > 0
    plot(time_hws,F_fa_LSF_1_LE_0.*conv_F_mks,'ok','MarkerSize',2)
    axis([0, Tstore(end), F_min*conv_F_mks, F_max*conv_F_mks])
end
ylabel('Force (N)');
xlabel(['LE: 0, Time (s), (n = ',num2str(NumOf_fa_LSF_1_LE_0_Files),')']);

%fa - LSF_1_LE_p2
subplot(4,1,3)
if NumOf_fa_LSF_1_LE_p2_Files > 0
    plot(time_hws,F_fa_LSF_1_LE_p2.*conv_F_mks,'ok','MarkerSize',2)
    axis([0, Tstore(end), F_min*conv_F_mks, F_max*conv_F_mks])
end
ylabel('Force (N)'); 
xlabel(['LE: 0.2, Time (s), (n = ',num2str(NumOf_fa_LSF_1_LE_p2_Files),')']);

%fa - LSF_1_LE_p4
subplot(4,1,4)
if NumOf_fa_LSF_1_LE_p4_Files > 0
    plot(time_hws,F_fa_LSF_1_LE_p4.*conv_F_mks,'ok','MarkerSize',2)
    axis([0, Tstore(end), F_min*conv_F_mks, F_max*conv_F_mks])
end
ylabel('Force (N)');
xlabel(['LE: 0.4, Time (s), (n = ',num2str(NumOf_fa_LSF_1_LE_p4_Files),')']);

set(fig11, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig11_F_fa_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig11_F_fa_LSF_1_v12h.jpg');

disp('Figure 11 done')

%% Direction of applied force plot w.r.t. time
fig12 = figure(12);

%y-motion
subplot(4,1,1)
plot(Tstore(1:(end-100)),y_g,'-k','LineWidth',2)
title(['fa, LSF: ',num2str(LSF_1)]);
ylabel('y goal (cm)');  
xlabel('Time (s)');
axis([0, Tstore(end), -11, 11])

%fa - LSF_1_LE_0
subplot(4,1,2)
if NumOf_fa_LSF_1_LE_0_Files > 0
    plot(time_hws,alpha_fa_LSF_1_LE_0,'ok','MarkerSize',2)
    axis([0, Tstore(end), 0, 2*pi])
end
ylabel('Alpha (rad)'); 
xlabel(['LE: 0, Time (s), (n = ',num2str(NumOf_fa_LSF_1_LE_0_Files),')']);

%fa - LSF_1_LE_p2
subplot(4,1,3)
if NumOf_fa_LSF_1_LE_p2_Files > 0
    plot(time_hws,alpha_fa_LSF_1_LE_p2,'ok','MarkerSize',2)
    axis([0, Tstore(end), 0, 2*pi])
end
ylabel('Alpha (rad)');
xlabel(['LE: 0.2, Time (s), (n = ',num2str(NumOf_fa_LSF_1_LE_p2_Files),')']);

%fa - LSF_1_LE_p4
subplot(4,1,4)
if NumOf_fa_LSF_1_LE_p4_Files > 0
    plot(time_hws,alpha_fa_LSF_1_LE_p4,'ok','MarkerSize',2)
    axis([0, Tstore(end), 0, 2*pi])
end
ylabel('Alpha (rad)');
xlabel(['LE: 0.4, Time (s), (n = ',num2str(NumOf_fa_LSF_1_LE_p4_Files),')']);

set(fig12, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig12_alpha_fa_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig12_alpha_fa_LSF_1_v12h.jpg');

disp('Figure 12 done')

%% Abdominal torque w.r.t. time
fig13 = figure(13);

%y-motion
subplot(4,1,1)
plot(Tstore(1:(end-100)),y_g,'-k','LineWidth',2)
ylabel('y goal (cm)'); 
title(['fa, LSF: ',num2str(LSF_1)]); 
xlabel('Time (s)');
axis([0, Tstore(end), -11, 11])

%fa - LSF_1_LE_0
subplot(4,1,2)
if NumOf_fa_LSF_1_LE_0_Files > 0
    plot(time_hws,tauAbdo_fa_LSF_1_LE_0.*conv_torque_mks,'ok','MarkerSize',2)
    axis([0, Tstore(end), tauAbdo_min*conv_torque_mks, tauAbdo_max*conv_torque_mks])
end
hold on;
plot(time_hws,zeroline,'-r','LineWidth',1)
ylabel({'Abdominal torque';'(N*m)'}); 
xlabel(['LE: 0, Time(s), (n = ',num2str(NumOf_fa_LSF_1_LE_0_Files),')']);

%fa - LSF_1_LE_p2
subplot(4,1,3)
if NumOf_fa_LSF_1_LE_p2_Files > 0
    plot(time_hws,tauAbdo_fa_LSF_1_LE_p2.*conv_torque_mks,'ok','MarkerSize',2)
    axis([0, Tstore(end), tauAbdo_min*conv_torque_mks, tauAbdo_max*conv_torque_mks])
end
hold on;
plot(time_hws,zeroline,'-r','LineWidth',1)
ylabel({'Abdominal torque';'(N*m)'}); 
xlabel(['LE: 0.2, Time(s), (n = ',num2str(NumOf_fa_LSF_1_LE_p2_Files),')']);

%fa - LSF_1_LE_p4
subplot(4,1,4)
if NumOf_fa_LSF_1_LE_p4_Files > 0
    plot(time_hws,tauAbdo_fa_LSF_1_LE_p4.*conv_torque_mks,'ok','MarkerSize',2)
    axis([0, Tstore(end), tauAbdo_min*conv_torque_mks, tauAbdo_max*conv_torque_mks])
end
hold on;
plot(time_hws,zeroline,'-r','LineWidth',1)
ylabel({'Abdominal torque';'(N*m)'}); 
xlabel(['LE: 0.4, Time(s), (n = ',num2str(NumOf_fa_LSF_1_LE_p4_Files),')']);

set(fig13, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig13_tauAbdo_fa_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig13_tauAbdo_fa_LSF_1_v12h.jpg');

disp('Figure 13 done')

%% Wing torque w.r.t. time
fig14 = figure(14);

%y-motion
subplot(4,1,1)
plot(Tstore(1:(end-100)),y_g,'-k','LineWidth',2)
ylabel('y goal (cm)'); 
title(['fa, LSF: ',num2str(LSF_1)]); 
xlabel('Time (s)');
axis([0, Tstore(end), -11, 11])

%fa - LSF_1_LE_0
subplot(4,1,2)
if NumOf_fa_LSF_1_LE_0_Files > 0
    plot(time_hws,tauWing_fa_LSF_1_LE_0.*conv_torque_mks,'ok','MarkerSize',2)
    axis([0, Tstore(end), tauWing_min*conv_torque_mks, tauWing_max*conv_torque_mks])
end
hold on;
plot(time_hws,zeroline,'-r','LineWidth',1)
ylabel({'Wing torque'; '(N*m)'}); 
xlabel(['LE: 0, Time(s), (n = ',num2str(NumOf_fa_LSF_1_LE_0_Files),')']);

%fa - LSF_1_LE_p2
subplot(4,1,3)
if NumOf_fa_LSF_1_LE_p2_Files > 0
    plot(time_hws,tauWing_fa_LSF_1_LE_p2.*conv_torque_mks,'ok','MarkerSize',2)
    axis([0, Tstore(end), tauWing_min*conv_torque_mks, tauWing_max*conv_torque_mks])
end
hold on;
plot(time_hws,zeroline,'-r','LineWidth',1)
ylabel({'Wing torque'; '(N*m)'}); 
xlabel(['LE: 0.2, Time(s), (n = ',num2str(NumOf_fa_LSF_1_LE_p2_Files),')']);

%fa - LSF_1_LE_p4
subplot(4,1,4)
if NumOf_fa_LSF_1_LE_p4_Files > 0
    plot(time_hws,tauWing_fa_LSF_1_LE_p4.*conv_torque_mks,'ok','MarkerSize',2)
    axis([0, Tstore(end), tauWing_min*conv_torque_mks, tauWing_max*conv_torque_mks])
end
hold on;
plot(time_hws,zeroline,'-r','LineWidth',1)
ylabel({'Wing torque'; '(N*m)'}); 
xlabel(['LE: 0.4, Time(s), (n = ',num2str(NumOf_fa_LSF_1_LE_p4_Files),')']);

set(fig14, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig14_tauWing_fa_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig14_tauWing_fa_LSF_1_v12h.jpg');

disp('Figure 14 done')

%% Force w.r.t. alpha
fig15 = figure(15);

%fa - LSF_1_LE_0
subplot(3,1,1)
if NumOf_fa_LSF_1_LE_0_Files > 0
    plot(alpha_fa_LSF_1_LE_0,F_fa_LSF_1_LE_0.*conv_F_mks,'ok','MarkerSize',2)
    axis([0, 2*pi, F_min*conv_F_mks, F_max*conv_F_mks])
end
title(['fa, LSF: ',num2str(LSF_1)]); 
ylabel('Force (N)'); 
xlabel(['LE: 0, Alpha (rad), (n = ',num2str(NumOf_fa_LSF_1_LE_0_Files),')']);

%fa - LSF_1_LE_p2
subplot(3,1,2)
if NumOf_fa_LSF_1_LE_p2_Files > 0
    plot(alpha_fa_LSF_1_LE_p2,F_fa_LSF_1_LE_p2.*conv_F_mks,'ok','MarkerSize',2)
    axis([0, 2*pi, F_min*conv_F_mks, F_max*conv_F_mks])
end
ylabel('Force (N)');
xlabel(['LE: 0.2, Alpha (rad), (n = ',num2str(NumOf_fa_LSF_1_LE_p2_Files),')']);

%fa - LSF_1_LE_p4
subplot(3,1,3)
if NumOf_fa_LSF_1_LE_p4_Files > 0
    plot(alpha_fa_LSF_1_LE_p4,F_fa_LSF_1_LE_p4.*conv_F_mks,'ok','MarkerSize',2)
    axis([0, 2*pi, F_min*conv_F_mks, F_max*conv_F_mks])
end
ylabel('Force (N)');
xlabel(['LE: 0.4, Alpha (rad), (n = ',num2str(NumOf_fa_LSF_1_LE_p4_Files),')']);

set(fig15, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig15_F_alpha_fa_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig15_F_alpha_fa_LSF_1_v12h.jpg');

disp('Figure 15 done')

%% Applied force w.r.t. abdominal torque
fig16 = figure(16);

%fa - LSF_1_LE_0
subplot(3,1,1)
if NumOf_fa_LSF_1_LE_0_Files > 0
    plot(tauAbdo_fa_LSF_1_LE_0.*conv_torque_mks,F_fa_LSF_1_LE_0.*conv_F_mks,'ok','MarkerSize',2)
    axis([tauAbdo_min*conv_torque_mks, tauAbdo_max*conv_torque_mks, F_min*conv_F_mks, F_max*conv_F_mks])
end
title(['fa, LSF: ',num2str(LSF_1)]);
ylabel('Force (N)'); 
xlabel(['LE: 0, Abdominal torque (N*m), (n = ',num2str(NumOf_fa_LSF_1_LE_0_Files),')']);

%fa - LSF_1_LE_p2
subplot(3,1,2)
if NumOf_fa_LSF_1_LE_p2_Files > 0
    plot(tauAbdo_fa_LSF_1_LE_p2.*conv_torque_mks,F_fa_LSF_1_LE_p2.*conv_F_mks,'ok','MarkerSize',2)
    axis([tauAbdo_min*conv_torque_mks, tauAbdo_max*conv_torque_mks, F_min*conv_F_mks, F_max*conv_F_mks])
end
ylabel('Force (N)');
xlabel(['LE: 0.2, Abdominal torque (N*m), (n = ',num2str(NumOf_fa_LSF_1_LE_p2_Files),')']);

%fa - LSF_1_LE_p4
subplot(3,1,3)
if NumOf_fa_LSF_1_LE_p4_Files > 0
    plot(tauAbdo_fa_LSF_1_LE_p4.*conv_torque_mks,F_fa_LSF_1_LE_p4.*conv_F_mks,'ok','MarkerSize',2)
    axis([tauAbdo_min*conv_torque_mks, tauAbdo_max*conv_torque_mks, F_min*conv_F_mks, F_max*conv_F_mks])
end
ylabel('Force (N)');
xlabel(['LE: 0.4, Abdominal torque (N*m), (n = ',num2str(NumOf_fa_LSF_1_LE_p4_Files),')']);

set(fig16, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig16_F_tauAbdo_fa_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig16_F_tauAbdo_fa_LSF_1_v12h.jpg');

disp('Figure 16 done')

%% Applied force w.r.t. wing torque
fig17 = figure(17);

%fa - LSF_1_LE_0
subplot(3,1,1)
if NumOf_fa_LSF_1_LE_0_Files > 0
    plot(tauWing_fa_LSF_1_LE_0.*conv_torque_mks,F_fa_LSF_1_LE_0.*conv_F_mks,'ok','MarkerSize',2)
    axis([tauWing_min*conv_torque_mks, tauWing_max*conv_torque_mks, F_min*conv_F_mks, F_max*conv_F_mks])
end
title(['fa, LSF: ',num2str(LSF_1)]);
ylabel('Force (N)'); 
xlabel(['LE: 0, Wing torque (N*m), (n = ',num2str(NumOf_fa_LSF_1_LE_0_Files),')']);

%fa - LSF_1_LE_p2
subplot(3,1,2)
if NumOf_fa_LSF_1_LE_p2_Files > 0
    plot(tauWing_fa_LSF_1_LE_p2.*conv_torque_mks,F_fa_LSF_1_LE_p2.*conv_F_mks,'ok','MarkerSize',2)
    axis([tauWing_min*conv_torque_mks, tauWing_max*conv_torque_mks, F_min*conv_F_mks, F_max*conv_F_mks])
end
ylabel('Force (N)'); 
xlabel(['LE: 0.2, Wing torque (N*m), (n = ',num2str(NumOf_fa_LSF_1_LE_p2_Files),')']);

%fa - LSF_1_LE_p4
subplot(3,1,3)
if NumOf_fa_LSF_1_LE_p4_Files > 0
    plot(tauWing_fa_LSF_1_LE_p4.*conv_torque_mks,F_fa_LSF_1_LE_p4.*conv_F_mks,'ok','MarkerSize',2)
    axis([tauWing_min*conv_torque_mks, tauWing_max*conv_torque_mks, F_min*conv_F_mks, F_max*conv_F_mks])
end
ylabel('Force (N)'); 
xlabel(['LE: 0.4, Wing torque (N*m), (n = ',num2str(NumOf_fa_LSF_1_LE_p4_Files),')']);

set(fig17, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig17_F_tauWing_fa_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig17_F_tauWing_fa_LSF_1_v12h.jpg');

disp('Figure 17 done')

%% Abdominal torque w.r.t. wing torque
fig18 = figure(18);

%fa - LSF_1_LE_0
subplot(3,1,1)
if NumOf_fa_LSF_1_LE_0_Files > 0
    plot(tauWing_fa_LSF_1_LE_0.*conv_torque_mks,tauAbdo_fa_LSF_1_LE_0.*conv_torque_mks,'ok','MarkerSize',2)
    axis([tauWing_min*conv_torque_mks, tauWing_max*conv_torque_mks, ...
        tauAbdo_min*conv_torque_mks, tauAbdo_max*conv_torque_mks])
end
title(['fa, LSF: ',num2str(LSF_1)]); 
ylabel({'Abdominal torque';'(N*m)'}); 
xlabel(['LE: 0, Wing torque (N*m), (n = ',num2str(NumOf_fa_LSF_1_LE_0_Files),')']);


%fa - LSF_1_LE_p2
subplot(3,1,2)
if NumOf_fa_LSF_1_LE_p2_Files > 0
    plot(tauWing_fa_LSF_1_LE_p2.*conv_torque_mks,tauAbdo_fa_LSF_1_LE_p2.*conv_torque_mks,'ok','MarkerSize',2)
    axis([tauWing_min*conv_torque_mks, tauWing_max*conv_torque_mks, ...
        tauAbdo_min*conv_torque_mks, tauAbdo_max*conv_torque_mks])
end
ylabel({'Abdominal torque';'(N*m)'}); 
xlabel(['LE: 0.2, Wing torque (N*m), (n = ',num2str(NumOf_fa_LSF_1_LE_p2_Files),')']);


%fa - LSF_1_LE_p4
subplot(3,1,3)
if NumOf_fa_LSF_1_LE_p4_Files > 0
    plot(tauWing_fa_LSF_1_LE_p4.*conv_torque_mks,tauAbdo_fa_LSF_1_LE_p4.*conv_torque_mks,'ok','MarkerSize',2)
    axis([tauWing_min*conv_torque_mks, tauWing_max*conv_torque_mks, ...
        tauAbdo_min*conv_torque_mks, tauAbdo_max*conv_torque_mks])
end
ylabel({'Abdominal torque';'(N*m)'}); 
xlabel(['LE: 0.4, Wing torque (N*m), (n = ',num2str(NumOf_fa_LSF_1_LE_p4_Files),')']);

set(fig18, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig18_tauAbdo_tauWing_fa_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig18_tauAbdo_tauWing_fa_LSF_1_v12h.jpg');

disp('Figure 18 done')

%% Abdominal torque w.r.t. alpha
fig19 = figure(19);

%fa - LSF_1_LE_0
subplot(3,1,1)
if NumOf_fa_LSF_1_LE_0_Files > 0
    plot(alpha_fa_LSF_1_LE_0,tauAbdo_fa_LSF_1_LE_0.*conv_torque_mks,'ok','MarkerSize',2)
    axis([0, 2*pi, tauAbdo_min*conv_torque_mks, tauAbdo_max*conv_torque_mks])
end
title(['fa, LSF: ',num2str(LSF_1)]); 
ylabel({'Abdominal torque';'(N*m)'}); 
xlabel(['LE: 0, Alpha (rad), (n = ',num2str(NumOf_fa_LSF_1_LE_0_Files),')']);

%fa - LSF_1_LE_p2
subplot(3,1,2)
if NumOf_fa_LSF_1_LE_p2_Files > 0
    plot(alpha_fa_LSF_1_LE_p2,tauAbdo_fa_LSF_1_LE_p2.*conv_torque_mks,'ok','MarkerSize',2)
    axis([0, 2*pi, tauAbdo_min*conv_torque_mks, tauAbdo_max*conv_torque_mks])
end
ylabel({'Abdominal torque';'(N*m)'});
xlabel(['LE: 0.2, Alpha (rad), (n = ',num2str(NumOf_fa_LSF_1_LE_p2_Files),')']);

%fa - LSF_1_LE_p4
subplot(3,1,3)
if NumOf_fa_LSF_1_LE_p4_Files > 0
    plot(alpha_fa_LSF_1_LE_p4,tauAbdo_fa_LSF_1_LE_p4.*conv_torque_mks,'ok','MarkerSize',2)
    axis([0, 2*pi, tauAbdo_min*conv_torque_mks, tauAbdo_max*conv_torque_mks])
end
ylabel({'Abdominal torque';'(N*m)'});
xlabel(['LE: 0.4, Alpha (rad), (n = ',num2str(NumOf_fa_LSF_1_LE_p4_Files),')']);

set(fig19, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig19_tauAbdo_alpha_fa_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig19_tauAbdo_alpha_fa_LSF_1_v12h.jpg');

disp('Figure 19 done')

%% Wing torque w.r.t. alpha
fig20 = figure(20);

%fa - LSF_1_LE_0
subplot(3,1,1)
if NumOf_fa_LSF_1_LE_0_Files > 0
    plot(alpha_fa_LSF_1_LE_0,tauWing_fa_LSF_1_LE_0.*conv_torque_mks,'ok','MarkerSize',2)
    axis([0, 2*pi, tauWing_min*conv_torque_mks, tauWing_max*conv_torque_mks])
end
title(['fa, LSF: ',num2str(LSF_1)]); 
ylabel({'Wing torque'; '(N*m)'}); 
xlabel(['LE: 0, Alpha (rad), (n = ',num2str(NumOf_fa_LSF_1_LE_0_Files),')']);

%fa - LSF_1_LE_p2
subplot(3,1,2)
if NumOf_fa_LSF_1_LE_p2_Files > 0
    plot(alpha_fa_LSF_1_LE_p2,tauWing_fa_LSF_1_LE_p2.*conv_torque_mks,'ok','MarkerSize',2)
    axis([0, 2*pi, tauWing_min*conv_torque_mks, tauWing_max*conv_torque_mks])
end
ylabel({'Wing torque'; '(N*m)'}); 
xlabel(['LE: 0.2, Alpha (rad), (n = ',num2str(NumOf_fa_LSF_1_LE_p2_Files),')']);

%fa - LSF_1_LE_p4
subplot(3,1,3)
if NumOf_fa_LSF_1_LE_p4_Files > 0
    plot(alpha_fa_LSF_1_LE_p4,tauWing_fa_LSF_1_LE_p4.*conv_torque_mks,'ok','MarkerSize',2)
    axis([0, 2*pi, tauWing_min*conv_torque_mks, tauWing_max*conv_torque_mks])
end
ylabel({'Wing torque'; '(N*m)'}); 
xlabel(['LE: 0.4, Alpha (rad), (n = ',num2str(NumOf_fa_LSF_1_LE_p4_Files),')']);

set(fig20, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig20_tauWing_alpha_fa_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig20_tauWing_alpha_fa_LSF_1_v12h.jpg');

disp('Figure 20 done')

%% 3 axis plot - Force, alpha, Abdominal torque
fig21 = figure(21);

%fa - LSF_1_LE_0
subplot(2,2,1)
if NumOf_fa_LSF_1_LE_0_Files > 0
    plot3(alpha_fa_LSF_1_LE_0,F_fa_LSF_1_LE_0.*conv_F_mks,...
        tauAbdo_fa_LSF_1_LE_0.*conv_torque_mks,'ok','MarkerSize',2)
    axis([0, 2*pi, F_min*conv_F_mks, F_max*conv_F_mks,...
        tauAbdo_min*conv_torque_mks, tauAbdo_max*conv_torque_mks])
end
title({'fa'; ['LSF: ',num2str(LSF_1)];'LE: 0'});
zlabel('Abdominal torque (N*m)'); 
ylabel('Force (N)');
xlabel(['Alpha (rad), (n = ',num2str(NumOf_fa_LSF_1_LE_0_Files),')']);

%fa - LSF_1_LE_p2
subplot(2,2,2)
if NumOf_fa_LSF_1_LE_p2_Files > 0
    plot3(alpha_fa_LSF_1_LE_p2,F_fa_LSF_1_LE_p2.*conv_F_mks,...
        tauAbdo_fa_LSF_1_LE_p2.*conv_torque_mks,'ok','MarkerSize',2)
    axis([0, 2*pi, F_min*conv_F_mks, F_max*conv_F_mks,...
        tauAbdo_min*conv_torque_mks, tauAbdo_max*conv_torque_mks])
end
title({'fa'; ['LSF: ',num2str(LSF_1)];'LE: 0.2'});
zlabel('Abdominal torque (N*m)'); 
ylabel('Force (N)');
xlabel(['Alpha (rad), (n = ',num2str(NumOf_fa_LSF_1_LE_p2_Files),')']);


%fa - LSF_1_LE_p4
subplot(2,2,4)
if NumOf_fa_LSF_1_LE_p4_Files > 0
    plot3(alpha_fa_LSF_1_LE_p4,F_fa_LSF_1_LE_p4.*conv_F_mks,...
        tauAbdo_fa_LSF_1_LE_p4.*conv_torque_mks,'ok','MarkerSize',2)
    axis([0, 2*pi, F_min*conv_F_mks, F_max*conv_F_mks,...
        tauAbdo_min*conv_torque_mks, tauAbdo_max*conv_torque_mks])
end
title({'fa'; ['LSF: ',num2str(LSF_1)];'LE: 0.4'});
zlabel('Abdominal torque (N*m)'); 
ylabel('Force (N)');
xlabel(['Alpha (rad), (n = ',num2str(NumOf_fa_LSF_1_LE_p4_Files),')']);

set(fig21, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig21_force_alpha_tauAbdo_fa_LSF_1_v12h.fig');
% saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig21_force_alpha_tauAbdo_fa_LSF_1_v12h.jpg');

disp('Figure 21 done')

%% 3 axis plot - Force, alpha, Wing torque
fig22 = figure(22);

%fa - LSF_1_LE_0
subplot(2,2,1)
if NumOf_fa_LSF_1_LE_0_Files > 0
    plot3(alpha_fa_LSF_1_LE_0,F_fa_LSF_1_LE_0.*conv_F_mks,...
        tauWing_fa_LSF_1_LE_0.*conv_torque_mks,'ok','MarkerSize',2)
    axis([0, 2*pi, F_min*conv_F_mks, F_max*conv_F_mks,...
        tauWing_min*conv_torque_mks, tauWing_max*conv_torque_mks])
end
title({'fa'; ['LSF: ',num2str(LSF_1)];'LE: 0'});
zlabel('Wing torque (N*m)'); 
ylabel('Force (N)');
xlabel(['Alpha (rad), (n = ',num2str(NumOf_fa_LSF_1_LE_0_Files),')']);

%fa - LSF_1_LE_p2
subplot(2,2,2)
if NumOf_fa_LSF_1_LE_p2_Files > 0
    plot3(alpha_fa_LSF_1_LE_p2,F_fa_LSF_1_LE_p2.*conv_F_mks,...
        tauWing_fa_LSF_1_LE_p2.*conv_torque_mks,'ok','MarkerSize',2)
    axis([0, 2*pi, F_min*conv_F_mks, F_max*conv_F_mks,...
        tauWing_min*conv_torque_mks, tauWing_max*conv_torque_mks])
end
title({'fa'; ['LSF: ',num2str(LSF_1)];'LE: 0.2'});
zlabel('Wing torque (N*m)'); 
ylabel('Force (N)');
xlabel(['Alpha (rad), (n = ',num2str(NumOf_fa_LSF_1_LE_p2_Files),')']);

%fa - LSF_1_LE_p4
subplot(2,2,4)
if NumOf_fa_LSF_1_LE_p4_Files > 0
    plot3(alpha_fa_LSF_1_LE_p4,F_fa_LSF_1_LE_p4.*conv_F_mks,...
        tauWing_fa_LSF_1_LE_p4.*conv_torque_mks,'ok','MarkerSize',2)
    axis([0, 2*pi, F_min*conv_F_mks, F_max*conv_F_mks,...
        tauWing_min*conv_torque_mks, tauWing_max*conv_torque_mks])
end
title({'fa'; ['LSF: ',num2str(LSF_1)];'LE: 0.4'});
zlabel('Wing torque (N*m)'); 
ylabel('Force (N)');
xlabel(['Alpha (rad), (n = ',num2str(NumOf_fa_LSF_1_LE_p4_Files),')']);

set(fig22, 'Position', [120,120,850,550])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig22_force_alpha_tauWing_fa_LSF_1_v12h.fig');
% saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig22_force_alpha_tauWing_fa_LSF_1_v12h.jpg');

disp('Figure 22 done')

%% Individualized abdominal torque w.r.t. time, LE: 0
fig23 = figure(23);

%fa - LSF_1_LE_0
for i = 1:NumOf_fa_LSF_1_LE_0_Files
    subplot(5,8,i)
    plot(time_hws,tauAbdo_fa_LSF_1_LE_0(i,:).*conv_torque_mks,'ok','MarkerSize',2)
    hold on;
    plot(time_hws,zeroline,'-r','LineWidth',1) 
    xlabel(['#',num2str(i),', Time(s)']);
    if i == 1
        title({'LSF: 1';'fa, LE: 0'})
        ylabel({'Abdominal torque';'(N*m)'});
    end
    axis([0, Tstore(end), tauAbdo_min*conv_torque_mks, tauAbdo_max*conv_torque_mks])
end

set(fig23, 'Position', [120,120,1200,700])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig23_tauAbdoIndiv_fa_LSF_1_LE_0_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig23_tauAbdoIndiv_fa_LSF_1_LE_0_v12h.jpg');

disp('Figure 23 done')

%% Individualized abdominal torque w.r.t. time, LE: 0.2
fig24 = figure(24);

%fa - LSF_1_LE_p2
for i = 1:NumOf_fa_LSF_1_LE_p2_Files
    subplot(5,8,i)
    plot(time_hws,tauAbdo_fa_LSF_1_LE_p2(i,:).*conv_torque_mks,'ok','MarkerSize',2)
    hold on;
    plot(time_hws,zeroline,'-r','LineWidth',1) 
    xlabel(['#',num2str(i),', Time(s)']);
    if i == 1
        title({'LSF: 1';'fa, LE: 0.2'})
        ylabel({'Abdominal torque';'(N*m)'});
    end
    axis([0, Tstore(end), tauAbdo_min*conv_torque_mks, tauAbdo_max*conv_torque_mks])
end

set(fig24, 'Position', [120,120,1200,700])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig24_tauAbdoIndiv_fa_LSF_1_LE_p2_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig24_tauAbdoIndiv_fa_LSF_1_LE_p2_v12h.jpg');

disp('Figure 24 done')

%% Individualized abdominal torque w.r.t. time, LE: 0.4
fig25 = figure(25);

%fa - LSF_1_LE_p4
for i = 1:NumOf_fa_LSF_1_LE_p4_Files
    subplot(5,8,i)
    plot(time_hws,tauAbdo_fa_LSF_1_LE_p4(i,:).*conv_torque_mks,'ok','MarkerSize',2)
    hold on;
    plot(time_hws,zeroline,'-r','LineWidth',1) 
    xlabel(['#',num2str(i),', Time(s)']);
    if i == 1
        title({'LSF: 1';'fa, LE: 0.4'})
        ylabel({'Abdominal torque';'(N*m)'});
    end
    axis([0, Tstore(end), tauAbdo_min*conv_torque_mks, tauAbdo_max*conv_torque_mks])
end

set(fig25, 'Position', [120,120,1200,700])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig25_tauAbdoIndiv_fa_LSF_1_LE_p4_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig25_tauAbdoIndiv_fa_LSF_1_LE_p4_v12h.jpg');

disp('Figure 25 done')

%% Final window
disp('Completely done, son')
disp(['Time at end: ',datestr(now)])