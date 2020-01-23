%1/30/2018 Script developed to calculate the bode plot
%8/20/18 Script modified for LengthScaleFactor (LSF)
%10/17/18 Script modified for size extension (LE)
%1/11/19 Script modified for Model Predictive Control (MPC)
    %1/16/19 Script corrected previous indexing error.
%6/18/19 Script modified to compare between and within treatments 

%% We the People, in Order to form a perfect Union,...
clc;
clearvars;
close all;
disp(['Time at start: ',datestr(now)])

%% Assign moth variables IF NECESSARY
disp('Assign moth parameters')

%Length Scale Factors
LSF_1 = 1;

%% List the directory stuff 
listOf_fa_LSF_1_LE_0_Files = dir('../SimData_MPC/fa_Sims/LSF_1/Winstore_*_fa_LSF_1_LE_0.mat'); 
listOf_fa_LSF_1_LE_p2_Files = dir('../SimData_MPC/fa_Sims/LSF_1/Winstore_*_fa_LSF_1_LE_p2.mat'); 
listOf_fa_LSF_1_LE_p4_Files = dir('../SimData_MPC/fa_Sims/LSF_1/Winstore_*_fa_LSF_1_LE_p4.mat'); 
%the asterisk is a wildcard
%The dir function returns a "listing" of an M x 1 "structure." The 
%structure has five fields in this case listing: name, date, byte, isdir, 
%datenum.
%I used the wildcard because I know the number of text files will
%definitely increase as we gather more data.
%For more information enter   help dir   into MATLAB mainframe
disp('Load time vector')
load('../SimData_MPC/Tstore_MPC_hws_sp.mat')

NumOf_fa_LSF_1_LE_0_Files = numel(listOf_fa_LSF_1_LE_0_Files);
NumOf_fa_LSF_1_LE_p2_Files = numel(listOf_fa_LSF_1_LE_p2_Files);
NumOf_fa_LSF_1_LE_p4_Files = numel(listOf_fa_LSF_1_LE_p4_Files);

%% Assign the imported numbers to internal variables
disp('Assign internal variables')

hws = 500; %THIS IS 100 FOR 1c, 1d, AND 1f
timestep = numel(Tstore(1:(end-100)))/hws;
tderiv = Tstore(2);

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

%% Import all the variables
disp('Importing curated data')

%fa - LSF_1_LE_0
load('../CuratedData_MPC/fa/LSF_1/Winstore_y_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/trackingerror_fa_LSF_1_LE_0.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_beta_fa_LSF_1_LE_0.mat')

%fa - LSF_1_LE_p2
load('../CuratedData_MPC/fa/LSF_1/Winstore_y_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/trackingerror_fa_LSF_1_LE_p2.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_beta_fa_LSF_1_LE_p2.mat')

%fa - LSF_1_LE_p4
load('../CuratedData_MPC/fa/LSF_1/Winstore_y_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/trackingerror_fa_LSF_1_LE_p4.mat')
load('../CuratedData_MPC/fa/LSF_1/Winstore_beta_fa_LSF_1_LE_p4.mat')

disp('Curated data is loaded')

%% FFT values needed for the figures
%The following five lines of code are from 
%http://www.mathworks.com/help/matlab/ref/fft.html 
N = length(Tstore(1:(end-100))); %Total number of samples
Fs = N/10; %Tstore((end-100)) is 10, so I will hard code 10 seconds
%Sample frequency = total number of samples/total time
f_omega = Fs*((0:N-1))/N; %Actual frequency of the transformed signal
ampscale = length(y_g)/2; %Scaling the amplitude as necessary
maxfreq = 20; %in Hz; Maximum frequency you want to view in the figure
maxrange = find(f_omega == maxfreq);

disp('Assign relevant FFT values')

%% Fast Fourier transform of EACH signal
disp('Calculate the FFT of EACH signal')

y_sim_fa_LSF_1_LE_0 = zeros(NumOf_fa_LSF_1_LE_0_Files,hws*timestep);
y_sim_fa_LSF_1_LE_p2 = zeros(NumOf_fa_LSF_1_LE_p2_Files,hws*timestep);
y_sim_fa_LSF_1_LE_p4 = zeros(NumOf_fa_LSF_1_LE_p4_Files,hws*timestep);

te_sim_fa_LSF_1_LE_0 = zeros(NumOf_fa_LSF_1_LE_0_Files,hws*timestep);
te_sim_fa_LSF_1_LE_p2 = zeros(NumOf_fa_LSF_1_LE_p2_Files,hws*timestep);
te_sim_fa_LSF_1_LE_p4 = zeros(NumOf_fa_LSF_1_LE_p4_Files,hws*timestep);

beta_sim_fa_LSF_1_LE_0 = zeros(NumOf_fa_LSF_1_LE_0_Files,hws*timestep);
beta_sim_fa_LSF_1_LE_p2 = zeros(NumOf_fa_LSF_1_LE_p2_Files,hws*timestep);
beta_sim_fa_LSF_1_LE_p4 = zeros(NumOf_fa_LSF_1_LE_p4_Files,hws*timestep);

%y_goal FFT
y_theory = fft(y_g); 
te_theory = zeros(1,numel(y_g));

%FFT of the individual full runs
for i = 1:NumOf_fa_LSF_1_LE_0_Files
    y_sim_fa_LSF_1_LE_0(i,:) = fft(Winstore_y_fa_LSF_1_LE_0(i,1:hws*timestep)...
        - nanmean(Winstore_y_fa_LSF_1_LE_0(i,1:hws*timestep)));
    te_sim_fa_LSF_1_LE_0(i,:) = fft(trackingerror_fa_LSF_1_LE_0(i,1:hws*timestep)...
        - nanmean(trackingerror_fa_LSF_1_LE_0(i,1:hws*timestep)));
    beta_sim_fa_LSF_1_LE_0(i,:) = fft(Winstore_beta_fa_LSF_1_LE_0(i,1:hws*timestep)...
        - nanmean(Winstore_beta_fa_LSF_1_LE_0(i,1:hws*timestep)));
end

for i = 1:NumOf_fa_LSF_1_LE_p2_Files
    y_sim_fa_LSF_1_LE_p2(i,:) = fft(Winstore_y_fa_LSF_1_LE_p2(i,1:hws*timestep)...
        - nanmean(Winstore_y_fa_LSF_1_LE_p2(i,1:hws*timestep)));
    te_sim_fa_LSF_1_LE_p2(i,:) = fft(trackingerror_fa_LSF_1_LE_p2(i,1:hws*timestep)...
        - nanmean(trackingerror_fa_LSF_1_LE_p2(i,1:hws*timestep)));
    beta_sim_fa_LSF_1_LE_p2(i,:) = fft(Winstore_beta_fa_LSF_1_LE_p2(i,1:hws*timestep)...
        - nanmean(Winstore_beta_fa_LSF_1_LE_p2(i,1:hws*timestep)));
end

for i = 1:NumOf_fa_LSF_1_LE_p4_Files
    y_sim_fa_LSF_1_LE_p4(i,:) = fft(Winstore_y_fa_LSF_1_LE_p4(i,1:hws*timestep)...
        - nanmean(Winstore_y_fa_LSF_1_LE_p4(i,1:hws*timestep)));
    te_sim_fa_LSF_1_LE_p4(i,:) = fft(trackingerror_fa_LSF_1_LE_p4(i,1:hws*timestep)...
        - nanmean(trackingerror_fa_LSF_1_LE_p4(i,1:hws*timestep)));
    beta_sim_fa_LSF_1_LE_p4(i,:) = fft(Winstore_beta_fa_LSF_1_LE_p4(i,1:hws*timestep)...
        - nanmean(Winstore_beta_fa_LSF_1_LE_p4(i,1:hws*timestep)));
end

%% Calculate the mean & st. dev. FFT of each treatment

disp('Calculate mean & stdev of each treatment')

y_sim_fa_LSF_1_LE_0_mean = zeros(1,hws*timestep);
y_sim_fa_LSF_1_LE_p2_mean = zeros(1,hws*timestep);
y_sim_fa_LSF_1_LE_p4_mean = zeros(1,hws*timestep);

y_sim_fa_LSF_1_LE_0_std = zeros(1,hws*timestep);
y_sim_fa_LSF_1_LE_p2_std = zeros(1,hws*timestep);
y_sim_fa_LSF_1_LE_p4_std = zeros(1,hws*timestep);

te_sim_fa_LSF_1_LE_0_mean = zeros(1,hws*timestep);
te_sim_fa_LSF_1_LE_p2_mean = zeros(1,hws*timestep);
te_sim_fa_LSF_1_LE_p4_mean = zeros(1,hws*timestep);

te_sim_fa_LSF_1_LE_0_std = zeros(1,hws*timestep);
te_sim_fa_LSF_1_LE_p2_std = zeros(1,hws*timestep);
te_sim_fa_LSF_1_LE_p4_std = zeros(1,hws*timestep);

beta_sim_fa_LSF_1_LE_0_mean = zeros(1,hws*timestep);
beta_sim_fa_LSF_1_LE_p2_mean = zeros(1,hws*timestep);
beta_sim_fa_LSF_1_LE_p4_mean = zeros(1,hws*timestep);

beta_sim_fa_LSF_1_LE_0_std = zeros(1,hws*timestep);
beta_sim_fa_LSF_1_LE_p2_std = zeros(1,hws*timestep);
beta_sim_fa_LSF_1_LE_p4_std = zeros(1,hws*timestep);

for i = 1:hws*timestep
    y_sim_fa_LSF_1_LE_0_mean(1,i) = nanmean(y_sim_fa_LSF_1_LE_0(:,i));
    y_sim_fa_LSF_1_LE_p2_mean(1,i) = nanmean(y_sim_fa_LSF_1_LE_p2(:,i));
    y_sim_fa_LSF_1_LE_p4_mean(1,i) = nanmean(y_sim_fa_LSF_1_LE_p4(:,i));
    
    y_sim_fa_LSF_1_LE_0_std(1,i) = nanstd(y_sim_fa_LSF_1_LE_0(:,i));
    y_sim_fa_LSF_1_LE_p2_std(1,i) = nanstd(y_sim_fa_LSF_1_LE_p2(:,i));
    y_sim_fa_LSF_1_LE_p4_std(1,i) = nanstd(y_sim_fa_LSF_1_LE_p4(:,i));

    te_sim_fa_LSF_1_LE_0_mean(1,i) = nanmean(te_sim_fa_LSF_1_LE_0(:,i));
    te_sim_fa_LSF_1_LE_p2_mean(1,i) = nanmean(te_sim_fa_LSF_1_LE_p2(:,i));
    te_sim_fa_LSF_1_LE_p4_mean(1,i) = nanmean(te_sim_fa_LSF_1_LE_p4(:,i));
     
    te_sim_fa_LSF_1_LE_0_std(1,i) = nanstd(te_sim_fa_LSF_1_LE_0(:,i));
    te_sim_fa_LSF_1_LE_p2_std(1,i) = nanstd(te_sim_fa_LSF_1_LE_p2(:,i));
    te_sim_fa_LSF_1_LE_p4_std(1,i) = nanstd(te_sim_fa_LSF_1_LE_p4(:,i));
    
    beta_sim_fa_LSF_1_LE_0_mean(1,i) = nanmean(beta_sim_fa_LSF_1_LE_0(:,i));
    beta_sim_fa_LSF_1_LE_p2_mean(1,i) = nanmean(beta_sim_fa_LSF_1_LE_p2(:,i));
    beta_sim_fa_LSF_1_LE_p4_mean(1,i) = nanmean(beta_sim_fa_LSF_1_LE_p4(:,i));
     
    beta_sim_fa_LSF_1_LE_0_std(1,i) = nanstd(beta_sim_fa_LSF_1_LE_0(:,i));
    beta_sim_fa_LSF_1_LE_p2_std(1,i) = nanstd(beta_sim_fa_LSF_1_LE_p2(:,i));
    beta_sim_fa_LSF_1_LE_p4_std(1,i) = nanstd(beta_sim_fa_LSF_1_LE_p4(:,i));
 end

%% Calculate Gain and Phase of EACH signal
disp('Calculate the Gain and Phase of EACH signal')

Gain_ysim_fa_LSF_1_LE_0 = zeros(NumOf_fa_LSF_1_LE_0_Files,numel(prime_f));
Gain_ysim_fa_LSF_1_LE_p2 = zeros(NumOf_fa_LSF_1_LE_p2_Files,numel(prime_f));
Gain_ysim_fa_LSF_1_LE_p4 = zeros(NumOf_fa_LSF_1_LE_p4_Files,numel(prime_f));

Phase_ysim_fa_LSF_1_LE_0 = zeros(NumOf_fa_LSF_1_LE_0_Files,numel(prime_f));
Phase_ysim_fa_LSF_1_LE_p2 = zeros(NumOf_fa_LSF_1_LE_p2_Files,numel(prime_f));
Phase_ysim_fa_LSF_1_LE_p4 = zeros(NumOf_fa_LSF_1_LE_p4_Files,numel(prime_f));

Gain_te_sim_fa_LSF_1_LE_0 = zeros(NumOf_fa_LSF_1_LE_0_Files,numel(prime_f));
Gain_te_sim_fa_LSF_1_LE_p2 = zeros(NumOf_fa_LSF_1_LE_p2_Files,numel(prime_f));
Gain_te_sim_fa_LSF_1_LE_p4 = zeros(NumOf_fa_LSF_1_LE_p4_Files,numel(prime_f));

Phase_te_sim_fa_LSF_1_LE_0 = zeros(NumOf_fa_LSF_1_LE_0_Files,numel(prime_f));
Phase_te_sim_fa_LSF_1_LE_p2 = zeros(NumOf_fa_LSF_1_LE_p2_Files,numel(prime_f));
Phase_te_sim_fa_LSF_1_LE_p4 = zeros(NumOf_fa_LSF_1_LE_p4_Files,numel(prime_f));

Gain_beta_sim_fa_LSF_1_LE_0 = zeros(NumOf_fa_LSF_1_LE_0_Files,numel(prime_f));
Gain_beta_sim_fa_LSF_1_LE_p2 = zeros(NumOf_fa_LSF_1_LE_p2_Files,numel(prime_f));
Gain_beta_sim_fa_LSF_1_LE_p4 = zeros(NumOf_fa_LSF_1_LE_p4_Files,numel(prime_f));

Phase_beta_sim_fa_LSF_1_LE_0 = zeros(NumOf_fa_LSF_1_LE_0_Files,numel(prime_f));
Phase_beta_sim_fa_LSF_1_LE_p2 = zeros(NumOf_fa_LSF_1_LE_p2_Files,numel(prime_f));
Phase_beta_sim_fa_LSF_1_LE_p4 = zeros(NumOf_fa_LSF_1_LE_p4_Files,numel(prime_f));

disp('Location(s) of prime frequencies identified')
loci(1:numel(prime_f)) = [find(f_omega == prime_f(1)),...
    find(f_omega == prime_f(2)), find(f_omega == prime_f(3)),...
    find(f_omega == prime_f(4)), find(f_omega == prime_f(5)),...
    find(f_omega == prime_f(6)), find(f_omega == prime_f(7)),...
    find(f_omega == prime_f(8)), find(f_omega == prime_f(9)),...
    find(f_omega == prime_f(10)), find(f_omega == prime_f(11))];

%Gain & phase of input
Gain_ytheory = abs(y_theory(loci))./abs(y_theory(loci));
Phase_ytheory = atan(imag(y_theory(loci))./real(y_theory(loci)));

% Gain_te_theory = ones(1,numel(prime_f)); %abs(te_theory(loci))./abs(te_theory(loci));
Phase_te_theory = zeros(1,numel(prime_f)); %atan(imag(te_theory(loci))./real(te_theory(loci)));

%Gain and phase of output
%fa - LSF_1_LE_0
for i = 1:NumOf_fa_LSF_1_LE_0_Files
    for k = 1:numel(loci)
        Gain_ysim_fa_LSF_1_LE_0(i,k) = ...
            abs(y_sim_fa_LSF_1_LE_0(i,loci(k)))./abs(y_theory(loci(k)));
        Phase_ysim_fa_LSF_1_LE_0(i,k) = ...
            atan(imag(y_sim_fa_LSF_1_LE_0(i,loci(k)))./real(y_sim_fa_LSF_1_LE_0(i,loci(k))));
        Gain_te_sim_fa_LSF_1_LE_0(i,k) = ...
            abs(te_sim_fa_LSF_1_LE_0(i,loci(k)))./abs(y_theory(loci(k)));
        Phase_te_sim_fa_LSF_1_LE_0(i,k) = ...
            atan(imag(te_sim_fa_LSF_1_LE_0(i,loci(k)))./real(te_sim_fa_LSF_1_LE_0(i,loci(k))));
        Gain_beta_sim_fa_LSF_1_LE_0(i,k) = ...
            abs(beta_sim_fa_LSF_1_LE_0(i,loci(k)))./abs(y_theory(loci(k)));
        Phase_beta_sim_fa_LSF_1_LE_0(i,k) = ...
            atan(imag(beta_sim_fa_LSF_1_LE_0(i,loci(k)))./real(beta_sim_fa_LSF_1_LE_0(i,loci(k))));
    end
end

%fa - LSF_1_LE_p2
for i = 1:NumOf_fa_LSF_1_LE_p2_Files
    for k = 1:numel(loci)
        Gain_ysim_fa_LSF_1_LE_p2(i,k) = ...
            abs(y_sim_fa_LSF_1_LE_p2(i,loci(k)))./abs(y_theory(loci(k)));
        Phase_ysim_fa_LSF_1_LE_p2(i,k) = ...
            atan(imag(y_sim_fa_LSF_1_LE_p2(i,loci(k)))./real(y_sim_fa_LSF_1_LE_p2(i,loci(k))));
        Gain_te_sim_fa_LSF_1_LE_p2(i,k) = ...
            abs(te_sim_fa_LSF_1_LE_p2(i,loci(k)))./abs(y_theory(loci(k)));
        Phase_te_sim_fa_LSF_1_LE_p2(i,k) = ...
            atan(imag(te_sim_fa_LSF_1_LE_p2(i,loci(k)))./real(te_sim_fa_LSF_1_LE_p2(i,loci(k))));
        Gain_beta_sim_fa_LSF_1_LE_p2(i,k) = ...
            abs(beta_sim_fa_LSF_1_LE_p2(i,loci(k)))./abs(y_theory(loci(k)));
        Phase_beta_sim_fa_LSF_1_LE_p2(i,k) = ...
            atan(imag(beta_sim_fa_LSF_1_LE_p2(i,loci(k)))./real(beta_sim_fa_LSF_1_LE_p2(i,loci(k))));
    end
end

%fa - LSF_1_LE_p4
for i = 1:NumOf_fa_LSF_1_LE_p4_Files
    for k = 1:numel(loci)
        Gain_ysim_fa_LSF_1_LE_p4(i,k) = ...
            abs(y_sim_fa_LSF_1_LE_p4(i,loci(k)))./abs(y_theory(loci(k)));
        Phase_ysim_fa_LSF_1_LE_p4(i,k) = ...
            atan(imag(y_sim_fa_LSF_1_LE_p4(i,loci(k)))./real(y_sim_fa_LSF_1_LE_p4(i,loci(k))));
        Gain_te_sim_fa_LSF_1_LE_p4(i,k) = ...
            abs(te_sim_fa_LSF_1_LE_p4(i,loci(k)))./abs(y_theory(loci(k)));
        Phase_te_sim_fa_LSF_1_LE_p4(i,k) = ...
            atan(imag(te_sim_fa_LSF_1_LE_p4(i,loci(k)))./real(te_sim_fa_LSF_1_LE_p4(i,loci(k))));
        Gain_beta_sim_fa_LSF_1_LE_p4(i,k) = ...
            abs(beta_sim_fa_LSF_1_LE_p4(i,loci(k)))./abs(y_theory(loci(k)));
        Phase_beta_sim_fa_LSF_1_LE_p4(i,k) = ...
            atan(imag(beta_sim_fa_LSF_1_LE_p4(i,loci(k)))./real(beta_sim_fa_LSF_1_LE_p4(i,loci(k))));
    end
end

%% Calculate the mean & st. dev. of each treatment's Gain and Phase

disp('Calculate mean & stdev of Phase and Gain of each treatment')

%Gain and phase of y motion
Gain_ysim_fa_LSF_1_LE_0_mean = zeros(1,numel(prime_f));
Gain_ysim_fa_LSF_1_LE_p2_mean = zeros(1,numel(prime_f));
Gain_ysim_fa_LSF_1_LE_p4_mean = zeros(1,numel(prime_f));

Gain_ysim_fa_LSF_1_LE_0_std = zeros(1,numel(prime_f));
Gain_ysim_fa_LSF_1_LE_p2_std = zeros(1,numel(prime_f));
Gain_ysim_fa_LSF_1_LE_p4_std = zeros(1,numel(prime_f));

Phase_ysim_fa_LSF_1_LE_0_mean = zeros(1,numel(prime_f));
Phase_ysim_fa_LSF_1_LE_p2_mean = zeros(1,numel(prime_f));
Phase_ysim_fa_LSF_1_LE_p4_mean = zeros(1,numel(prime_f));

Phase_ysim_fa_LSF_1_LE_0_std = zeros(1,numel(prime_f));
Phase_ysim_fa_LSF_1_LE_p2_std = zeros(1,numel(prime_f));
Phase_ysim_fa_LSF_1_LE_p4_std = zeros(1,numel(prime_f));

%Gain and phase of tracking error
Gain_te_sim_fa_LSF_1_LE_0_mean = zeros(1,numel(prime_f));
Gain_te_sim_fa_LSF_1_LE_p2_mean = zeros(1,numel(prime_f));
Gain_te_sim_fa_LSF_1_LE_p4_mean = zeros(1,numel(prime_f));

Gain_te_sim_fa_LSF_1_LE_0_std = zeros(1,numel(prime_f));
Gain_te_sim_fa_LSF_1_LE_p2_std = zeros(1,numel(prime_f));
Gain_te_sim_fa_LSF_1_LE_p4_std = zeros(1,numel(prime_f));

Phase_te_sim_fa_LSF_1_LE_0_mean = zeros(1,numel(prime_f));
Phase_te_sim_fa_LSF_1_LE_p2_mean = zeros(1,numel(prime_f));
Phase_te_sim_fa_LSF_1_LE_p4_mean = zeros(1,numel(prime_f));

Phase_te_sim_fa_LSF_1_LE_0_std = zeros(1,numel(prime_f));
Phase_te_sim_fa_LSF_1_LE_p2_std = zeros(1,numel(prime_f));
Phase_te_sim_fa_LSF_1_LE_p4_std = zeros(1,numel(prime_f));

%Gain and phase of beta
Gain_beta_sim_fa_LSF_1_LE_0_mean = zeros(1,numel(prime_f));
Gain_beta_sim_fa_LSF_1_LE_p2_mean = zeros(1,numel(prime_f));
Gain_beta_sim_fa_LSF_1_LE_p4_mean = zeros(1,numel(prime_f));

Gain_beta_sim_fa_LSF_1_LE_0_std = zeros(1,numel(prime_f));
Gain_beta_sim_fa_LSF_1_LE_p2_std = zeros(1,numel(prime_f));
Gain_beta_sim_fa_LSF_1_LE_p4_std = zeros(1,numel(prime_f));

Phase_beta_sim_fa_LSF_1_LE_0_mean = zeros(1,numel(prime_f));
Phase_beta_sim_fa_LSF_1_LE_p2_mean = zeros(1,numel(prime_f));
Phase_beta_sim_fa_LSF_1_LE_p4_mean = zeros(1,numel(prime_f));

Phase_beta_sim_fa_LSF_1_LE_0_std = zeros(1,numel(prime_f));
Phase_beta_sim_fa_LSF_1_LE_p2_std = zeros(1,numel(prime_f));
Phase_beta_sim_fa_LSF_1_LE_p4_std = zeros(1,numel(prime_f));

for i = 1:numel(prime_f)
    %y-motion
    Gain_ysim_fa_LSF_1_LE_0_mean(1,i) = nanmean(Gain_ysim_fa_LSF_1_LE_0(:,i));
    Gain_ysim_fa_LSF_1_LE_p2_mean(1,i) = nanmean(Gain_ysim_fa_LSF_1_LE_p2(:,i));
    Gain_ysim_fa_LSF_1_LE_p4_mean(1,i) = nanmean(Gain_ysim_fa_LSF_1_LE_p4(:,i));
    
    Gain_ysim_fa_LSF_1_LE_0_std(1,i) = nanstd(Gain_ysim_fa_LSF_1_LE_0(:,i));
    Gain_ysim_fa_LSF_1_LE_p2_std(1,i) = nanstd(Gain_ysim_fa_LSF_1_LE_p2(:,i));
    Gain_ysim_fa_LSF_1_LE_p4_std(1,i) = nanstd(Gain_ysim_fa_LSF_1_LE_p4(:,i));
    
    Phase_ysim_fa_LSF_1_LE_0_mean(1,i) = nanmean(Phase_ysim_fa_LSF_1_LE_0(:,i));
    Phase_ysim_fa_LSF_1_LE_p2_mean(1,i) = nanmean(Phase_ysim_fa_LSF_1_LE_p2(:,i));
    Phase_ysim_fa_LSF_1_LE_p4_mean(1,i) = nanmean(Phase_ysim_fa_LSF_1_LE_p4(:,i));
    
    Phase_ysim_fa_LSF_1_LE_0_std(1,i) = nanstd(Phase_ysim_fa_LSF_1_LE_0(:,i));
    Phase_ysim_fa_LSF_1_LE_p2_std(1,i) = nanstd(Phase_ysim_fa_LSF_1_LE_p2(:,i));
    Phase_ysim_fa_LSF_1_LE_p4_std(1,i) = nanstd(Phase_ysim_fa_LSF_1_LE_p4(:,i));
    
    %Tracking error
    Gain_te_sim_fa_LSF_1_LE_0_mean(1,i) = nanmean(Gain_te_sim_fa_LSF_1_LE_0(:,i));
    Gain_te_sim_fa_LSF_1_LE_p2_mean(1,i) = nanmean(Gain_te_sim_fa_LSF_1_LE_p2(:,i));
    Gain_te_sim_fa_LSF_1_LE_p4_mean(1,i) = nanmean(Gain_te_sim_fa_LSF_1_LE_p4(:,i));
    
    Gain_te_sim_fa_LSF_1_LE_0_std(1,i) = nanstd(Gain_te_sim_fa_LSF_1_LE_0(:,i));
    Gain_te_sim_fa_LSF_1_LE_p2_std(1,i) = nanstd(Gain_te_sim_fa_LSF_1_LE_p2(:,i));
    Gain_te_sim_fa_LSF_1_LE_p4_std(1,i) = nanstd(Gain_te_sim_fa_LSF_1_LE_p4(:,i));
   
    Phase_te_sim_fa_LSF_1_LE_0_mean(1,i) = nanmean(Phase_te_sim_fa_LSF_1_LE_0(:,i));
    Phase_te_sim_fa_LSF_1_LE_p2_mean(1,i) = nanmean(Phase_te_sim_fa_LSF_1_LE_p2(:,i));
    Phase_te_sim_fa_LSF_1_LE_p4_mean(1,i) = nanmean(Phase_te_sim_fa_LSF_1_LE_p4(:,i));
   
    Phase_te_sim_fa_LSF_1_LE_0_std(1,i) = nanstd(Phase_te_sim_fa_LSF_1_LE_0(:,i));
    Phase_te_sim_fa_LSF_1_LE_p2_std(1,i) = nanstd(Phase_te_sim_fa_LSF_1_LE_p2(:,i));
    Phase_te_sim_fa_LSF_1_LE_p4_std(1,i) = nanstd(Phase_te_sim_fa_LSF_1_LE_p4(:,i));
    
    %Beta
    Gain_beta_sim_fa_LSF_1_LE_0_mean(1,i) = nanmean(Gain_beta_sim_fa_LSF_1_LE_0(:,i));
    Gain_beta_sim_fa_LSF_1_LE_p2_mean(1,i) = nanmean(Gain_beta_sim_fa_LSF_1_LE_p2(:,i));
    Gain_beta_sim_fa_LSF_1_LE_p4_mean(1,i) = nanmean(Gain_beta_sim_fa_LSF_1_LE_p4(:,i));
    
    Gain_beta_sim_fa_LSF_1_LE_0_std(1,i) = nanstd(Gain_beta_sim_fa_LSF_1_LE_0(:,i));
    Gain_beta_sim_fa_LSF_1_LE_p2_std(1,i) = nanstd(Gain_beta_sim_fa_LSF_1_LE_p2(:,i));
    Gain_beta_sim_fa_LSF_1_LE_p4_std(1,i) = nanstd(Gain_beta_sim_fa_LSF_1_LE_p4(:,i));
   
    Phase_beta_sim_fa_LSF_1_LE_0_mean(1,i) = nanmean(Phase_beta_sim_fa_LSF_1_LE_0(:,i));
    Phase_beta_sim_fa_LSF_1_LE_p2_mean(1,i) = nanmean(Phase_beta_sim_fa_LSF_1_LE_p2(:,i));
    Phase_beta_sim_fa_LSF_1_LE_p4_mean(1,i) = nanmean(Phase_beta_sim_fa_LSF_1_LE_p4(:,i));
   
    Phase_beta_sim_fa_LSF_1_LE_0_std(1,i) = nanstd(Phase_beta_sim_fa_LSF_1_LE_0(:,i));
    Phase_beta_sim_fa_LSF_1_LE_p2_std(1,i) = nanstd(Phase_beta_sim_fa_LSF_1_LE_p2(:,i));
    Phase_beta_sim_fa_LSF_1_LE_p4_std(1,i) = nanstd(Phase_beta_sim_fa_LSF_1_LE_p4(:,i));
end

%% Find the min & max values for...

%fft y_sim
fft_y_Max_vec = [max(abs(y_sim_fa_LSF_1_LE_0_mean(1:maxrange))/ampscale),...
    max(abs(y_sim_fa_LSF_1_LE_p2_mean(1:maxrange))/ampscale),...
    max(abs(y_sim_fa_LSF_1_LE_p4_mean(1:maxrange))/ampscale)];
fft_y_max = max(fft_y_Max_vec);

%fft te_sim
fft_te_Max_vec = [max(abs(te_sim_fa_LSF_1_LE_0_mean(1:maxrange))/ampscale),...
    max(abs(te_sim_fa_LSF_1_LE_p2_mean(1:maxrange))/ampscale),...
    max(abs(te_sim_fa_LSF_1_LE_p4_mean(1:maxrange))/ampscale)];
fft_te_max = max(fft_te_Max_vec);

%fft beta_sim
fft_beta_Max_vec = [max(abs(beta_sim_fa_LSF_1_LE_0_mean(1:maxrange))/ampscale),...
    max(abs(beta_sim_fa_LSF_1_LE_p2_mean(1:maxrange))/ampscale),...
    max(abs(beta_sim_fa_LSF_1_LE_p4_mean(1:maxrange))/ampscale)];
fft_beta_max = max(fft_beta_Max_vec);

%Gain of y-motion
gain_y_Max_vec = [max(max(Gain_ysim_fa_LSF_1_LE_0_mean)),...
    max(max(Gain_ysim_fa_LSF_1_LE_p2_mean)),...
    max(max(Gain_ysim_fa_LSF_1_LE_p4_mean))];
gain_y_max = max(gain_y_Max_vec);

gain_y_Min_vec = [min(min(Gain_ysim_fa_LSF_1_LE_0_mean)),...
    min(min(Gain_ysim_fa_LSF_1_LE_p2_mean)),...
    min(min(Gain_ysim_fa_LSF_1_LE_p4_mean))];
gain_y_min = min(gain_y_Min_vec);

%Phase of y-motion
phase_y_Max_vec = [max(max(Phase_ysim_fa_LSF_1_LE_0_mean-Phase_ytheory)),...
    max(max(Phase_ysim_fa_LSF_1_LE_p2_mean-Phase_ytheory)),...
    max(max(Phase_ysim_fa_LSF_1_LE_p4_mean-Phase_ytheory))];
phase_y_max = max(phase_y_Max_vec);

phase_y_Min_vec = [min(min(Phase_ysim_fa_LSF_1_LE_0_mean-Phase_ytheory)),...
    min(min(Phase_ysim_fa_LSF_1_LE_p2_mean-Phase_ytheory)),...
    min(min(Phase_ysim_fa_LSF_1_LE_p4_mean-Phase_ytheory))];
phase_y_min = min(phase_y_Min_vec);

%Gain of tracking error
gain_te_Max_vec = [max(max(Gain_te_sim_fa_LSF_1_LE_0_mean)),...
    max(max(Gain_te_sim_fa_LSF_1_LE_p2_mean)),...
    max(max(Gain_te_sim_fa_LSF_1_LE_p4_mean))];
gain_te_max = max(gain_te_Max_vec);

gain_te_Min_vec = [min(min(Gain_te_sim_fa_LSF_1_LE_0_mean)),...
    min(min(Gain_te_sim_fa_LSF_1_LE_p2_mean)),...
    min(min(Gain_te_sim_fa_LSF_1_LE_p4_mean))];
gain_te_min = min(gain_te_Min_vec);

%Phase of tracking error
phase_te_Max_vec = [max(max(Phase_te_sim_fa_LSF_1_LE_0_mean)),...
    max(max(Phase_te_sim_fa_LSF_1_LE_p2_mean)),...
    max(max(Phase_te_sim_fa_LSF_1_LE_p4_mean))];
phase_te_max = max(phase_te_Max_vec);

phase_te_Min_vec = [min(min(Phase_te_sim_fa_LSF_1_LE_0_mean)),...
    min(min(Phase_te_sim_fa_LSF_1_LE_p2_mean)),...
    min(min(Phase_te_sim_fa_LSF_1_LE_p4_mean))];
phase_te_min = min(phase_te_Min_vec);

%Gain of beta
gain_beta_Max_vec = [max(max(Gain_beta_sim_fa_LSF_1_LE_0_mean)),...
    max(max(Gain_beta_sim_fa_LSF_1_LE_p2_mean)),...
    max(max(Gain_beta_sim_fa_LSF_1_LE_p4_mean))];
gain_beta_max = max(gain_beta_Max_vec);

gain_beta_Min_vec = [min(min(Gain_beta_sim_fa_LSF_1_LE_0_mean)),...
    min(min(Gain_beta_sim_fa_LSF_1_LE_p2_mean)),...
    min(min(Gain_beta_sim_fa_LSF_1_LE_p4_mean))];
gain_beta_min = min(gain_beta_Min_vec);

%Phase of beta
phase_beta_Max_vec = [max(max(Phase_beta_sim_fa_LSF_1_LE_0_mean)),...
    max(max(Phase_beta_sim_fa_LSF_1_LE_p2_mean)),...
    max(max(Phase_beta_sim_fa_LSF_1_LE_p4_mean))];
phase_beta_max = max(phase_beta_Max_vec);

phase_beta_Min_vec = [min(min(Phase_beta_sim_fa_LSF_1_LE_0_mean)),...
    min(min(Phase_beta_sim_fa_LSF_1_LE_p2_mean)),...
    min(min(Phase_beta_sim_fa_LSF_1_LE_p4_mean))];
phase_beta_min = min(phase_beta_Min_vec);

%% Plot the FFT of y
fig7 = figure(7);
%Y motion fft
subplot(3,3,1)
plot(f_omega(1:maxrange),abs(y_sim_fa_LSF_1_LE_0_mean(1:maxrange))/ampscale,'k','LineWidth',2);
hold on;
plot(prime_f, prime_a, 'r*')
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0'})
ylabel('fft(y)')
axis([0, 20, -0.01, fft_y_max])

subplot(3,3,2)
plot(f_omega(1:maxrange),abs(y_sim_fa_LSF_1_LE_p2_mean(1:maxrange))/ampscale,'k','LineWidth',2);
hold on;
plot(prime_f, prime_a, 'r*')
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0.2'})
axis([0, 20, -0.01, fft_y_max])

subplot(3,3,3)
plot(f_omega(1:maxrange),abs(y_sim_fa_LSF_1_LE_p4_mean(1:maxrange))/ampscale,'k','LineWidth',2);
hold on;
plot(prime_f, prime_a, 'r*')
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0.4'})
axis([0, 20, -0.01, fft_y_max])

%Tracking error fft
subplot(3,3,4)
plot(f_omega(1:maxrange),abs(te_sim_fa_LSF_1_LE_0_mean(1:maxrange))/ampscale,'k','LineWidth',2);
hold on;
plot(prime_f, zeros(1,numel(prime_a)), 'r*')
ylabel('fft(tracking error)')
axis([0, 20, -0.01, fft_te_max])

subplot(3,3,5)
plot(f_omega(1:maxrange),abs(te_sim_fa_LSF_1_LE_p2_mean(1:maxrange))/ampscale,'k','LineWidth',2);
hold on;
plot(prime_f, zeros(1,numel(prime_a)), 'r*')
axis([0, 20, -0.01, fft_te_max])

subplot(3,3,6)
plot(f_omega(1:maxrange),abs(te_sim_fa_LSF_1_LE_p4_mean(1:maxrange))/ampscale,'k','LineWidth',2);
hold on;
plot(prime_f, zeros(1,numel(prime_a)), 'r*')
axis([0, 20, -0.01, fft_te_max])

%Beta fft
subplot(3,3,7)
plot(f_omega(1:maxrange),abs(beta_sim_fa_LSF_1_LE_0_mean(1:maxrange))/ampscale,'k','LineWidth',2);
hold on;
plot(prime_f, zeros(1,numel(prime_a)), 'r*')
xlabel({'f (Hz)';['(n=',num2str(NumOf_fa_LSF_1_LE_0_Files),')']});
ylabel('fft(beta)')
axis([0, 20, -0.01, fft_beta_max])

subplot(3,3,8)
plot(f_omega(1:maxrange),abs(beta_sim_fa_LSF_1_LE_p2_mean(1:maxrange))/ampscale,'k','LineWidth',2);
hold on;
plot(prime_f, zeros(1,numel(prime_a)), 'r*')
xlabel({'f (Hz)';['(n=',num2str(NumOf_fa_LSF_1_LE_p2_Files),')']});
axis([0, 20, -0.01, fft_beta_max])

subplot(3,3,9)
plot(f_omega(1:maxrange),abs(beta_sim_fa_LSF_1_LE_p4_mean(1:maxrange))/ampscale,'k','LineWidth',2);
hold on;
plot(prime_f, zeros(1,numel(prime_a)), 'r*')
xlabel({'f (Hz)';['(n=',num2str(NumOf_fa_LSF_1_LE_p4_Files),')']});
axis([0, 20, -0.01, fft_beta_max])

set(fig7, 'Position', [120,120,850,525])

saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig7_fa_fft_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig7_fa_fft_LSF_1_v12h.jpg');

disp('Figure 7 done')

%% Gain and Phase of the y-motion
fig8 = figure(8);

%Gain y-motion
subplot(2,3,1)
plot(prime_f, ones, 'ko')
hold on;
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0'})
plot(f_omega(loci),Gain_ysim_fa_LSF_1_LE_0_mean,'k','LineWidth',2)
ylabel('Gain(y)')
axis([0, 20, gain_y_min, gain_y_max])

subplot(2,3,2)
plot(prime_f, ones, 'ko')
hold on;
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0.2'})
plot(f_omega(loci),Gain_ysim_fa_LSF_1_LE_p2_mean,'k','LineWidth',2)
axis([0, 20, gain_y_min, gain_y_max])

subplot(2,3,3)
plot(prime_f, ones, 'ko')
hold on;
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0.4'})
plot(f_omega(loci),Gain_ysim_fa_LSF_1_LE_p4_mean,'k','LineWidth',2)
axis([0, 20, gain_y_min, gain_y_max])

%Phase y-motion
subplot(2,3,4)
plot(f_omega(loci),180/pi*(Phase_ysim_fa_LSF_1_LE_0_mean...
    -Phase_ytheory),'k','LineWidth',2)
xlabel({'f (Hz)';['(n=',num2str(NumOf_fa_LSF_1_LE_0_Files),')']});
ylabel({'y-motion';'Phase difference (deg)'; '(output-input)'})
axis([0, 20, (180/pi)*phase_y_min, (180/pi)*phase_y_max])

subplot(2,3,5)
plot(f_omega(loci),180/pi*(Phase_ysim_fa_LSF_1_LE_p2_mean...
    -Phase_ytheory),'k','LineWidth',2)
xlabel({'f (Hz)';['(n=',num2str(NumOf_fa_LSF_1_LE_p2_Files),')']});
axis([0, 20, (180/pi)*phase_y_min, (180/pi)*phase_y_max])

subplot(2,3,6)
plot(f_omega(loci),180/pi*(Phase_ysim_fa_LSF_1_LE_p4_mean...
    -Phase_ytheory),'k','LineWidth',2)
xlabel({'f (Hz)';['(n=',num2str(NumOf_fa_LSF_1_LE_p4_Files),')']});
axis([0, 20, (180/pi)*phase_y_min, (180/pi)*phase_y_max])

set(fig8, 'Position', [120,120,850,525])
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig8_fa_bode_y_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig8_fa_bode_y_LSF_1_v12h.jpg');

disp('Figure 8 done')

%% Gain and Phase of the tracking error
fig9 = figure(9);

%Gain of tracking error
subplot(2,3,1)
plot(prime_f, ones, 'ko')
hold on;
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0'})
plot(f_omega(loci),Gain_te_sim_fa_LSF_1_LE_0_mean,'k','LineWidth',2)
ylabel('Gain(tracking error)')
axis([0, 20, gain_te_min, gain_te_max])

subplot(2,3,2)
plot(prime_f, ones, 'ko')
hold on;
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0.2'})
plot(f_omega(loci),Gain_te_sim_fa_LSF_1_LE_p2_mean,'k','LineWidth',2)
axis([0, 20, gain_te_min, gain_te_max])

subplot(2,3,3)
plot(prime_f, ones, 'ko')
hold on;
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0.4'})
plot(f_omega(loci),Gain_te_sim_fa_LSF_1_LE_p4_mean,'k','LineWidth',2)
axis([0, 20, gain_te_min, gain_te_max])

%Phase of tracking error
subplot(2,3,4)
plot(f_omega(loci),180/pi*(Phase_te_sim_fa_LSF_1_LE_0_mean),'k','LineWidth',2)
xlabel({'f (Hz)';['(n=',num2str(NumOf_fa_LSF_1_LE_0_Files),')']});
ylabel({'Tracking error';'Phase difference (deg)'; '(output-input)'})
axis([0, 20, (180/pi)*phase_te_min, (180/pi)*phase_te_max])

subplot(2,3,5)
plot(f_omega(loci),180/pi*(Phase_te_sim_fa_LSF_1_LE_p2_mean),'k','LineWidth',2)
xlabel({'f (Hz)';['(n=',num2str(NumOf_fa_LSF_1_LE_p2_Files),')']});
axis([0, 20, (180/pi)*phase_te_min, (180/pi)*phase_te_max])

subplot(2,3,6)
plot(f_omega(loci),180/pi*(Phase_te_sim_fa_LSF_1_LE_p4_mean),'k','LineWidth',2)
xlabel({'f (Hz)';['(n=',num2str(NumOf_fa_LSF_1_LE_p4_Files),')']});
axis([0, 20, (180/pi)*phase_te_min, (180/pi)*phase_te_max])

set(fig9, 'Position', [120,120,850,525])
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig9_fa_bode_te_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig9_fa_bode_te_LSF_1_v12h.jpg');

disp('Figure 9 done')

%% Gain and Phase of beta
fig10 = figure(10);

%Gain of beta
subplot(2,3,1)
plot(prime_f, ones, 'ko')
hold on;
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0'})
plot(f_omega(loci),Gain_beta_sim_fa_LSF_1_LE_0_mean,'k','LineWidth',2)
ylabel('Gain(beta)')
axis([0, 20, gain_beta_min, gain_beta_max])

subplot(2,3,2)
plot(prime_f, ones, 'ko')
hold on;
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0.2'})
plot(f_omega(loci),Gain_beta_sim_fa_LSF_1_LE_p2_mean,'k','LineWidth',2)
axis([0, 20, gain_beta_min, gain_beta_max])

subplot(2,3,3)
plot(prime_f, ones, 'ko')
hold on;
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0.4'})
plot(f_omega(loci),Gain_beta_sim_fa_LSF_1_LE_p4_mean,'k','LineWidth',2)
axis([0, 20, gain_beta_min, gain_beta_max])

%Phase of beta
subplot(2,3,4)
plot(f_omega(loci),180/pi*(Phase_beta_sim_fa_LSF_1_LE_0_mean),'k','LineWidth',2)
xlabel({'f (Hz)';['(n=',num2str(NumOf_fa_LSF_1_LE_0_Files),')']});
ylabel({'Beta';'Phase difference (deg)'; '(output-input)'})
axis([0, 20, (180/pi)*phase_beta_min, (180/pi)*phase_beta_max])

subplot(2,3,5)
plot(f_omega(loci),180/pi*(Phase_beta_sim_fa_LSF_1_LE_p2_mean),'k','LineWidth',2)
xlabel({'f (Hz)';['(n=',num2str(NumOf_fa_LSF_1_LE_p2_Files),')']});
axis([0, 20, (180/pi)*phase_beta_min, (180/pi)*phase_beta_max])

subplot(2,3,6)
plot(f_omega(loci),180/pi*(Phase_beta_sim_fa_LSF_1_LE_p4_mean),'k','LineWidth',2)
xlabel({'f (Hz)';['(n=',num2str(NumOf_fa_LSF_1_LE_p4_Files),')']});
axis([0, 20, (180/pi)*phase_beta_min, (180/pi)*phase_beta_max])

set(fig10, 'Position', [120,120,850,525])
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig10_fa_bode_beta_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig10_fa_bode_beta_LSF_1_v12h.jpg');

disp('Figure 10 done')

%% Gain and Phase of the y-motion
fig8a = figure(11);

%Gain y-motion
subplot(2,1,1)
plot(prime_f, ones, 'ko')
hold on;
title({'fa';['LSF: ',num2str(LSF_1)];'LE: 0'})
plot(f_omega(loci),Gain_ysim_fa_LSF_1_LE_0_mean,'k','LineWidth',2)
ylabel('Gain(y)')
axis([0, 20, gain_y_min, gain_y_max])

%Phase y-motion
subplot(2,1,2)
plot(f_omega(loci),180/pi*(Phase_ysim_fa_LSF_1_LE_0_mean...
    -Phase_ytheory),'k','LineWidth',2)
xlabel({'f (Hz)';['(n=',num2str(NumOf_fa_LSF_1_LE_0_Files),')']});
ylabel({'y-motion';'Phase difference (deg)'; '(output-input)'})
axis([0, 20, (180/pi)*phase_y_min, (180/pi)*phase_y_max])

set(fig8, 'Position', [120,120,850,525])
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig8a_fa_bode_y_LSF_1_v12h.fig');
saveas(gcf,'../Figures_MPC/fa/LSF_1/Fig8a_fa_bode_y_LSF_1_v12h.jpg');

disp('Figure 8a done')

%% Final window
disp('Completely done, son')
disp(['Time at end: ',datestr(now)])