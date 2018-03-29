% Quick simulation - shows the fits on noisy data
% Simulates curves for tissue of interest and reference region with noise
% then shows the fits by the ERRM.
% This figure was adapted into Supplementary Figure S1 in the manuscript.

% Estimated runtime: < 1 second

%% Pre-setup
clearvars
addpath('./mfiles')
tic

%% Settings

% Noise - manuscript version uses sigmaCt=0.02, sigmaCrr=0.002
sigmaCt = 0.02;
sigmaCrr = 0.1*sigmaCt;

% Parameters for tissue of interest
kt = 0.25;
kep = 0.625;
ve = kt/kep;
vpList = [0, 0.025, 0.05, 0.10];

% Parameters for reference tissue
ktRR = 0.07;
kepRR = 0.5;
veRR = ktRR/kepRR;

% Properties of the simulation
sSize= .1; % Temporal resolution for generating simulated curves (in seconds)
tRes = 1; % Desired temporal resolution (in seconds), obtained by downsampling simulated data
t=sSize:sSize:600; % Time, in seconds
t=t'/60; % Convert time into minutes

% Generate the arterial input function using literature-based model, bolus arrival at 60s
Cp = ParkerAif(t,t(60/sSize));

%% Simulate concentration in Tissue of Interest
for i=1:length(vpList)
    Ct(:,i) = ToftsKety(Cp,[kt kep vpList(i)],t,1);
end

%% Simulate concentration in Reference Region
Crr = ToftsKety(Cp,[ktRR kepRR],t,0);
trueCrr = Crr;

%% Add noise
noiselessCt = Ct;
noiselessCrr = Crr;
noiselessT = t;
Ct = Ct + sigmaCt * randn(size(Ct));
Crr = Crr + sigmaCrr * randn(size(Crr));

%% Downsample

dFactor = tRes/sSize;
t=downsample(t,dFactor);
Ct=downsample(Ct,dFactor);
Crr=downsample(Crr,dFactor);
Cp = downsample(Cp,dFactor);

%% Fit models to simulated data

[pkRRM, fittedRRM] = LRRM(Ct,Crr,t); % RRM
[pkERRM, fittedERRM] = ERRM(Ct,Crr,t); % ERRM
[pkERRMn, fittedERRMn] = ERRM(Ct,Crr,t,true); % ERRM - using lsqnonneg
[pkCERRM, fittedCERRM] = CERRM(Ct,Crr,t); % CERRM

% Doing the Reference Tissue Method
% Results are not shown
[pkRTM, estCpRTM] = doRTM(Ct,Crr,t,[ktRR veRR],0);
[pkERTM, estCpERTM] = doRTM(Ct,Crr,t,[ktRR veRR],1);
fittedRTM = ToftsKety(estCpRTM,pkRTM,t,0);
fittedERTM = ToftsKety(estCpERTM,pkERTM,t,1);

%% Plot results

% Define colours for tissue curve plots
listColours{5} = [55,200,113]./255; % Ct vp=0
listColours{4} = [0,102,128]./255; % Ct vp=0.025
listColours{3} = [0,136,170]./255; % Ct vp=0.05
listColours{2} = [0,204,255]./255; % Ct vp=0.10
listColours{1} = [181,238,252]./255; % Crr

% Show the early part where most of interesting things happen
opt_xlim = [0.5 4.0];

% Figure might be too wide for users with smaller screen?
figure('Position',[300 300 1200 700])

% Plot Cp
subplot(2,3,1)
plot(t,Cp,'r', 'LineWidth',2)
xlabel('Time [min]')
ylabel('Concentration [mM]')
title('Blood plasma (C_p)')
legend('C_p(t) - Blood Plasma')
xlim(opt_xlim)

% Plot noiseless curves
subplot(2,3,2)
iterPlot(noiselessT,[noiselessCt noiselessCrr],listColours)
xlabel('Time [min]')
ylabel('Concentration [mM]')
title('Simulated tissue curves (C_t and C_{RR})')
legend('vp=0','vp=0.025','vp=0.05','vp=0.10','C_{RR}(t)')
xlim(opt_xlim)
ylim([-0.01 1.4])

% Plot noisy curves
subplot(2,3,3)
iterPlot(t,[Ct Crr],listColours)
xlabel('Time [min]')
ylabel('Concentration [mM]')
title('Simulated tissue curves with noise')
legend('vp=0','vp=0.025','vp=0.05','vp=0.10','C_{RR}(t)')
xlim(opt_xlim)
ylim([-0.01 1.4])

% Plot ERRM fits on noisy data
subplot(2,3,4)
iterPlot(t,fittedERRM,listColours)
iterScatter(t,cumtrapz(t,Ct),5,listColours)
xlabel('Time [min]')
ylabel('Integral of Concentration [mM*min]')
title('ERRM fit on noisy curves')
xlim(opt_xlim)
ylim([-0.01 1.8])

% Plot ERRM fits on noisy data - zoomed in
subplot(2,3,5)
iterPlot(t,fittedERRM,listColours)
iterScatter(t,cumtrapz(t,Ct),40,listColours)
xlabel('Time [min]')
ylabel('Integral of Concentration [mM*min]')
title('ERRM fit on noisy curves - zoomed in')
xlim([1 1.5])
ylim([-0.01 0.3])

% Plot residuals of ERRM fit
subplot(2,3,6)
iterPlot(t,cumtrapz(t,Ct)-fittedERRM,listColours)
xlabel('Time [min]')
ylabel('Residuals [mM*min]')
title('Residuals of ERRM fit on noisy curves')
xlim(opt_xlim)
%%
toc
disp('Finished quick simulation showing ERRM fit on noisy data')