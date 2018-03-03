% [OBSOLETE] Quick simulation
% Simulates a noiseless curve for tissue of interest and reference region
% Plasma volume component is included in tissue of interest
% with vp = 0.001 to 0.1
% Simulated data is fitted with reference region models

% Only some of the results are shown at the end.

% Estimated runtime: <10 seconds

%% Pre-setup

% Start with a clean slate
clearvars 

% Add auxiliary files
addpath('./mfiles')

% Start timing the script
tic

%% Settings

% Sigma (stdDev) for noise level in tissue of interest and reference region
% This script works well if sigmaCt=0, sigmaCrr=0
% otherwise, the y-axis limits will have to be modified
sigmaCt = 0.0;
sigmaCrr = 0.1*sigmaCt;

% Parameters for tissue of interest
kt = 0.25; % Ktrans in 1/min
kep = 0.625; % kep in 1/min
ve = kt/kep;

% Parameters for reference tissue
ktRR = 0.07; % units: 1/min
kepRR = 0.5; % units: 1/min
veRR = ktRR/kepRR;

% Range of plasma volumes for tissue of interest
vpRange = 0:0.001:0.1;

% Properties of the simulation
sSize= 0.1; % Temporal resolution for generating simulated curves (in seconds)
tRes = 1; % Desired temporal resolution (in seconds), obtained by downsampling simulated data
t= sSize:sSize:600; % Time, in seconds
t= t'/60; % Convert time from seconds to minutes

% Generate the arterial input function using literature-based model
% Bolus arrival at 60 seconds
Cp = ParkerAif(t,t(60/sSize));

%% Simulate concentration in Tissue of Interest
Ct = zeros(length(t),length(vpRange));
for i=1:length(vpRange)
    Ct(:,i) = ToftsKety(Cp,[kt kep vpRange(i)],t,1);
end
% Alternative method: simulate Ct with vp=0
% then make copies of Ct and add vp component to each copy

%% Simulate concentration in Reference Region
Crr = ToftsKety(Cp,[ktRR kepRR],t,0);

%% Add noise
% Has no effect if sigmaC* = 0
Ct = Ct + sigmaCt * randn(size(Ct));
Crr = Crr + sigmaCrr * randn(size(Crr));

%% Downsample
% Has no effect if tRes=sSize;
dFactor = tRes/sSize;
t=downsample(t,dFactor);
Ct=downsample(Ct,dFactor);
Crr=downsample(Crr,dFactor);
Cp = downsample(Cp,dFactor);

%% Fit models to simulated data

pkRRM = LRRM(Ct,Crr,t); % RRM
pkERRM = ERRM(Ct,Crr,t); % ERRM
pkERRMn = ERRM(Ct,Crr,t,true); % ERRM - using lsqnonneg
pkCERRM = CERRM(Ct,Crr,t); % CERRM
pkCERRMn = CERRM(Ct,Crr,t,[],true); % CERRM using lsqnonneg

% RTM = Reference Tissue Method - Not detailed in revised manuscript
pkRTM = doRTM(Ct,Crr,t,[ktRR veRR],0); % RTM with Tofts Model
pkERTM = doRTM(Ct,Crr,t,[ktRR veRR],1); % RTM with extended Tofts Model

%% Calculate the percent errors

errKtRRM = PercentError(pkRRM(:,1)*ktRR,kt);
errVeRRM = PercentError(pkRRM(:,2)*veRR,ve);
errKepRRM = PercentError(pkRRM(:,3),kep);

errKtERRM = PercentError(pkERRM(:,1)*ktRR,kt);
errVeERRM = PercentError(pkERRM(:,2)*veRR,ve);
errKepERRM = PercentError(pkERRM(:,3),kep);
errVpERRM = PercentError(pkERRM(:,4)*ktRR,vpRange);

errKtERRMn = PercentError(pkERRMn(:,1)*ktRR,kt);
errVeERRMn = PercentError(pkERRMn(:,2)*veRR,ve);
errKepERRMn = PercentError(pkERRMn(:,3),kep);
errVpERRMn = PercentError(pkERRMn(:,4)*ktRR,vpRange);

errKtCERRM = PercentError(pkCERRM(:,1)*ktRR,kt);
errVeCERRM = PercentError(pkCERRM(:,2)*veRR,ve);
errKepCERRM = PercentError(pkCERRM(:,3),kep);
errVpCERRM = PercentError(pkCERRM(:,4)*ktRR,vpRange);

errKtCERRMn = PercentError(pkCERRMn(:,1)*ktRR,kt);
errVeCERRMn = PercentError(pkCERRMn(:,2)*veRR,ve);
errKepCERRMn = PercentError(pkCERRMn(:,3),kep);
errVpCERRMn = PercentError(pkCERRMn(:,4)*ktRR,vpRange);

errKtRTM = PercentError(pkRTM(:,1),kt);
errVeRTM = PercentError(pkRTM(:,1)./pkRTM(:,2),ve);
errKepRTM = PercentError(pkRTM(:,2),kep);

errKtERTM = PercentError(pkERTM(:,1),kt);
errVeERTM = PercentError(pkERTM(:,1)./pkERTM(:,2),ve);
errKepERTM = PercentError(pkERTM(:,2),kep);
errVpERTM = PercentError(pkERTM(:,3),vpRange);

%% Plot results

% Specify the line colours
lineColours{1} = [255,44,127]./255; % RRM - pink/red
lineColours{2} = [255,196,68]./255; % ERRM - orange/gold
lineColours{3} = [126,47,142]./255; % CERRM - purple

% Define the line styles
lineStyles{1} = '-'; % RRM - solid
lineStyles{2} = '-'; % ERRM - solid
lineStyles{3} = '--'; % CERRM - dashed

% x-axis is plasma volume
x = vpRange;

% Plot everything as subplots in one figure
figure('Position',[300 300 700 700])

% Percent error in Ktrans
subplot(2,2,1)
iterPlot(x,[errKtRRM,errKtERRM,errKtCERRM],lineColours,lineStyles)
legend('RRM','ERRM','CERRM','Location','northwest')
xlabel('Plasma volume')
ylabel('Percent Error')
title('K^{trans}')

% Percent error in kep
subplot(2,2,2)
iterPlot(x,[errKepRRM,errKepERRM,errKepCERRM],lineColours,lineStyles)
legend('RRM','ERRM','CERRM','Location','northwest')
xlabel('Plasma volume')
title('k_{ep}')

% Percent error in ve
subplot(2,2,3)
iterPlot(x,[errVeRRM,errVeERRM,errVeCERRM],lineColours,lineStyles)
legend('RRM','ERRM','CERRM','Location','northwest')
xlabel('Plasma volume')
ylabel('Percent Error')
title('v_e')

% Percent error in vp (ERRM and CERRM only, no RRM)
subplot(2,2,4)
iterPlot(x,[errVpERRM,errVpCERRM],lineColours(2:3),lineStyles(2:3))
legend('ERRM','CERRM')
xlabel('Plasma volume')
ylabel('Percent Error')
title('v_p')
ylim([-0.05 0.2])
%%
toc
disp('Finished quick simulation')
