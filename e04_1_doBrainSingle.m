% Process the GBM data with RRM, ERRM, CERRM, and ETM

% Note, there is additional code at the end for producing figures.
% They are not produced by default because it leads to a flood of figures.
% Change setting in the next cell to show these figures.

% Estimated runtime: ~5 seconds

%% Pre-setup
clearvars
addpath('./mfiles')

% Select which dataset to process
i=1; % Options: 1, 2, 3 (called A, B, C, in manuscripts)
% Figure 4 in manuscript uses i=1

% Decide if you want a flood of figures
openFloodgates = 0; % Options: 0 (false) or 1 (true)

%% Load data
matList = {'TCGA-06-0185-1.mat','TCGA-06-0185-2.mat','TCGA-06-0881-1.mat'};
load(fullfile('./data/gbm',matList{i}));
% The loaded data contains:
% t   [70x1] - time (in minutes) at each dce frame
% Cp  [70x1] - concentration in blood plasma
% Crr [70x1] - concentration in reference region
% Ct  [70xN] - concentration in N tumour voxels
% mask[256x256x16] - mask of the tumour region (used to reconstruct maps)
% Unnecessary items also included in loaded data:
% sT - number of time points
% sX,sY,sZ - dimensions of mask

%% Calculate standard deviation of pre-contrast frames
preContrastFrames = Ct(1:6,:);
stdDev_preContrast = std(preContrastFrames(:)); % ~0.015

stdDev_perContrastRR = std(Crr(1:6));

%% Fit with ERRM and CERRM

% CERRM() provides paramaters for both ERRM (pkE) and CERRM (pkCE)
disp('Fitting ERRM and CERRM')
tic
[pkCE, ~, estKepRR_CE, pkE]=CERRM(Ct,Crr,t);
toc
%% Fit with non-extended RRM
disp('Fitting Non-extended RRM')
tic
pkR=LRRM(Ct,Crr,t);
toc

%% Fit with ExtToftsModel
disp('Fitting Extended Tofts Model')
tic
[pkT residT fitETM] = LLSQ(Ct,Cp,t,1);
toc

%% Fit the reference region with Tofts models

% Tofts model fit
pkCrr = LLSQ(Crr,Cp,t,0);
ktRR = pkCrr(1);
kepRR = pkCrr(2);
veRR = ktRR/kepRR;
% Extended Tofts model fit
pkCrr = LLSQ(Crr,Cp,t,1);
vpRR = pkCrr(3);

%% Correlations

fprintf('\n')
disp('CORRELATIONS:');

% Only use voxels where max concentration is between 0.05 and 1 mM
corrMask = max(Ct)'>0.05 & max(Ct)'<1;

% Ktrans correlation corefficients
cKtR=getCorrCoef(pkT(corrMask,1),pkR(corrMask,1));
cKtE=getCorrCoef(pkT(corrMask,1),pkE(corrMask,1));
cKtCE=getCorrCoef(pkT(corrMask,1),pkCE(corrMask,1));
disp('Correlation with Ktrans [RRM ERRM CERRM]')
disp([cKtR cKtE cKtCE])

% kep correlation coefficients
cKepR=getCorrCoef(pkT(corrMask,2),pkR(corrMask,3));
cKepE=getCorrCoef(pkT(corrMask,2),pkE(corrMask,3));
cKepCE=getCorrCoef(pkT(corrMask,2),pkCE(corrMask,3));
disp('Correlation with kep [RRM ERRM CERRM]')
disp([cKepR cKepE cKepCE])

% ve correlation coefficients - requires extra bit of filtering
% There are sometimes extreme values (outliers) in ve, and the correlation
% cofficient is sensitive to these outliers. So, we only use data between
% the 1-99 percentile as a way to filter out the outliers.
qtRange = [.01 .99];
estVe = pkT(:,1)./pkT(:,2);
qtETM = quantile(estVe(corrMask),qtRange);
qtR = quantile(pkR(corrMask,2),qtRange);
qtE = quantile(pkE(corrMask,2),qtRange);
qtCE = quantile(pkCE(corrMask,2),qtRange);
veMask = (estVe > qtETM(1) & estVe < qtETM(2)) ... 
    & (pkR(:,2) > qtR(1) & pkR(:,2) < qtR(2)) ...
    & (pkE(:,2) > qtE(1) & pkE(:,2) < qtE(2)) ...
    & (pkCE(:,2) > qtCE(1) & pkCE(:,2) < qtCE(2));
veMask = veMask & corrMask;

cVeR=getCorrCoef(estVe(veMask),pkR(veMask,2));
cVeE=getCorrCoef(estVe(veMask),pkE(veMask,2));
cVeCE=getCorrCoef(estVe(veMask),pkCE(veMask,2));

disp('Correlation with ve [RRM ERRM CERRM]')
disp([cVeR cVeE cVeCE])

% vp correlation coefficients
cVpE=getCorrCoef(pkT(corrMask,3),pkE(corrMask,4));
cVpCE=getCorrCoef(pkT(corrMask,3),pkCE(corrMask,4));
disp('Correlation with vp [RRM ERRM CERRM]')
disp([NaN cVpE cVpCE])

%% Make the maps
% To display the maps, scroll down to find code for producing maps
% Alternatively, use `showMe()` in the terminal, e.g:
% > showMe(kepC) 
% will show the kep maps produced by the CERRM
ktT = zeros(size(mask));
ktE = ktT;
ktC = ktT;
ktR = ktT;

kepT = ktT;
kepE = ktT;
kepC = ktT;
kepR = ktT;

vpT = ktT;
vpE = ktT;
vpC = ktT;

veE = ktT;
veT = ktT;
veC = ktT;
veR = ktT;

% ETM
ktT(mask(:))=pkT(:,1);
kepT(mask(:))=pkT(:,2);
vpT(mask(:))=pkT(:,3);
veT(mask(:))=pkT(:,1)./pkT(:,2);

% ERRM
ktE(mask(:))=pkE(:,1);
kepE(mask(:))=pkE(:,3);
vpE(mask(:))=pkE(:,4);
veE(mask(:))=pkE(:,2);

% CERRM
ktC(mask(:))=pkCE(:,1);
veC(mask(:))=pkCE(:,2);
kepC(mask(:))=pkCE(:,3);
vpC(mask(:))=pkCE(:,4);

% RRM
ktR(mask(:))=pkR(:,1);
veR(mask(:))=pkR(:,2);
kepR(mask(:))=pkR(:,3);

% Crop the maps
kepR = AutoCrop(kepR);
kepE = AutoCrop(kepE);
kepC = AutoCrop(kepC);
kepT = AutoCrop(kepT);

ktR = AutoCrop(ktR);
ktE = AutoCrop(ktE);
ktC = AutoCrop(ktC);
ktT = AutoCrop(ktT);

veR = AutoCrop(veR);
veE = AutoCrop(veE);
veC = AutoCrop(veC);
veT = AutoCrop(veT);

vpE = AutoCrop(vpE);
vpC = AutoCrop(vpC);
vpT = AutoCrop(vpT);
vpR = zeros(size(vpE)); % RRM doesn't have vp map

%% POTENTIAL END OF FILE
% Continue and show figures if that's what user wants
if ~openFloodgates
    return
end

%% Plot Cp and Crr

figure
subplot(2,1,1)
plot(t,Cp,'r','LineWidth',2)
xlabel('Time [min]')
ylabel('Concentration [mM]')
title('Blood plasma')

subplot(2,1,2)
plot(t,Crr,'g','LineWidth',2)
xlabel('Time [min]')
ylabel('Concentration [mM]')
title('Reference Region')

%% Produce scatter plots of estimates

% Depending on user's screen, these figures might appear squished.
% Resize manually to unsquish them.

% Size of the scatter plots
sz = 5;

% Colours
listColours{1} = [255,44,127]./255; % RRM - pink/red
listColours{2} = [255,196,68]./255; % ERRM - orange/gold
listColours{3} = [126,47,142]./255; % CERRM - purple

% %%% KTrans
% Define x- and y-axis range
mX = 0.5;
mY = 20;
x = 0:mX:mX;

% Find points that are within the x- and y-axis range
curM1 = pkT(:,1)<mX & pkT(:,1) > 0;
curM2a = pkR(:,1)<mY & pkR(:,1) > 0;
curM2b = pkE(:,1)<mY & pkE(:,1) > 0;
curM2c = pkCE(:,1)<mY & pkCE(:,1) > 0;

% Calculate slope using linear least squares
m1=pkT(curM1&curM2a,1)\pkR(curM1&curM2a,1);
m2=pkT(curM1&curM2b,1)\pkE(curM1&curM2b,1);
m3=pkT(curM1&curM2c,1)\pkCE(curM1&curM2c,1);

% Plot figure
figure
subplot(1,3,1)
scatter(pkT(curM1&curM2a,1),pkR(curM1&curM2a,1),sz,listColours{1})
hold on; plot(x,x*m1);
ylim([0 mY])
xlim([0 mX])
title('Scatter - Ktrans - RRM vs ETM')
xlabel('Ktrans - ETM')
ylabel('Ktrans - RRM')

subplot(1,3,2)
scatter(pkT(curM1&curM2b,1),pkE(curM1&curM2b,1),sz,listColours{2})
hold on; plot(x,x*m2);
ylim([0 mY])
xlim([0 mX])
title('Scatter - Ktrans - ERRM vs ETM')
xlabel('Ktrans - ETM')
ylabel('Ktrans - ERRM')

subplot(1,3,3)
scatter(pkT(curM1&curM2c,1),pkCE(curM1&curM2c,1),sz,listColours{3})
hold on; plot(x,x*m3);
ylim([0 mY])
xlim([0 mX])
title('Scatter - Ktrans - CERRM vs ETM')
xlabel('Ktrans - ETM')
ylabel('Ktrans - CERRM')

% %%% kep
% Same steps as for Ktrans
mX = 2;
mY = 2;
x = 0:mX:mX;

curM1 = pkT(:,2)<mX & pkT(:,2) > 0;
curM2a = pkR(:,3)<mY & pkR(:,3) > 0;
curM2b = pkE(:,3)<mY & pkE(:,3) > 0;
curM2c = pkCE(:,3)<mY & pkCE(:,3) > 0;

m1=pkT(curM1&curM2a,2)\pkR(curM1&curM2a,3);
m2=pkT(curM1&curM2b,2)\pkE(curM1&curM2b,3);
m3=pkT(curM1&curM2c,2)\pkCE(curM1&curM2c,3);

figure
subplot(1,3,1)
scatter(pkT(curM1&curM2a,2),pkR(curM1&curM2a,3),sz,listColours{1})
hold on; plot(x,x*m1);
ylim([0 mY])
xlim([0 mX])
title('Scatter - kep - RRM vs ETM')
xlabel('kep - ETM')
ylabel('kep - RRM')

subplot(1,3,2)
scatter(pkT(curM1&curM2b,2),pkE(curM1&curM2b,3),sz,listColours{2})
hold on; plot(x,x*m2);
ylim([0 mY])
xlim([0 mX])
title('Scatter - kep - ERRM vs ETM')
xlabel('kep - ETM')
ylabel('kep - ERRM')

subplot(1,3,3)
scatter(pkT(curM1&curM2c,2),pkCE(curM1&curM2c,3),sz,listColours{3})
hold on; plot(x,x*m3);
ylim([0 mY])
xlim([0 mX])
title('Scatter - kep - CERRM vs ETM')
xlabel('kep - ETM')
ylabel('kep - CERRM')

% %%% ve
mX = 0.5;
mY = 10;
x = 0:mX:mX;

% Calculate ve from Extended Tofts model (Ktrans/kep)
estVeT = pkT(:,1)./pkT(:,2);

curM1 = estVeT<mX & estVeT > 0;
curM2a = pkR(:,2)<mY & pkR(:,2) > 0;
curM2b = pkE(:,2)<mY & pkE(:,2) > 0;
curM2c = pkCE(:,2)<mY & pkCE(:,2) > 0;

m1=estVeT(curM1&curM2a)\pkR(curM1&curM2a,2);
m2=estVeT(curM1&curM2b)\pkE(curM1&curM2b,2);
m3=estVeT(curM1&curM2c)\pkCE(curM1&curM2c,2);

figure
subplot(1,3,1)
scatter(estVeT(curM1&curM2a),pkR(curM1&curM2a,2),sz,listColours{1})
hold on; plot(x,x*m1);
ylim([0 mY])
xlim([0 mX])
title('Scatter - ve - RRM vs ETM')
xlabel('ve - ETM')
ylabel('ve - RRM')

subplot(1,3,2)
scatter(estVeT(curM1&curM2b),pkE(curM1&curM2b,2),sz,listColours{2})
hold on; plot(x,x*m2);
ylim([0 mY])
xlim([0 mX])
title('Scatter - ve - ERRM vs ETM')
xlabel('ve - ETM')
ylabel('ve - ERRM')

subplot(1,3,3)
scatter(estVeT(curM1&curM2c),pkCE(curM1&curM2c,2),sz,listColours{3})
hold on; plot(x,x*m3);
ylim([0 mY])
xlim([0 mX])
title('Scatter - ve - CERRM vs ETM')
xlabel('ve - ETM')
ylabel('ve - CERRM')

% %%% vp
mX = 0.1;
mY = 2;
x = 0:mX:mX;

curM1 = pkT(:,3)<mX & pkT(:,3) > 0;
curM2b = pkE(:,4)<mY & pkE(:,4) > 0;
curM2c = pkCE(:,4)<mY & pkCE(:,4) > 0;

m2=pkT(curM1&curM2b,3)\pkE(curM1&curM2b,4);
m3=pkT(curM1&curM2c,3)\pkCE(curM1&curM2c,4);

figure
subplot(1,2,1)
scatter(pkT(curM1&curM2b,3),pkE(curM1&curM2b,4),sz,listColours{2})
hold on; plot(x,x*m2);
ylim([0 2])
xlim([0 0.1])
title('Scatter - vp - ERRM vs ETM')
xlabel('vp - ETM')
ylabel('vp - ERRM')

subplot(1,2,2)
scatter(pkT(curM1&curM2c,3),pkCE(curM1&curM2c,4),sz,listColours{3})
hold on; plot(x,x*m3);
ylim([0 2])
xlim([0 0.1])
title('Scatter - vp - CERRM vs ETM')
xlabel('vp - ETM')
ylabel('vp - CERRM')

%% Produce in-vivo maps of kep

% These settings are optimized for the first patient
% They'll likely have to be changed for the other two

% Pick which slices to show
mySlices = 7:9; % Manuscript shows slices 7 to 9 of patient 1 (or A)
% Set the range of the colorbar - Manuscript uses range of 0 to 1 for kep
climR = [0 1]; % % For RRM, ERRM and CERRM
climT = [0 1]; % For ExtToftsModel
% Pick which colormap to use - manuscript uses 'jet '
cMap = 'jet'; % 

% This function should take care of the rest
makeGridMaps({kepR,kepE,kepC,kepT},mySlices,{'RRM','ERRM','CERRM','ETM'},...
    {climR,climR,climR,climT},cMap)

% Plots title - might not work for some users if they're missing toolboxes
suptitle('k_{ep}')

%% Produce in-vivo maps of Ktrans
% Refer to previous cell (for kep) for comments

mySlices = 7:9;
climR = [0 3];
climT = [0 0.1];
cMap = 'jet'; 

makeGridMaps({ktR,ktE,ktC,ktT},mySlices,{'RRM','ERRM','CERRM','ETM'},...
    {climR,climR,climR,climT},cMap)
suptitle('K^{trans}')

%% Produce in-vivo maps of ve
% Refer to previous cell (for kep) for comments

mySlices = 7:9;
climR = [0 5];
climT = [0 0.35];

makeGridMaps({veR,veE,veC,veT},mySlices,{'RRM','ERRM','CERRM','ETM'},...
    {climR,climR,climR,climT},cMap)
suptitle('v_e')

%% Produce in-vivo maps of vp
% Refer to previous cell (for kep) for comments

mySlices = 7:9;
climR = [0 1];
climT = [0 0.05];

makeGridMaps({vpR,vpE,vpC,vpT},mySlices,{'RRM','ERRM','CERRM','ETM'},...
    {climR,climR,climR,climT},cMap)
suptitle('v_p')

%% Histogram of kepRR estimates
% This wasn't show in the manuscript

% Define the max/min of the x-axis for the histogram
binLim = 3;

% The rest should work on its own
% Histogram should show negative estimates as a different colour.
x=pkE(corrMask,5);

figure
histogram(x,70,'BinLimits',[-binLim binLim])
hold on
histogram(x,35,'BinLimits',[-binLim 0])
title('Distribution of kepRR estimates from ERRM')
xlabel('kepRR estimate from ERRM [min^{-1}]')
ylabel('Counts')