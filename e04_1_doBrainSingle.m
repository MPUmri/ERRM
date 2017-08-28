% Process the GBM data with RRM, ERRM, CERRM, and ETM

% Estimated runtime: ~5 seconds

%% Pre-setup
clearvars
addpath('./mfiles')

%5 Load data
load('./data/gbm/gbmData.mat');
% The loaded data contains:
% t   [70x1] - time (in minutes) at each dce frame
% Cp  [70x1] - concentration in blood plasma
% Crr [70x1] - concentration in reference region
% Ct  [70x10813] - concentration in tumour voxels (10813 voxels)
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
[pkCE, fitE, estKepRR_CE, pkE, fitE]=CERRM(Ct,Crr,t);

% % ERRM
% tic
% pkE=ERRM(Ct,Crr,t);
% toc
% % CERRM
% tic
% [pkCE r estKepRR_CE]=CERRM(Ct,Crr,t);
% toc


%% Fit with non-extended RRM
tic
pkR=LRRM(Ct,Crr,t);
toc

%% Fit with ETM
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
vpRR = pkCrr(3)

%% Produce the maps
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

%% Correlations

% Calculate ve from Extended Tofts model (Ktrans/kep)
estVeT = pkT(:,1)./pkT(:,2);

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

% Calculate slow using linear least squares
m1=pkT(curM1&curM2a,1)\pkR(curM1&curM2a,1);
m2=pkT(curM1&curM2b,1)\pkE(curM1&curM2b,1);
m3=pkT(curM1&curM2c,1)\pkCE(curM1&curM2c,1);

% Calculate Pearson correlation coefficient
[c1,p1] = getCorrCoef(pkT(curM1&curM2a,1),pkR(curM1&curM2a,1));
[c2,p2] = getCorrCoef(pkT(curM1&curM2b,1),pkE(curM1&curM2b,1));
[c3,p3] = getCorrCoef(pkT(curM1&curM2c,1),pkCE(curM1&curM2c,1));

% Plot figure
figure
scatter(pkT(curM1&curM2a,1),pkR(curM1&curM2a,1),sz,listColours{1})
hold on; plot(x,x*m1);
ylim([0 mY])
xlim([0 mX])
title('Scatter - Ktrans - RRM vs ETM')
xlabel('Ktrans - ETM')
ylabel('Ktrans - RRM')

figure
scatter(pkT(curM1&curM2b,1),pkE(curM1&curM2b,1),sz,listColours{2})
hold on; plot(x,x*m2);
ylim([0 mY])
xlim([0 mX])
title('Scatter - Ktrans - ERRM vs ETM')
xlabel('Ktrans - ETM')
ylabel('Ktrans - ERRM')

figure
scatter(pkT(curM1&curM2c,1),pkCE(curM1&curM2c,1),sz,listColours{3})
hold on; plot(x,x*m3);
ylim([0 mY])
xlim([0 mX])
title('Scatter - Ktrans - CERRM vs ETM')
xlabel('Ktrans - ETM')
ylabel('Ktrans - CERRM')

slopesKt = [m1 m2 m3]
corrKt = [c1 c2 c3]
pOfCorrKt = [p1 p2 p3]

% [polyfit(pkT(curM1&curM2a,1),pkR(curM1&curM2a,1),1)...
%     polyfit(pkT(curM1&curM2a,1),pkE(curM1&curM2a,1),1)...
%     polyfit(pkT(curM1&curM2a,1),pkCE(curM1&curM2a,1),1)]

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

[c1,p1] = getCorrCoef(pkT(curM1&curM2a,2),pkR(curM1&curM2a,3));
[c2,p2] = getCorrCoef(pkT(curM1&curM2b,2),pkE(curM1&curM2b,3));
[c3,p3] = getCorrCoef(pkT(curM1&curM2c,2),pkCE(curM1&curM2c,3));

figure
scatter(pkT(curM1&curM2a,2),pkR(curM1&curM2a,3),sz,listColours{1})
hold on; plot(x,x*m1);
ylim([0 mY])
xlim([0 mX])
title('Scatter - kep - RRM vs ETM')
xlabel('kep - ETM')
ylabel('kep - RRM')

figure
scatter(pkT(curM1&curM2b,2),pkE(curM1&curM2b,3),sz,listColours{2})
hold on; plot(x,x*m2);
ylim([0 mY])
xlim([0 mX])
title('Scatter - kep - ERRM vs ETM')
xlabel('kep - ETM')
ylabel('kep - ERRM')

figure
scatter(pkT(curM1&curM2c,2),pkCE(curM1&curM2c,3),sz,listColours{3})
hold on; plot(x,x*m3);
ylim([0 mY])
xlim([0 mX])
title('Scatter - kep - CERRM vs ETM')
xlabel('kep - ETM')
ylabel('kep - CERRM')

slopesKep = [m1 m2 m3]
corrKep = [c1 c2 c3]
pOfCorrKep = [p1 p2 p3]

% [polyfit(pkT(curM1&curM2a,2),pkR(curM1&curM2a,3),1)...
%     polyfit(pkT(curM1&curM2a,2),pkE(curM1&curM2a,3),1)...
%     polyfit(pkT(curM1&curM2a,2),pkCE(curM1&curM2a,3),1)]

% %%% ve
mX = 0.5;
mY = 10;
x = 0:mX:mX;

curM1 = estVeT<mX & estVeT > 0;
curM2a = pkR(:,2)<mY & pkR(:,2) > 0;
curM2b = pkE(:,2)<mY & pkE(:,2) > 0;
curM2c = pkCE(:,2)<mY & pkCE(:,2) > 0;

m1=estVeT(curM1&curM2a)\pkR(curM1&curM2a,2);
m2=estVeT(curM1&curM2b)\pkE(curM1&curM2b,2);
m3=estVeT(curM1&curM2c)\pkCE(curM1&curM2c,2);

[c1,p1] = getCorrCoef(estVeT(curM1&curM2a),pkR(curM1&curM2a,2));
[c2,p2] = getCorrCoef(estVeT(curM1&curM2b),pkE(curM1&curM2b,2));
[c3,p3] = getCorrCoef(estVeT(curM1&curM2c),pkCE(curM1&curM2c,2));

figure
scatter(estVeT(curM1&curM2a),pkR(curM1&curM2a,2),sz,listColours{1})
hold on; plot(x,x*m1);
ylim([0 mY])
xlim([0 mX])
title('Scatter - ve - RRM vs ETM')
xlabel('ve - ETM')
ylabel('ve - RRM')

figure
scatter(estVeT(curM1&curM2b),pkE(curM1&curM2b,2),sz,listColours{2})
hold on; plot(x,x*m2);
ylim([0 mY])
xlim([0 mX])
title('Scatter - ve - ERRM vs ETM')
xlabel('ve - ETM')
ylabel('ve - ERRM')

figure
scatter(estVeT(curM1&curM2c),pkCE(curM1&curM2c,2),sz,listColours{3})
hold on; plot(x,x*m3);
ylim([0 mY])
xlim([0 mX])
title('Scatter - ve - CERRM vs ETM')
xlabel('ve - ETM')
ylabel('ve - CERRM')

slopesVe = [m1 m2 m3]
corrVe = [c1 c2 c3]
pOfCorrVe = [p1 p2 p3]

% [polyfit(estVeT(curM1&curM2a),pkR(curM1&curM2a,2),1)...
%     polyfit(estVeT(curM1&curM2a),pkE(curM1&curM2a,2),1)...
%     polyfit(estVeT(curM1&curM2a),pkCE(curM1&curM2a,2),1)]

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
scatter(pkT(curM1&curM2b,3),pkE(curM1&curM2b,4),sz,listColours{2})
hold on; plot(x,x*m2);
ylim([0 2])
xlim([0 0.1])
title('Scatter - vp - ERRM vs ETM')
xlabel('vp - ETM')
ylabel('vp - ERRM')

figure
scatter(pkT(curM1&curM2c,3),pkCE(curM1&curM2c,4),sz,listColours{3})
hold on; plot(x,x*m3);
ylim([0 2])
xlim([0 0.1])
title('Scatter - vp - CERRM vs ETM')
xlabel('vp - ETM')
ylabel('vp - CERRM')

%%
mySlices = 7:9;

%% kep
kepR = AutoCrop(kepR);
kepE = AutoCrop(kepE);
kepC = AutoCrop(kepC);
kepT = AutoCrop(kepT);

climR = [0 1];
climT = [0 1];

figure
subplot(4,3,1)
imagesc(kepR(:,:,mySlices(1)),climR);
axis square
ylabel('RRM')
subplot(4,3,2)
imagesc(kepR(:,:,mySlices(2)),climR);
axis square
subplot(4,3,3)
imagesc(kepR(:,:,mySlices(3)),climR);
axis square

subplot(4,3,4)
imagesc(kepE(:,:,mySlices(1)),climR);
axis square
ylabel('ERRM')
subplot(4,3,5)
imagesc(kepE(:,:,mySlices(2)),climR);
axis square
subplot(4,3,6)
imagesc(kepE(:,:,mySlices(3)),climR);
axis square

subplot(4,3,7)
imagesc(kepC(:,:,mySlices(1)),climR);
axis square
ylabel('CERRM')
subplot(4,3,8)
imagesc(kepC(:,:,mySlices(2)),climR);
axis square
subplot(4,3,9)
imagesc(kepC(:,:,mySlices(3)),climR);
axis square

subplot(4,3,10)
imagesc(kepT(:,:,mySlices(1)),climT);
axis square
ylabel('ETM')
subplot(4,3,11)
imagesc(kepT(:,:,mySlices(2)),climT);
axis square
subplot(4,3,12)
imagesc(kepT(:,:,mySlices(3)),climT);
axis square

% Plots title - might not work for some users if they're missing toolboxes
suptitle('k_{ep}')
%% KTrans
ktR = AutoCrop(ktR);
ktE = AutoCrop(ktE);
ktC = AutoCrop(ktC);
ktT = AutoCrop(ktT);

climR = [0 3];
climT = [0 .1];

figure
subplot(4,3,1)
imagesc(ktR(:,:,mySlices(1)),climR);
axis square
ylabel('RRM')
subplot(4,3,2)
imagesc(ktR(:,:,mySlices(2)),climR);
axis square
subplot(4,3,3)
imagesc(ktR(:,:,mySlices(3)),climR);
axis square

subplot(4,3,4)
imagesc(ktE(:,:,mySlices(1)),climR);
axis square
ylabel('ERRM')
subplot(4,3,5)
imagesc(ktE(:,:,mySlices(2)),climR);
axis square
subplot(4,3,6)
imagesc(ktE(:,:,mySlices(3)),climR);
axis square

subplot(4,3,7)
imagesc(ktC(:,:,mySlices(1)),climR);
axis square
ylabel('CERRM')
subplot(4,3,8)
imagesc(ktC(:,:,mySlices(2)),climR);
axis square
subplot(4,3,9)
imagesc(ktC(:,:,mySlices(3)),climR);
axis square

subplot(4,3,10)
imagesc(ktT(:,:,mySlices(1)),climT);
axis square
ylabel('ETM')
subplot(4,3,11)
imagesc(ktT(:,:,mySlices(2)),climT);
axis square
subplot(4,3,12)
imagesc(ktT(:,:,mySlices(3)),climT);
axis square

suptitle('K^{trans}')
%% ve
veR = AutoCrop(veR);
veE = AutoCrop(veE);
veC = AutoCrop(veC);
veT = AutoCrop(veT);

climR = [0 5];
climT = [0 0.35];

figure
subplot(4,3,1)
imagesc(veR(:,:,mySlices(1)),climR);
axis square
ylabel('RRM')
subplot(4,3,2)
imagesc(veR(:,:,mySlices(2)),climR);
axis square
subplot(4,3,3)
imagesc(veR(:,:,mySlices(3)),climR);
axis square

subplot(4,3,4)
imagesc(veE(:,:,mySlices(1)),climR);
axis square
ylabel('ERRM')
subplot(4,3,5)
imagesc(veE(:,:,mySlices(2)),climR);
axis square
subplot(4,3,6)
imagesc(veE(:,:,mySlices(3)),climR);
axis square

subplot(4,3,7)
imagesc(veC(:,:,mySlices(1)),climR);
axis square
ylabel('CERRM')
subplot(4,3,8)
imagesc(veC(:,:,mySlices(2)),climR);
axis square
subplot(4,3,9)
imagesc(veC(:,:,mySlices(3)),climR);
axis square

subplot(4,3,10)
imagesc(veT(:,:,mySlices(1)),climT);
axis square
ylabel('ETM')
subplot(4,3,11)
imagesc(veT(:,:,mySlices(2)),climT);
axis square
subplot(4,3,12)
imagesc(veT(:,:,mySlices(3)),climT);
axis square

suptitle('v_e')
%% vp
vpE = AutoCrop(vpE);
vpC = AutoCrop(vpC);
vpT = AutoCrop(vpT);
vpR = zeros(size(vpE));

climR = [0 1];
climT = [0 0.05];

figure
subplot(4,3,1)
imagesc(vpR(:,:,mySlices(1)),climR);
axis square
ylabel('RRM')
subplot(4,3,2)
imagesc(vpR(:,:,mySlices(2)),climR);
axis square
subplot(4,3,3)
imagesc(vpR(:,:,mySlices(3)),climR);
axis square

subplot(4,3,4)
imagesc(vpE(:,:,mySlices(1)),climR);
axis square
ylabel('ERRM')
subplot(4,3,5)
imagesc(vpE(:,:,mySlices(2)),climR);
axis square
subplot(4,3,6)
imagesc(vpE(:,:,mySlices(3)),climR);
axis square

subplot(4,3,7)
imagesc(vpC(:,:,mySlices(1)),climR);
axis square
ylabel('CERRM')
subplot(4,3,8)
imagesc(vpC(:,:,mySlices(2)),climR);
axis square
subplot(4,3,9)
imagesc(vpC(:,:,mySlices(3)),climR);
axis square

subplot(4,3,10)
imagesc(vpT(:,:,mySlices(1)),climT);
axis square
ylabel('ETM')
subplot(4,3,11)
imagesc(vpT(:,:,mySlices(2)),climT);
axis square
subplot(4,3,12)
imagesc(vpT(:,:,mySlices(3)),climT);
axis square

suptitle('v_p')
%% Plot Cp and Crr

figure('Position',[300 300 500 600])
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
