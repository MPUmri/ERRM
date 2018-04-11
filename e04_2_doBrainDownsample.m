% Downsample GBM data to compare stability of CERRM and ETM

% Estimated runtime: 65 seconds

%% Pre-setup

clearvars
addpath('./mfiles')
rng(12345)

%% Load data
load('./data/gbm/TCGA-06-0185-1.mat');
t1=t;
Cp1=Cp;
Crr1=Crr;
Ct1=Ct(:,max(Ct)'>0.05 & max(Ct)'<1);
load('./data/gbm/TCGA-06-0185-2.mat');
t2=t;
Cp2=Cp;
Crr2=Crr;
Ct2=Ct(:,max(Ct)'>0.05 & max(Ct)'<1);
load('./data/gbm/TCGA-06-0881-1.mat');
t3=t;
Cp3=Cp;
Crr3=Crr;
Ct3=Ct(:,max(Ct)'>0.05 & max(Ct)'<1);
% The loaded data contains:
% t   [70x1] - time (in minutes) at each dce frame
% Cp  [70x1] - concentration in blood plasma
% Crr [70x1] - concentration in reference region
% Ct  [70x10813] - concentration in tumour voxels (10813 voxels)

tic

%% Fit to high temporal resolution (5.4 s - original) data

% CERRM
[tempPkCE1, ~, estKepRR1(1)]=CERRM(Ct1,Crr1,t1);
[tempPkCE2, ~, estKepRR2(1)]=CERRM(Ct2,Crr2,t2);
[tempPkCE3, ~, estKepRR3(1)]=CERRM(Ct3,Crr3,t3);
% ETM
tempPkET1=LLSQ(Ct1,Cp1,t1,1);
tempPkET2=LLSQ(Ct2,Cp2,t2,1);
tempPkET3=LLSQ(Ct3,Cp3,t3,1);

% 'true' parameters are fits at original temporal resolution
trueKtCE = [tempPkCE1(:,1); tempPkCE2(:,1); tempPkCE3(:,1)];
trueVeCE = [tempPkCE1(:,2); tempPkCE2(:,2); tempPkCE3(:,2)];
trueKepCE = [tempPkCE1(:,3); tempPkCE2(:,3); tempPkCE3(:,3)];
trueVpCE = [tempPkCE1(:,4); tempPkCE2(:,4); tempPkCE3(:,4)];

trueKtET = [tempPkET1(:,1); tempPkET2(:,1); tempPkET3(:,1)];
trueKepET = [tempPkET1(:,2); tempPkET2(:,2); tempPkET3(:,2)];
trueVeET = trueKtET ./ trueKepET;
trueVpET = [tempPkET1(:,3); tempPkET2(:,3); tempPkET3(:,3)];

TRes(1) = 60*(t(2)-t(1));

% Save the fitted parameters (for first patient - Figure 6)
rawFitsET = zeros([size(tempPkET1) 12]);
rawFitsCE = zeros([size(tempPkCE1) 12]);

rawFitsET(:,:,1) = tempPkET1;
rawFitsCE(:,:,1) = tempPkCE1;

%% Downsample the data and refit
nVox1 = size(Ct1,2);
nVox2 = size(Ct2,2);
nVox3 = size(Ct3,2);
for i=2:12  
    % Randomize the first frame when downsampling
    phaseValues1 = randi([0 i-1],nVox1,1);
    phaseValues2 = randi([0 i-1],nVox2,1);
    phaseValues3 = randi([0 i-1],nVox3,1);
    clearvars tempPkCE tempPkET pkERRM
    for j=1:nVox1
        curCt = downsample(Ct1(:,j),i,phaseValues1(j));
        curCrr = downsample(Crr1,i,phaseValues1(j));
        curCp = downsample(Cp1,i,phaseValues1(j));
        curT = downsample(t1,i,phaseValues1(j));
        
        pkERRM1(j,:)=ERRM(curCt,curCrr,curT,0,0,0); % ERRM
        tempPkET1(j,:)=LLSQ(curCt,curCp,curT,1); % ETM
    end
    for j=1:nVox2
        curCt = downsample(Ct2(:,j),i,phaseValues2(j));
        curCrr = downsample(Crr2,i,phaseValues2(j));
        curCp = downsample(Cp2,i,phaseValues2(j));
        curT = downsample(t2,i,phaseValues2(j));
        
        pkERRM2(j,:)=ERRM(curCt,curCrr,curT,0,0,0); % ERRM
        tempPkET2(j,:)=LLSQ(curCt,curCp,curT,1); % ETM
    end
    for j=1:nVox3
        curCt = downsample(Ct3(:,j),i,phaseValues3(j));
        curCrr = downsample(Crr3,i,phaseValues3(j));
        curCp = downsample(Cp3,i,phaseValues3(j));
        curT = downsample(t3,i,phaseValues3(j));
        
        pkERRM3(j,:)=ERRM(curCt,curCrr,curT,0,0,0); % ERRM
        tempPkET3(j,:)=LLSQ(curCt,curCp,curT,1); % ETM
    end
    
    % Get estimate for kepRR - from ERRM
    pkERRM = pkERRM1;
    rawKepRR = pkERRM(:,5);
    % Find voxels where all estimates are real and positive
    goodVals = pkERRM(:,1)>0 & pkERRM(:,2)>0 & pkERRM(:,3)>0 & pkERRM(:,4)>0 & pkERRM(:,5)>0 & imag(pkERRM(:,5))==0;
    estKepRR1(i) = iqrMean(rawKepRR(goodVals));
    
    pkERRM = pkERRM2;
    rawKepRR = pkERRM(:,5);
    % Find voxels where all estimates are real and positive
    goodVals = pkERRM(:,1)>0 & pkERRM(:,2)>0 & pkERRM(:,3)>0 & pkERRM(:,4)>0 & pkERRM(:,5)>0 & imag(pkERRM(:,5))==0;
    estKepRR2(i) = iqrMean(rawKepRR(goodVals));
    
    pkERRM = pkERRM3;
    rawKepRR = pkERRM(:,5);
    % Find voxels where all estimates are real and positive
    goodVals = pkERRM(:,1)>0 & pkERRM(:,2)>0 & pkERRM(:,3)>0 & pkERRM(:,4)>0 & pkERRM(:,5)>0 & imag(pkERRM(:,5))==0;
    estKepRR3(i) = iqrMean(rawKepRR(goodVals));
    
    % Now fit CERRM
    for j=1:nVox1
        curCt = downsample(Ct1(:,j),i,phaseValues1(j));
        curCrr = downsample(Crr1,i,phaseValues1(j));
        curCp = downsample(Cp1,i,phaseValues1(j));
        curT = downsample(t1,i,phaseValues1(j));
        
        tempPkCE1(j,:)=CERRM(curCt,curCrr,curT,estKepRR1(i));
    end
    for j=1:nVox2
        curCt = downsample(Ct2(:,j),i,phaseValues2(j));
        curCrr = downsample(Crr2,i,phaseValues2(j));
        curCp = downsample(Cp2,i,phaseValues2(j));
        curT = downsample(t2,i,phaseValues2(j));
        
        tempPkCE2(j,:)=CERRM(curCt,curCrr,curT,estKepRR2(i));
    end
    for j=1:nVox3
        curCt = downsample(Ct3(:,j),i,phaseValues3(j));
        curCrr = downsample(Crr3,i,phaseValues3(j));
        curCp = downsample(Cp3,i,phaseValues3(j));
        curT = downsample(t3,i,phaseValues3(j));
        
        tempPkCE3(j,:)=CERRM(curCt,curCrr,curT,estKepRR3(i));
    end
    
    rawFitsET(:,:,i) = tempPkET1;
    rawFitsCE(:,:,i) = tempPkCE1;
    
    % Calculate percent change
    tempErrKtCE=PercentError([tempPkCE1(:,1); tempPkCE2(:,1); tempPkCE3(:,1)],trueKtCE);
    tempErrVeCE=PercentError([tempPkCE1(:,2); tempPkCE2(:,2); tempPkCE3(:,2)],trueVeCE);
    tempErrKepCE=PercentError([tempPkCE1(:,3); tempPkCE2(:,3); tempPkCE3(:,3)],trueKepCE);
    tempErrVpCE=PercentError([tempPkCE1(:,4); tempPkCE2(:,4); tempPkCE3(:,4)],trueVpCE);
        
    tempErrKtET=PercentError([tempPkET1(:,1); tempPkET2(:,1); tempPkET3(:,1)],trueKtET);
    tempErrVeET=PercentError([tempPkET1(:,1); tempPkET2(:,1); tempPkET3(:,1)]./[tempPkET1(:,2); tempPkET2(:,2); tempPkET3(:,2)],trueVeET);
    tempErrKepET=PercentError([tempPkET1(:,2); tempPkET2(:,2); tempPkET3(:,2)],trueKepET);
    tempErrVpET=PercentError([tempPkET1(:,3); tempPkET2(:,3); tempPkET3(:,3)],trueVpET);

    % Get quartiles
    % (Order is not efficient - could preallocate to improve speed)
    qKtCE(i,:)=quantile(tempErrKtCE(:),[.25 .5 .75]);
    qVeCE(i,:)=quantile(tempErrVeCE(:),[.25 .5 .75]);
    qKepCE(i,:)=quantile(tempErrKepCE(:),[.25 .5 .75]);
    qVpCE(i,:)=quantile(tempErrVpCE(:),[.25 .5 .75]);
  
    qKtET(i,:)=quantile(tempErrKtET(:),[.25 .5 .75]);
    qVeET(i,:)=quantile(tempErrVeET(:),[.25 .5 .75]);
    qKepET(i,:)=quantile(tempErrKepET(:),[.25 .5 .75]);
    qVpET(i,:)=quantile(tempErrVpET(:),[.25 .5 .75]);
    
    TRes(i) = 60*(curT(2)-curT(1));
end
toc

%% Plot: Median +/- interquartile range vs Temporal Resolution
lErr = 1;
hErr = 3;

cCE = [.5 .2 .55];
cET = [.3 .7 .6];

figure('Position',[300 300 1000 300])
subplot(1,4,1)
hold on
shadedErrorPlot(TRes',qKtET(:,2),qKtET(:,lErr),qKtET(:,hErr),cET,1);
shadedErrorPlot(TRes',qKtCE(:,2),qKtCE(:,lErr),qKtCE(:,hErr),cCE,1);
hold off
ylim([-100 100])
xlabel('Temporal Resolution [s]')
ylabel('Percent Error')
title('Ktrans')

subplot(1,4,2)
hold on
shadedErrorPlot(TRes',qKepET(:,2),qKepET(:,lErr),qKepET(:,hErr),cET,1);
shadedErrorPlot(TRes',qKepCE(:,2),qKepCE(:,lErr),qKepCE(:,hErr),cCE,1);
hold off
ylim([-100 100])
xlabel('Temporal Resolution [s]')
ylabel('Percent Error')
title('kep')

subplot(1,4,3)
hold on
shadedErrorPlot(TRes',qVeET(:,2),qVeET(:,lErr),qVeET(:,hErr),cET,1);
shadedErrorPlot(TRes',qVeCE(:,2),qVeCE(:,lErr),qVeCE(:,hErr),cCE,1);
hold off
ylim([-100 100])
xlabel('Temporal Resolution [s]')
ylabel('Percent Error')
title('ve')

subplot(1,4,4)
hold on
shadedErrorPlot(TRes',qVpET(:,2),qVpET(:,lErr),qVpET(:,hErr),cET,1);
shadedErrorPlot(TRes',qVpCE(:,2),qVpCE(:,lErr),qVpCE(:,hErr),cCE,1);
hold off
ylim([-100 300])
xlabel('Temporal Resolution [s]')
ylabel('Percent Error')
title('vp')
%% Plot downsampled maps (Fig 6) - Preparation

load('./data/gbm/TCGA-06-0185-1.mat');
mask = mask>0;
mask(mask) = max(Ct)'>0.05 & max(Ct)'<1;

[mapKtET, mapKepET, mapVeET, mapVpET, ...
    mapKtCE, mapKepCE, mapVeCE, mapVpCE] = deal(zeros([length(mask(:)) 12]));

for i=1:12
   mapKtET(mask(:),i) = rawFitsET(:,1,i); 
   mapVeET(mask(:),i) = rawFitsET(:,1,i)./rawFitsET(:,2,i);
   mapKepET(mask(:),i) = rawFitsET(:,2,i); 
   mapVpET(mask(:),i) = rawFitsET(:,3,i); 
   
   mapKtCE(mask(:),i) = rawFitsCE(:,1,i); 
   mapVeCE(mask(:),i) = rawFitsCE(:,2,i); 
   mapKepCE(mask(:),i) = rawFitsCE(:,3,i); 
   mapVpCE(mask(:),i) = rawFitsCE(:,4,i); 
end

mapKtET = reshape(mapKtET,[size(mask) 12]);
mapVeET = reshape(mapVeET,[size(mask) 12]);
mapKepET = reshape(mapKepET,[size(mask) 12]);
mapVpET = reshape(mapVpET,[size(mask) 12]);

mapKtCE = reshape(mapKtCE,[size(mask) 12]);
mapVeCE = reshape(mapVeCE,[size(mask) 12]);
mapKepCE = reshape(mapKepCE,[size(mask) 12]);
mapVpCE = reshape(mapVpCE,[size(mask) 12]);
%% Plot downsampled maps

% Pick which slice to plot
chosenSlice = 7; % Manuscript uses slice 7

% The temporal resolutions will be TRes(tresList)
tresList = [1,6,12];
% Code is designed to work for [1,6,12] - changing it will cause issues
% Specifically, the subplots will not show correctly.

for i=tresList
    curMapKtET = AutoCrop(mapKtET(:,:,:,i));
    curMapVeET = AutoCrop(mapVeET(:,:,:,i));
    curMapKepET = AutoCrop(mapKepET(:,:,:,i));
    curMapVpET = AutoCrop(mapVpET(:,:,:,i));

    curMapKtCE = AutoCrop(mapKtCE(:,:,:,i));
    curMapVeCE = AutoCrop(mapVeCE(:,:,:,i));
    curMapKepCE = AutoCrop(mapKepCE(:,:,:,i));
    curMapVpCE = AutoCrop(mapVpCE(:,:,:,i));
    
    figure(1)
    subplot(2,3,round(i/6)+1)
    imagesc(curMapKtET(:,:,chosenSlice))
    caxis([0 0.1]); colormap('jet')
    axis image; set(gca,'xtick',[]); set(gca,'ytick',[])
    if i==1
        ylabel('ETM-Ktrans')
    end
    xlabel([num2str(TRes(i)) ' s'])
    subplot(2,3,round(i/6)+4)
    imagesc(curMapKtCE(:,:,chosenSlice))
    caxis([0 3]); colormap('jet')
    axis image; set(gca,'xtick',[]); set(gca,'ytick',[])
    if i==1
        ylabel('CERRM-Ktrans')
    end
    xlabel([num2str(TRes(i)) ' s'])
    
    figure(2)
    subplot(2,3,round(i/6)+1)
    imagesc(curMapKepET(:,:,chosenSlice))
    caxis([0 1]); colormap('jet')
    axis image; set(gca,'xtick',[]); set(gca,'ytick',[])
    if i==1
        ylabel('ETM-kep')
    end
    xlabel([num2str(TRes(i)) ' s'])
    subplot(2,3,round(i/6)+4)
    imagesc(curMapKepCE(:,:,chosenSlice))
    caxis([0 1]); colormap('jet')
    axis image; set(gca,'xtick',[]); set(gca,'ytick',[])
    if i==1
        ylabel('CERRM-kep')
    end
    xlabel([num2str(TRes(i)) ' s'])
    
    figure(3)
    subplot(2,3,round(i/6)+1)
    imagesc(curMapVeET(:,:,chosenSlice))
    caxis([0 0.35]); colormap('jet')
    axis image; set(gca,'xtick',[]); set(gca,'ytick',[])
    if i==1
        ylabel('ETM-ve')
    end
    xlabel([num2str(TRes(i)) ' s'])
    subplot(2,3,round(i/6)+4)
    imagesc(curMapVeCE(:,:,chosenSlice))
    caxis([0 5]); colormap('jet')
    axis image; set(gca,'xtick',[]); set(gca,'ytick',[])
    if i==1
        ylabel('CERRM-ve')
    end
    xlabel([num2str(TRes(i)) ' s'])
    
    figure(4)
    subplot(2,3,round(i/6)+1)
    imagesc(curMapVpET(:,:,chosenSlice))
    caxis([0 0.05]); colormap('jet')
    axis image; set(gca,'xtick',[]); set(gca,'ytick',[])
    if i==1
        ylabel('ETM-vp')
    end
    xlabel([num2str(TRes(i)) ' s'])
    subplot(2,3,round(i/6)+4)
    imagesc(curMapVpCE(:,:,chosenSlice))
    caxis([0 1]); colormap('jet')
    axis image; set(gca,'xtick',[]); set(gca,'ytick',[])
    if i==1
        ylabel('CERRM-vp')
    end
    xlabel([num2str(TRes(i)) ' s'])
end

%% Plot all the input functions and reference region curves (SuppFigS4)

% Colours might be different between these and manuscript version

figure;

subplot(1,2,1)
plot(t,[Cp1 Cp2 Cp3], 'linewidth', 3)
ylim([0 16])
xlabel('Time [min]')
ylabel('Concentration [mM]')

subplot(1,2,2)
plot(t,[Crr1 Crr2 Crr3], 'linewidth', 3)
ylim([0 0.3])
xlabel('Time [min]')
ylabel('Concentration [mM]')

%% Export results to CSV
% This exports the median error (err), and 25th and 75th quantiles (qt25, qt75)
% to a CSV file.

% I don't think this ended up being used. But it's here incase someone
% finds it useful.

% Headers are:
% - errKt, errVe, etc - Median percent change
% - qt25Kt, qt25Ve, etc - 25th percentile (1st quartile) of percent change
% - qt75Kt, qt75Ve, etc - 75th percentile (3rd quartile) of percent change

outFile = fullfile('./dataResults',['e04-downsampleGbmResults.csv']);
hdr=['FitMethod,TemporalRes,errKt,errVe,errKep,errVp,qt25Kt,qt75Kt,qt25Ve,qt75Ve,qt25Kep,qt75Kep,qt25Vp,qt75Vp'];
outID = fopen(outFile, 'w+');
fprintf(outID, '%s\n', hdr); % Print header into csv file
for i=1:12
   outLine = {'ETM',TRes(i),qKtET(i,2),qVeET(i,2),qKepET(i,2),qVpET(i,2),...
       qKtET(i,1),qKtET(i,3),qVeET(i,1),qVeET(i,3),qKepET(i,1),qKepET(i,3),...
       qVpET(i,1),qVpET(i,3)};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'CERRM',TRes(i),qKtCE(i,2),qVeCE(i,2),qKepCE(i,2),qVpCE(i,2),...
       qKtCE(i,1),qKtCE(i,3),qVeCE(i,1),qVeCE(i,3),qKepCE(i,1),qKepCE(i,3),...
       qVpCE(i,1),qVpCE(i,3)};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
end
fclose(outID)