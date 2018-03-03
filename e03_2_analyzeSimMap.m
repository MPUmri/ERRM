% Adds noise and downsamples data of 120 perfusion parameter scombinations
% and then fits using reference region and tofts models.

% This is the main bulk of the simulation.
% It can take a couple of hours to run, and requires ~7 Gb of RAM.
% The RAM usage will probably decrease by replacing 'parfor' with 'for'

% The code is a bit wasteful because it does a couple of things which
% weren't actually included in the manuscript.
% For example, it fits the Reference Tissue Method, and the non-extended
% constrained reference region model. It also simulates temporal
% resolutions from 1-60s, but only 1-30s is shown in manuscript.

% Estimated runtime: ~ 1 to 3 hours, depending on settings

%% Pre-setup
clearvars
addpath('./mfiles')
rng(12345) % Arbitrary noise seed

%% Initial steps

% Try non-negative linear least fitting? (doubles the run time)
doNonNeg = false;

outDir = './data/mapDownsample'; % Output directory
outDirNN = './data/mapDownsample-NN'; % Output direction for non-negative fits

% Create output directories
if ~exist(outDir,'dir')
    mkdir(outDir);
    mkdir(fullfile(outDir,'CERRM'));
    mkdir(fullfile(outDir,'ETM'));
end
if doNonNeg && ~exist(outDirNN,'dir')
    mkdir(outDirNN);
    mkdir(fullfile(outDirNN,'CERRM'));
    mkdir(fullfile(outDirNN,'ETM'));
end

% Load noiseless simulated data
load('./data/simMap.mat');

% Get initial temporal resolution
sSize=initTRes;

% Range of noise to use (in units mM)
listSigmaC = 0:0.01:0.05;

% Number of replications for each noise level
repF = 100;

% Desired range of temporal resolutions (in seconds)
TRes = 1:60;

% Define and simulate the reference tissue
KtRR = 0.07; % units: 1/min
kepRR = 0.5; % units: 1/min
veRR = KtRR/kepRR;
CrrClean = ToftsKety(Cp,[KtRR,kepRR],t,0); % Noiseless reference region curve

% Replicate the noiseless data (don't add noise yet)
[sT, sX, sY] = size(simCt);
CtClean = reshape(simCt,[sT sX*sY]);
CtClean = repmat(CtClean,[1 repF]);
nVox=sX*sY*repF;

CpClean = repmat(Cp,[1 nVox]);
CrrClean = repmat(CrrClean,[1 nVox]);
% We now have 100 replications of noiseless data.
% Next, we'll add the required noise, apply downsampling, and performs fits

%% Main process
tic
% Cycle through all desired noise levels
for q=1:length(listSigmaC) 
    q % Let the console know something is being done
    sigmaC = listSigmaC(q);
    % Add noise to replicated noiseless data 
    Ct = CtClean + sigmaC*randn(size(CtClean));
    Crr = CrrClean + 0.1*sigmaC*randn(size(CrrClean));
    Cp = CpClean + sigmaC*randn(size(CpClean));
% Choosing not to indent the 2nd for-loop
for i=1:length(TRes)
    i % Let the console know where we're at
    dFactor=TRes(i)/sSize; % Amount to downsample by
    % Randomized starting frame for downsampling
    phaseValues = randi([0 dFactor-1],nVox,1); 
    
    % Build name for output file
    curFile = ['Downsample-Noise-' num2str(q) '-TRes-' num2str(i) '.mat'];
    
    % Pre-allocate
    [pkTM, pkRTM] = deal(zeros(nVox,2)); 
    [pkETM, pkERTM, pkLRRM, pkLRRM_raw, pkCLRRM] = ...
        deal(zeros(nVox,3));
    [pkERRM, pkCERRM] = deal(zeros(nVox,5));
    
    % Downsample and fit Extended / Tofts Model to all voxels
    parfor j=1:nVox
        curT = downsample(t, dFactor, phaseValues(j));
        curCt = downsample(Ct(:,j), dFactor, phaseValues(j));
        curCp = downsample(Cp(:,j), dFactor, phaseValues(j));
        pkTM(j,:)=LLSQ(curCt,curCp,curT,0);
        pkETM(j,:)=LLSQ(curCt,curCp,curT,1);
    end
    outFile = fullfile(outDir, 'ETM', curFile);
    save(outFile, 'pkTM','pkETM','phaseValues');
    
    % Do ERRM, RRM, RTM, ERTM
    parfor j=1:nVox
        curT = downsample(t, dFactor, phaseValues(j));
        curCt = downsample(Ct(:,j), dFactor, phaseValues(j));
        curCrr = downsample(Crr(:,j), dFactor, phaseValues(j));
        pkERRM(j,:)= ERRM(curCt,curCrr,curT,0,0,0);
        [pkLRRM(j,:), ~, pkLRRM_raw(j,:)]=LRRM(curCt,curCrr,curT);
        pkRTM(j,:)= doRTM(curCt,curCrr,curT,[KtRR veRR],0);
        pkERTM(j,:)= doRTM(curCt,curCrr,curT,[KtRR veRR],1);
    end
    
    % Get estimate for kepRR - from ERRM
    rawKepRR = pkERRM(:,5); % Estimated kepRR from all fits
    if std(rawKepRR)<1e-3
       % The proposed approach uses the interquartile mean to estimate kepRR
       % but this doesn't work in noiseless case if interquartile ~ 0.
       % So this 'if statement' is to handle the noiseless case. 
       % Here, the median is used instead of interquartile mean
       estKepRR = nanmedian(rawKepRR);
    else
       % Find voxels where all estimates are positive and kepRR estimate is real
       goodVals = pkERRM(:,1)>0 & pkERRM(:,2)>0 & pkERRM(:,3)>0 & pkERRM(:,4)>0 & pkERRM(:,5)>0 & imag(pkERRM(:,5))==0;
       estKepRR = iqrMean(rawKepRR(goodVals)); % Take the interquartile mean
    end
    
    % Get estimate for kepRR - from LRRM
    % This is used by the non-extended CLRRM (not mentioned in manuscript)
    p1 = pkLRRM_raw(:,1); % First fitted parameter from LRRM
    p2 = pkLRRM_raw(:,2); % Second fitted parameter from LRRM
    goodVals = (p1>0) & (p2>0) & pkLRRM_raw(:,3)>0; % Good values are when all parameters are positive
    % Get kepRR, which is ratio of p2/p1, for voxels where all fitted parameters are positive
    x= p2(goodVals)./p1(goodVals);
    % Using the median, as detailed in previous paper, doi: 10.1002/mrm.26530
    estKepRR_LRRM = nanmedian(x(x>0));

    % Do CERRM and CLRRM
    parfor j=1:nVox
        curT = downsample(t, dFactor, phaseValues(j));
        curCt = downsample(Ct(:,j), dFactor, phaseValues(j));
        curCrr = downsample(Crr(:,j), dFactor, phaseValues(j));
        pkCERRM(j,:) = CERRM(curCt,curCrr,curT,estKepRR);
        pkCLRRM(j,:) = CLRRM(curCt,curCrr,curT,estKepRR_LRRM);
    end
    outFile = fullfile(outDir, 'CERRM', curFile);
    save(outFile, 'pkCERRM','pkERRM','pkLRRM','pkCLRRM','pkRTM','pkERTM',...
        'estKepRR','estKepRR_LRRM', 'sigmaC', 'phaseValues');
    
    if doNonNeg
        % Re-do everything using non-negative least squares
        
        % Pre-allocate - don't know why I didn't use deal()
        pkTM = zeros(nVox,2);
        pkETM = zeros(nVox,3);
        pkRTM = zeros(nVox,2);
        pkERTM = zeros(nVox,3);
        pkERRM = zeros(nVox,5);
        pkLRRM = zeros(nVox,3);
        pkLRRM_raw = zeros(nVox,3);
        pkCERRM = zeros(nVox,5);
        pkCLRRM = zeros(nVox,3);

        % Downsample and fit Tofts Model and Extended Tofts Model
        parfor j=1:nVox
            curT = downsample(t, dFactor,phaseValues(j));
            curCt = downsample(Ct(:,j), dFactor,phaseValues(j));
            curCp = downsample(Cp(:,j), dFactor, phaseValues(j));
            pkTM(j,:)=LLSQ(curCt,curCp,curT,0,doNonNeg);
            pkETM(j,:)=LLSQ(curCt,curCp,curT,1,doNonNeg);
        end
        outFile = fullfile(outDirNN, 'ETM', curFile);
        save(outFile, 'pkTM','pkETM','phaseValues');

        % Do ERRM, RRM, RTM, ERTM
        parfor j=1:nVox
            curT = downsample(t, dFactor, phaseValues(j));
            curCt = downsample(Ct(:,j), dFactor, phaseValues(j));
            curCrr = downsample(Crr(:,j), dFactor, phaseValues(j));
            [pkLRRM(j,:), ~, pkLRRM_raw(j,:)]=LRRM(curCt,curCrr,curT);
            pkERRM(j,:)=ERRM(curCt,curCrr,curT,doNonNeg,0,0);
            pkRTM(j,:)=doRTM(curCt,curCrr,curT,[KtRR veRR],0,doNonNeg);
            pkERTM(j,:)=doRTM(curCt,curCrr,curT,[KtRR veRR],1,doNonNeg);
        end
        
        % Get estimate for kepRR - from ERRM
        rawKepRR = pkERRM(:,5);
        if std(rawKepRR)<1e-3
           estKepRR = nanmedian(rawKepRR);
        else
           % Find voxels where all estimates are real and positive
           goodVals = pkERRM(:,1)>0 & pkERRM(:,2)>0 & pkERRM(:,3)>0 & pkERRM(:,4)>0 & pkERRM(:,5)>0 & imag(pkERRM(:,5))==0;
           estKepRR = iqrMean(rawKepRR(goodVals));
        end
        
        % Get estimate for kepRR - from LRRM
        p1 = pkLRRM_raw(:,1);
        p2 = pkLRRM_raw(:,2);
        goodVals = (p1>0) & (p2>0) & pkLRRM_raw(:,3)>0;
        x= p2(goodVals)./p1(goodVals);
        estKepRR_LRRM = nanmedian(x(x>0));

        % Do CERRM and CLRRM
        parfor j=1:nVox
            curT = downsample(t, dFactor, phaseValues(j));
            curCt = downsample(Ct(:,j), dFactor, phaseValues(j));
            curCrr = downsample(Crr(:,j), dFactor, phaseValues(j));
            pkCERRM(j,:) = CERRM(curCt,curCrr,curT,estKepRR,doNonNeg);
            pkCLRRM(j,:) = CLRRM(curCt,curCrr,curT,estKepRR_LRRM);
        end
        outFile = fullfile(outDirNN, 'CERRM', curFile);
        save(outFile, 'pkCERRM','pkERRM','pkRTM','pkERTM',...
            'estKepRR','sigmaC', 'phaseValues'); 
    end
end % for TRes
end % for sigmaC
toc
disp('Finished fitting data')
