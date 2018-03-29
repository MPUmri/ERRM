% Collects results from previous step and exports summary as csv

% Estimated runtime: 60 seconds

%% Pre-setup

clearvars
addpath('./mfiles')

%% Setup

% Was non-negative least-squares performed?
didNonNeg = false;

% Output directories
outDir = './data/mapDownsample';
outDirNN = './data/mapDownsample-NN';

%%
tic

load('./data/simMap.mat');

% The settings below show be consistent with e03_2_analyzeSimMap.m
TRes = 1:60; % Range of temporal resolutions that were simulated
listSigmaC = 0:0.01:0.05;  % Range of noise levels

% Reference region parameters
ktRR = 0.07;
kepRR = 0.5;
veRR = ktRR/kepRR;

repF = 100; % Number of replications for each noise level
trueKt = repmat(trueKt(:),[repF 1]);
trueKep = repmat(trueKep(:),[repF 1]);
trueVe = repmat(trueVe(:),[repF 1]);
trueVp = repmat(trueVp(:),[repF 1]);

%% Vp cutoff - Not actually used
% A range of vp can be defined here so that statistics are only calculated
% for curves with specific vp values
% Settings to [-inf inf] will use all simulated vp values
vpCutoff = [-inf inf]; % Min/Max true vp

vpMask = trueVp(:)>=vpCutoff(1) & trueVp(:)<=vpCutoff(2);
trueKt(~vpMask)=[];
trueKep(~vpMask)=[];
trueVe(~vpMask)=[];
trueVp(~vpMask)=[];

%% Collect results
for j=1:length(listSigmaC)
for i=1:length(TRes)
    curFile = ['Downsample-Noise-' num2str(j) '-TRes-' num2str(i) '.mat'];
    load(fullfile(outDir, 'CERRM', curFile));
    load(fullfile(outDir, 'ETM', curFile));
    
    % Record number of imaginary estimates from ERRM
    numImagE(i,j) = sum( sum(abs(imag(pkERRM')))>0 );
    
    % Use only real part when calculating percent error for ERRM
    pkERRM = real(pkERRM);
    
    % OPTIONAL
    % For ERRM: we could exclude fits where tissue of interest and 
    % reference tissue have same kep values, as the ERRM fails to fit those cases
    goodERRM = (trueKep~=kepRR);
    % However, we will not remove those fits, so that the same number of
    % fits are used for all models
    % (It doesn't really change the results anyways)
    goodERRM(:) = true; % Comment this line to ignore voxels where kep=kepRR
    
    %% Apply vpCutOff
    % Does nothing if vpCutoff = [-inf inf]  (Default setting)
    pkCERRM(~vpMask,:)=[];
    pkCLRRM(~vpMask,:)=[];
    pkERRM(~vpMask,:)=[];
    pkERTM(~vpMask,:)=[];
    pkETM(~vpMask,:)=[];
    pkLRRM(~vpMask,:)=[];
    pkRTM(~vpMask,:)=[];
    pkTM(~vpMask,:)=[];
    %% Collect mean, stdDev, and quantiles of the percent errors
    
    % Tofts Model - TM
    [ktError, meanKtErrT(i,j), stdKtErrT(i,j)] = PercentError(pkTM(:,1),trueKt(:));
    [veError, meanVeErrT(i,j), stdVeErrT(i,j)] = PercentError(pkTM(:,1)./pkTM(:,2),trueVe(:));
    [kepError, meanKepErrT(i,j), stdKepErrT(i,j)] = PercentError(pkTM(:,2),trueKep(:));
    
    qKtT(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
    qVeT(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
    qKepT(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);
    
    % Extended Tofts Model - ETM
    [ktError, meanKtErrET(i,j), stdKtErrET(i,j)] = PercentError(pkETM(:,1),trueKt(:));
    [veError, meanVeErrET(i,j), stdVeErrET(i,j)] = PercentError(pkETM(:,1)./pkETM(:,2),trueVe(:));
    [kepError, meanKepErrET(i,j), stdKepErrET(i,j)] = PercentError(pkETM(:,2),trueKep(:));
    [vpError, meanVpErrET(i,j), stdVpErrET(i,j)] = PercentError(pkETM(:,3),trueVp(:));
    
    qKtET(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
    qVeET(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
    qKepET(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);
    qVpET(i,j,:) = quantile(vpError,[.05 .25 .5 .75 .95]);
    
    % Reference Tissue Method - RTM (not mentioned in manuscript)
    [ktError, meanKtErrRT(i,j), stdKtErrRT(i,j)] = PercentError(pkRTM(:,1),trueKt(:));
    [veError, meanVeErrRT(i,j), stdVeErrRT(i,j)] = PercentError(pkRTM(:,1)./pkRTM(:,2),trueVe(:));
    [kepError, meanKepErrRT(i,j), stdKepErrRT(i,j)] = PercentError(pkRTM(:,2),trueKep(:));
    
    qKtRT(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
    qVeRT(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
    qKepRT(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);
    
    % Extended Reference Tissue Method - ERTM (not mentioned in manuscript)
    [ktError, meanKtErrERT(i,j), stdKtErrERT(i,j)] = PercentError(pkERTM(:,1),trueKt(:));
    [veError, meanVeErrERT(i,j), stdVeErrERT(i,j)] = PercentError(pkERTM(:,1)./pkERTM(:,2),trueVe(:));
    [kepError, meanKepErrERT(i,j), stdKepErrERT(i,j)] = PercentError(pkERTM(:,2),trueKep(:));
    [vpError, meanVpErrERT(i,j), stdVpErrERT(i,j)] = PercentError(pkERTM(:,3),trueVp(:));
    
    qKtERT(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
    qVeERT(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
    qKepERT(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);
    qVpERT(i,j,:) = quantile(vpError,[.05 .25 .5 .75 .95]);
    
    % Constrained Extended Reference Region Model - CERRM
    [ktError, meanKtErrCE(i,j), stdKtErrCE(i,j), badValKtCE(i,j)] = PercentError(ktRR*pkCERRM(:,1),trueKt(:));
    [veError, meanVeErrCE(i,j), stdVeErrCE(i,j), badValVeCE(i,j)] = PercentError(veRR*pkCERRM(:,2),trueVe(:));
    [kepError, meanKepErrCE(i,j), stdKepErrCE(i,j), badValKepCE(i,j)] = PercentError(pkCERRM(:,3),trueKep(:));
    [vpError, meanVpErrCE(i,j), stdVpErrCE(i,j), badValVpCE(i,j)] = PercentError(ktRR*pkCERRM(:,4),trueVp(:));
    
    qKtCE(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
    qVeCE(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
    qKepCE(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);
    qVpCE(i,j,:) = quantile(vpError,[.05 .25 .5 .75 .95]);
    
    % Extended Reference Region Model - ERRM
    [ktError, meanKtErrE(i,j), stdKtErrE(i,j), badValKtE(i,j)] = PercentError(ktRR*pkERRM(:,1),trueKt(:));
    [veError, meanVeErrE(i,j), stdVeErrE(i,j), badValVeE(i,j)] = PercentError(veRR*pkERRM(:,2),trueVe(:));
    [kepError, meanKepErrE(i,j), stdKepErrE(i,j), badValKepE(i,j)] = PercentError(pkERRM(:,3),trueKep(:));
    [vpError, meanVpErrE(i,j), stdVpErrE(i,j), badValVpE(i,j)] = PercentError(ktRR*pkERRM(:,4),trueVp(:));
    
    qKtE(i,j,:) = quantile(ktError(goodERRM),[.05 .25 .5 .75 .95]);
    qVeE(i,j,:) = quantile(veError(goodERRM),[.05 .25 .5 .75 .95]);
    qKepE(i,j,:) = quantile(kepError(goodERRM),[.05 .25 .5 .75 .95]);
    qVpE(i,j,:) = quantile(vpError,[.05 .25 .5 .75 .95]);
    
    % Reference Region Model - RRM
    [ktError, meanKtErrR(i,j), stdKtErrR(i,j)] = PercentError(ktRR*pkLRRM(:,1),trueKt(:));
    [veError, meanVeErrR(i,j), stdVeErrR(i,j)] = PercentError(veRR*pkLRRM(:,2),trueVe(:));
    [kepError, meanKepErrR(i,j), stdKepErrR(i,j)] = PercentError(pkLRRM(:,3),trueKep(:));
    
    qKtR(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
    qVeR(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
    qKepR(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);
    
    % Non-extended Contrained Reference Region Model - CRRM (not mentioned in manuscript)
    [ktError, meanKtErrCR(i,j), stdKtErrCR(i,j)] = PercentError(ktRR*pkCLRRM(:,1),trueKt(:));
    [veError, meanVeErrCR(i,j), stdVeErrCR(i,j)] = PercentError(veRR*pkCLRRM(:,2),trueVe(:));
    [kepError, meanKepErrCR(i,j), stdKepErrCR(i,j)] = PercentError(pkCLRRM(:,3),trueKep(:));
    
    qKtCR(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
    qVeCR(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
    qKepCR(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);
    
    %% Load fits using NNLS - if the fits were performed
    % This is pretty much a copy/paste job of the code above, with minor
    % changes. Minimal comments ahead.
    if didNonNeg
        load(fullfile(outDirNN, 'CERRM', curFile));
        load(fullfile(outDirNN, 'ETM', curFile));
        
        numImagEnn(i,j) = sum( sum(abs(imag(pkERRM')))>0 );
        pkERRM = real(pkERRM);
        
        % Apply vpCutOff
        pkCERRM(~vpMask,:)=[];
        pkERRM(~vpMask,:)=[];
        pkERTM(~vpMask,:)=[];
        pkETM(~vpMask,:)=[];
        pkRTM(~vpMask,:)=[];
        pkTM(~vpMask,:)=[];
        
        % TM
        [ktError, meanKtErrTnn(i,j), stdKtErrTnn(i,j)] = PercentError(pkTM(:,1),trueKt(:));
        [veError, meanVeErrTnn(i,j), stdVeErrTnn(i,j)] = PercentError(pkTM(:,1)./pkTM(:,2),trueVe(:));
        [kepError, meanKepErrTnn(i,j), stdKepErrTnn(i,j)] = PercentError(pkTM(:,2),trueKep(:));

        qKtTnn(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
        qVeTnn(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
        qKepTnn(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);

        % ETM
        [ktError, meanKtErrETnn(i,j), stdKtErrETnn(i,j)] = PercentError(pkETM(:,1),trueKt(:));
        [veError, meanVeErrETnn(i,j), stdVeErrETnn(i,j)] = PercentError(pkETM(:,1)./pkETM(:,2),trueVe(:));
        [kepError, meanKepErrETnn(i,j), stdKepErrETnn(i,j)] = PercentError(pkETM(:,2),trueKep(:));
        [vpError, meanVpErrETnn(i,j), stdVpErrETnn(i,j)] = PercentError(pkETM(:,3),trueVp(:));

        qKtETnn(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
        qVeETnn(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
        qKepETnn(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);
        qVpETnn(i,j,:) = quantile(vpError,[.05 .25 .5 .75 .95]);

        % RTM
        [ktError, meanKtErrRTnn(i,j), stdKtErrRTnn(i,j)] = PercentError(pkRTM(:,1),trueKt(:));
        [veError, meanVeErrRTnn(i,j), stdVeErrRTnn(i,j)] = PercentError(pkRTM(:,1)./pkRTM(:,2),trueVe(:));
        [kepError, meanKepErrRTnn(i,j), stdKepErrRTnn(i,j)] = PercentError(pkRTM(:,2),trueKep(:));

        qKtRTnn(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
        qVeRTnn(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
        qKepRTnn(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);

        % ERTM
        [ktError, meanKtErrERTnn(i,j), stdKtErrERTnn(i,j)] = PercentError(pkERTM(:,1),trueKt(:));
        [veError, meanVeErrERTnn(i,j), stdVeErrERTnn(i,j)] = PercentError(pkERTM(:,1)./pkERTM(:,2),trueVe(:));
        [kepError, meanKepErrERTnn(i,j), stdKepErrERTnn(i,j)] = PercentError(pkERTM(:,2),trueKep(:));
        [vpError, meanVpErrERTnn(i,j), stdVpErrERTnn(i,j)] = PercentError(pkERTM(:,3),trueVp(:));

        qKtERTnn(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
        qVeERTnn(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
        qKepERTnn(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);
        qVpERTnn(i,j,:) = quantile(vpError,[.05 .25 .5 .75 .95]);

        % CERRM
        [ktError, meanKtErrCEnn(i,j), stdKtErrCEnn(i,j), badValKtCEnn(i,j)] = PercentError(ktRR*pkCERRM(:,1),trueKt(:));
        [veError, meanVeErrCEnn(i,j), stdVeErrCEnn(i,j), badValVeCEnn(i,j)] = PercentError(veRR*pkCERRM(:,2),trueVe(:));
        [kepError, meanKepErrCEnn(i,j), stdKepErrCEnn(i,j), badValKepCEnn(i,j)] = PercentError(pkCERRM(:,3),trueKep(:));
        [vpError, meanVpErrCEnn(i,j), stdVpErrCEnn(i,j), badValVpCEnn(i,j)] = PercentError(ktRR*pkCERRM(:,4),trueVp(:));

        qKtCEnn(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
        qVeCEnn(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
        qKepCEnn(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);
        qVpCEnn(i,j,:) = quantile(vpError,[.05 .25 .5 .75 .95]);

        % ERRM
        [ktError, meanKtErrEnn(i,j), stdKtErrEnn(i,j), badValKtEnn(i,j)] = PercentError(ktRR*pkERRM(:,1),trueKt(:));
        [veError, meanVeErrEnn(i,j), stdVeErrEnn(i,j), badValVeEnn(i,j)] = PercentError(veRR*pkERRM(:,2),trueVe(:));
        [kepError, meanKepErrEnn(i,j), stdKepErrEnn(i,j), badValKepEnn(i,j)] = PercentError(pkERRM(:,3),trueKep(:));
        [vpError, meanVpErrEnn(i,j), stdVpErrEnn(i,j), badValVpEnn(i,j)] = PercentError(ktRR*pkERRM(:,4),trueVp(:));

        qKtEnn(i,j,:) = quantile(ktError(goodERRM),[.05 .25 .5 .75 .95]);
        qVeEnn(i,j,:) = quantile(veError(goodERRM),[.05 .25 .5 .75 .95]);
        qKepEnn(i,j,:) = quantile(kepError(goodERRM),[.05 .25 .5 .75 .95]);
        qVpEnn(i,j,:) = quantile(vpError,[.05 .25 .5 .75 .95]);
    
    end
end
end

%% Export mean and std of percent errors as csv
outFile = fullfile('./dataResults','e03a-downsampleMapResultsMean.csv');
hdr='FitMethod,TemporalRes,sigmaC,errKt,errVe,errKep,errVp,stdKt,stdVe,stdKep,stdVp';
outID = fopen(outFile, 'w+');
fprintf(outID, '%s\n', hdr); % Print header into csv file
for j=1:length(listSigmaC)
for i=1:length(TRes)
   outLine = {'ETM',TRes(i),listSigmaC(j),meanKtErrET(i,j),meanVeErrET(i,j),meanKepErrET(i,j),meanVpErrET(i,j),...
       stdKtErrET(i,j),stdVeErrET(i,j),stdKepErrET(i,j),stdVpErrET(i,j)};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'ERTM',TRes(i),listSigmaC(j),meanKtErrERT(i,j),meanVeErrERT(i,j),meanKepErrERT(i,j),meanVpErrERT(i,j),...
       stdKtErrERT(i,j),stdVeErrERT(i,j),stdKepErrERT(i,j),stdVpErrERT(i,j)};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'RTM',TRes(i),listSigmaC(j),meanKtErrRT(i,j),meanVeErrRT(i,j),meanKepErrRT(i,j),NaN,...
       stdKtErrRT(i,j),stdVeErrRT(i,j),stdKepErrRT(i,j),NaN};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'ERRM',TRes(i),listSigmaC(j),meanKtErrE(i,j),meanVeErrE(i,j),meanKepErrE(i,j),meanVpErrE(i,j),...
       stdKtErrE(i,j),stdVeErrE(i,j),stdKepErrE(i,j),stdVpErrE(i,j)};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'CERRM',TRes(i),listSigmaC(j),meanKtErrCE(i,j),meanVeErrCE(i,j),meanKepErrCE(i,j),meanVpErrCE(i,j),...
       stdKtErrCE(i,j),stdVeErrCE(i,j),stdKepErrCE(i,j),stdVpErrCE(i,j)};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'TM',TRes(i),listSigmaC(j),meanKtErrT(i,j),meanVeErrT(i,j),meanKepErrT(i,j),NaN,...
       stdKtErrT(i,j),stdVeErrT(i,j),stdKepErrT(i,j),NaN};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'RRM',TRes(i),listSigmaC(j),meanKtErrR(i,j),meanVeErrR(i,j),meanKepErrR(i,j),NaN,...
       stdKtErrR(i,j),stdVeErrR(i,j),stdKepErrR(i,j),NaN};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'CRRM',TRes(i),listSigmaC(j),meanKtErrCR(i,j),meanVeErrCR(i,j),meanKepErrCR(i,j),NaN,...
       stdKtErrCR(i,j),stdVeErrCR(i,j),stdKepErrCR(i,j),NaN};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   % Export results from NNLS fits, if they were performed
   if ~didNonNeg
       continue
   end
   
   outLine = {'ETM-NN',TRes(i),listSigmaC(j),meanKtErrETnn(i,j),meanVeErrETnn(i,j),meanKepErrETnn(i,j),meanVpErrETnn(i,j),...
       stdKtErrETnn(i,j),stdVeErrETnn(i,j),stdKepErrETnn(i,j),stdVpErrETnn(i,j)};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'ERTM-NN',TRes(i),listSigmaC(j),meanKtErrERTnn(i,j),meanVeErrERTnn(i,j),meanKepErrERTnn(i,j),meanVpErrERTnn(i,j),...
       stdKtErrERTnn(i,j),stdVeErrERTnn(i,j),stdKepErrERTnn(i,j),stdVpErrERTnn(i,j)};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'RTM-NN',TRes(i),listSigmaC(j),meanKtErrRTnn(i,j),meanVeErrRTnn(i,j),meanKepErrRTnn(i,j),NaN,...
       stdKtErrRTnn(i,j),stdVeErrRTnn(i,j),stdKepErrRTnn(i,j),NaN};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'ERRM-NN',TRes(i),listSigmaC(j),meanKtErrEnn(i,j),meanVeErrEnn(i,j),meanKepErrEnn(i,j),meanVpErrEnn(i,j),...
       stdKtErrEnn(i,j),stdVeErrEnn(i,j),stdKepErrEnn(i,j),stdVpErrEnn(i,j)};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'CERRM-NN',TRes(i),listSigmaC(j),meanKtErrCEnn(i,j),meanVeErrCEnn(i,j),meanKepErrCEnn(i,j),meanVpErrCEnn(i,j),...
       stdKtErrCEnn(i,j),stdVeErrCEnn(i,j),stdKepErrCEnn(i,j),stdVpErrCEnn(i,j)};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'TM-NN',TRes(i),listSigmaC(j),meanKtErrTnn(i,j),meanVeErrTnn(i,j),meanKepErrTnn(i,j),NaN,...
       stdKtErrTnn(i,j),stdVeErrTnn(i,j),stdKepErrTnn(i,j),NaN};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
end
end
fclose(outID);

%% Export quartiles of percent error as csv
outFile = fullfile('./dataResults',['e03b-downsampleMapResultsQuantiles.csv']);
hdr=['FitMethod,TemporalRes,sigmaC,q5Kt,q25Kt,q50Kt,q75Kt,q95Kt,q5Ve,q25Ve,q50Ve,q75Ve,q95Ve,q5Kep,q25Kep,q50Kep,q75Kep,q95Kep,q5Vp,q25Vp,q50Vp,q75Vp,q95Vp'];
outID = fopen(outFile, 'w+');
fprintf(outID, '%s\n', hdr); % Print header into csv file
for j=1:length(listSigmaC)
for i=1:length(TRes)
   outLine = {'ETM',TRes(i),listSigmaC(j),...
       qKtET(i,j,1),qKtET(i,j,2),qKtET(i,j,3),qKtET(i,j,4),qKtET(i,j,5),...
       qVeET(i,j,1),qVeET(i,j,2),qVeET(i,j,3),qVeET(i,j,4),qVeET(i,j,5),...
       qKepET(i,j,1),qKepET(i,j,2),qKepET(i,j,3),qKepET(i,j,4),qKepET(i,j,5),...
       qVpET(i,j,1),qVpET(i,j,2),qVpET(i,j,3),qVpET(i,j,4),qVpET(i,j,5),...
       };
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'ERTM',TRes(i),listSigmaC(j),...
       qKtERT(i,j,1),qKtERT(i,j,2),qKtERT(i,j,3),qKtERT(i,j,4),qKtERT(i,j,5),...
       qVeERT(i,j,1),qVeERT(i,j,2),qVeERT(i,j,3),qVeERT(i,j,4),qVeERT(i,j,5),...
       qKepERT(i,j,1),qKepERT(i,j,2),qKepERT(i,j,3),qKepERT(i,j,4),qKepERT(i,j,5),...
       qVpERT(i,j,1),qVpERT(i,j,2),qVpERT(i,j,3),qVpERT(i,j,4),qVpERT(i,j,5),...
       };
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'RTM',TRes(i),listSigmaC(j),...
       qKtRT(i,j,1),qKtRT(i,j,2),qKtRT(i,j,3),qKtRT(i,j,4),qKtRT(i,j,5),...
       qVeRT(i,j,1),qVeRT(i,j,2),qVeRT(i,j,3),qVeRT(i,j,4),qVeRT(i,j,5),...
       qKepRT(i,j,1),qKepRT(i,j,2),qKepRT(i,j,3),qKepRT(i,j,4),qKepRT(i,j,5),...
       NaN,NaN,NaN,NaN,NaN...
       };
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'CERRM',TRes(i),listSigmaC(j),...
       qKtCE(i,j,1),qKtCE(i,j,2),qKtCE(i,j,3),qKtCE(i,j,4),qKtCE(i,j,5),...
       qVeCE(i,j,1),qVeCE(i,j,2),qVeCE(i,j,3),qVeCE(i,j,4),qVeCE(i,j,5),...
       qKepCE(i,j,1),qKepCE(i,j,2),qKepCE(i,j,3),qKepCE(i,j,4),qKepCE(i,j,5),...
       qVpCE(i,j,1),qVpCE(i,j,2),qVpCE(i,j,3),qVpCE(i,j,4),qVpCE(i,j,5),...
       };
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'ERRM',TRes(i),listSigmaC(j),...
       qKtE(i,j,1),qKtE(i,j,2),qKtE(i,j,3),qKtE(i,j,4),qKtE(i,j,5),...
       qVeE(i,j,1),qVeE(i,j,2),qVeE(i,j,3),qVeE(i,j,4),qVeE(i,j,5),...
       qKepE(i,j,1),qKepE(i,j,2),qKepE(i,j,3),qKepE(i,j,4),qKepE(i,j,5),...
       qVpE(i,j,1),qVpE(i,j,2),qVpE(i,j,3),qVpE(i,j,4),qVpE(i,j,5),...
       };
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'TM',TRes(i),listSigmaC(j),...
       qKtT(i,j,1),qKtT(i,j,2),qKtT(i,j,3),qKtT(i,j,4),qKtT(i,j,5),...
       qVeT(i,j,1),qVeT(i,j,2),qVeT(i,j,3),qVeT(i,j,4),qVeT(i,j,5),...
       qKepT(i,j,1),qKepT(i,j,2),qKepT(i,j,3),qKepT(i,j,4),qKepT(i,j,5),...
       NaN,NaN,NaN,NaN,NaN,...
       };
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'RRM',TRes(i),listSigmaC(j),...
       qKtR(i,j,1),qKtR(i,j,2),qKtR(i,j,3),qKtR(i,j,4),qKtR(i,j,5),...
       qVeR(i,j,1),qVeR(i,j,2),qVeR(i,j,3),qVeR(i,j,4),qVeR(i,j,5),...
       qKepR(i,j,1),qKepR(i,j,2),qKepR(i,j,3),qKepR(i,j,4),qKepR(i,j,5),...
       NaN,NaN,NaN,NaN,NaN,...
       };
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'CRRM',TRes(i),listSigmaC(j),...
       qKtCR(i,j,1),qKtCR(i,j,2),qKtCR(i,j,3),qKtCR(i,j,4),qKtCR(i,j,5),...
       qVeCR(i,j,1),qVeCR(i,j,2),qVeCR(i,j,3),qVeCR(i,j,4),qVeCR(i,j,5),...
       qKepCR(i,j,1),qKepCR(i,j,2),qKepCR(i,j,3),qKepCR(i,j,4),qKepCR(i,j,5),...
       NaN,NaN,NaN,NaN,NaN,...
       };
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:}); 
   
   % Export results from NNLS fits, if they were performed
   if ~didNonNeg
       continue
   end
   
   outLine = {'ETM-NN',TRes(i),listSigmaC(j),...
       qKtETnn(i,j,1),qKtETnn(i,j,2),qKtETnn(i,j,3),qKtETnn(i,j,4),qKtETnn(i,j,5),...
       qVeETnn(i,j,1),qVeETnn(i,j,2),qVeETnn(i,j,3),qVeETnn(i,j,4),qVeETnn(i,j,5),...
       qKepETnn(i,j,1),qKepETnn(i,j,2),qKepETnn(i,j,3),qKepETnn(i,j,4),qKepETnn(i,j,5),...
       qVpETnn(i,j,1),qVpETnn(i,j,2),qVpETnn(i,j,3),qVpETnn(i,j,4),qVpETnn(i,j,5),...
       };
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'ERTM-NN',TRes(i),listSigmaC(j),...
       qKtERTnn(i,j,1),qKtERTnn(i,j,2),qKtERTnn(i,j,3),qKtERTnn(i,j,4),qKtERTnn(i,j,5),...
       qVeERTnn(i,j,1),qVeERTnn(i,j,2),qVeERTnn(i,j,3),qVeERTnn(i,j,4),qVeERTnn(i,j,5),...
       qKepERTnn(i,j,1),qKepERTnn(i,j,2),qKepERTnn(i,j,3),qKepERTnn(i,j,4),qKepERTnn(i,j,5),...
       qVpERTnn(i,j,1),qVpERTnn(i,j,2),qVpERTnn(i,j,3),qVpERTnn(i,j,4),qVpERTnn(i,j,5),...
       };
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'RTM-NN',TRes(i),listSigmaC(j),...
       qKtRTnn(i,j,1),qKtRTnn(i,j,2),qKtRTnn(i,j,3),qKtRTnn(i,j,4),qKtRTnn(i,j,5),...
       qVeRTnn(i,j,1),qVeRTnn(i,j,2),qVeRTnn(i,j,3),qVeRTnn(i,j,4),qVeRTnn(i,j,5),...
       qKepRTnn(i,j,1),qKepRTnn(i,j,2),qKepRTnn(i,j,3),qKepRTnn(i,j,4),qKepRTnn(i,j,5),...
       NaN,NaN,NaN,NaN,NaN...
       };
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'CERRM-NN',TRes(i),listSigmaC(j),...
       qKtCEnn(i,j,1),qKtCEnn(i,j,2),qKtCEnn(i,j,3),qKtCEnn(i,j,4),qKtCEnn(i,j,5),...
       qVeCEnn(i,j,1),qVeCEnn(i,j,2),qVeCEnn(i,j,3),qVeCEnn(i,j,4),qVeCEnn(i,j,5),...
       qKepCEnn(i,j,1),qKepCEnn(i,j,2),qKepCEnn(i,j,3),qKepCEnn(i,j,4),qKepCEnn(i,j,5),...
       qVpCEnn(i,j,1),qVpCEnn(i,j,2),qVpCEnn(i,j,3),qVpCEnn(i,j,4),qVpCEnn(i,j,5),...
       };
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'ERRM-NN',TRes(i),listSigmaC(j),...
       qKtEnn(i,j,1),qKtEnn(i,j,2),qKtEnn(i,j,3),qKtEnn(i,j,4),qKtEnn(i,j,5),...
       qVeEnn(i,j,1),qVeEnn(i,j,2),qVeEnn(i,j,3),qVeEnn(i,j,4),qVeEnn(i,j,5),...
       qKepEnn(i,j,1),qKepEnn(i,j,2),qKepEnn(i,j,3),qKepEnn(i,j,4),qKepEnn(i,j,5),...
       qVpEnn(i,j,1),qVpEnn(i,j,2),qVpEnn(i,j,3),qVpEnn(i,j,4),qVpEnn(i,j,5),...
       };
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'TM-NN',TRes(i),listSigmaC(j),...
       qKtTnn(i,j,1),qKtTnn(i,j,2),qKtTnn(i,j,3),qKtTnn(i,j,4),qKtTnn(i,j,5),...
       qVeTnn(i,j,1),qVeTnn(i,j,2),qVeTnn(i,j,3),qVeTnn(i,j,4),qVeTnn(i,j,5),...
       qKepTnn(i,j,1),qKepTnn(i,j,2),qKepTnn(i,j,3),qKepTnn(i,j,4),qKepTnn(i,j,5),...
       NaN,NaN,NaN,NaN,NaN,...
       };
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
end
end
fclose(outID);
toc

%% Plots

%% Median +/- interquartile range vs StdDev of Noise

myX = listSigmaC;
myTres = 5;
lErr = 2;
hErr = 4;

clrERRM = [247,151,86]./255;
clrCERRM = [126,47,142]./255;

figure('Position',[300 300 1200 500])

subplot(1,4,1)
errorbar(myX,qKtE(myTres,:,3),abs(qKtE(myTres,:,lErr)-qKtE(myTres,:,3)),abs(qKtE(myTres,:,hErr)-qKtE(myTres,:,3)),'linewidth',3,'linestyle','--','Color',clrERRM)
hold on
errorbar(myX,qKtCE(myTres,:,3),abs(qKtCE(myTres,:,lErr)-qKtCE(myTres,:,3)),abs(qKtCE(myTres,:,hErr)-qKtCE(myTres,:,3)),'linewidth',3,'Color',clrCERRM)
hold off
ylim([-70 10])
xlabel('StdDev [mM]')
ylabel('Percent Error')
title('Ktrans')

subplot(1,4,2)
errorbar(myX,qKepE(myTres,:,3),abs(qKepE(myTres,:,lErr)-qKepE(myTres,:,3)),abs(qKepE(myTres,:,hErr)-qKepE(myTres,:,3)),'linewidth',3,'linestyle','--','Color',clrERRM)
hold on
errorbar(myX,qKepCE(myTres,:,3),abs(qKepCE(myTres,:,lErr)-qKepCE(myTres,:,3)),abs(qKepCE(myTres,:,hErr)-qKepCE(myTres,:,3)),'linewidth',3,'Color',clrCERRM)
hold off
ylim([-50 20])
xlabel('StdDev [mM]')
ylabel('Percent Error')
title('kep')

subplot(1,4,3)
errorbar(myX,qVeE(myTres,:,3),abs(qVeE(myTres,:,lErr)-qVeE(myTres,:,3)),abs(qVeE(myTres,:,hErr)-qVeE(myTres,:,3)),'linewidth',3,'linestyle','--','Color',clrERRM)
hold on
errorbar(myX,qVeCE(myTres,:,3),abs(qVeCE(myTres,:,lErr)-qVeCE(myTres,:,3)),abs(qVeCE(myTres,:,hErr)-qVeCE(myTres,:,3)),'linewidth',3,'Color',clrCERRM)
hold off
ylim([-50 20])
xlabel('StdDev [mM]')
ylabel('Percent Error')
title('ve')

subplot(1,4,4)
errorbar(myX,qVpE(myTres,:,3),abs(qVpE(myTres,:,lErr)-qVpE(myTres,:,3)),abs(qVpE(myTres,:,hErr)-qVpE(myTres,:,3)),'linewidth',3,'linestyle','--','Color',clrERRM)
hold on
errorbar(myX,qVpCE(myTres,:,3),abs(qVpCE(myTres,:,lErr)-qVpCE(myTres,:,3)),abs(qVpCE(myTres,:,hErr)-qVpCE(myTres,:,3)),'linewidth',3,'Color',clrCERRM)
hold off
ylim([-40 160])
xlabel('StdDev [mM]')
ylabel('Percent Error')
title('vp')
legend('ERRM','CERRM')

%% Ratios of interquartile range between ERRM and CERRM
quartRatio_Kt=(qKtE(myTres,:,hErr)-qKtE(myTres,:,lErr))./(qKtCE(myTres,:,hErr)-qKtCE(myTres,:,lErr));
quartRatio_Kep=(qKepE(myTres,:,hErr)-qKepE(myTres,:,lErr))./(qKepCE(myTres,:,hErr)-qKepCE(myTres,:,lErr));
quartRatio_Ve=(qVeE(myTres,:,hErr)-qVeE(myTres,:,lErr))./(qVeCE(myTres,:,hErr)-qVeCE(myTres,:,lErr));
quartRatio_Vp=(qVpE(myTres,:,hErr)-qVpE(myTres,:,lErr))./(qVpCE(myTres,:,hErr)-qVpCE(myTres,:,lErr));

figure
boxplot([quartRatio_Kt; quartRatio_Kep; quartRatio_Ve; quartRatio_Vp]', ...
    'Labels',{'Ktrans','kep','ve','vp'})
title('Ratio of Interquartile Range between ERRM and CERRM')
ylim([0 4])

% The above figure shows the ratio of the interquartile range, i.e. the error
% bars in Figure 2, between the ERRM and CERRM. For example, a value of 2
% means that the ERRM's error bars are twice as large as the CERRM.
% The boxplot shows the spread of this ratio for the different noise
% levels.

% The manuscript undersells the improvement of the CERRM by saying that its
% error bars were smaller by a factor of 1.5, 3, 1.5, and ~1 for Ktrans,
% kep, ve, and vp, respectively. 
% The figure shows a median improvement by a factor of 2, 3.2, 2, and ~1,
% respectively.

%% Median +/- interquartile range vs Temporal Resolution
mySig = 2;
lErr = 2;
hErr = 4;

xRange = 1:30;

cCE = [.5 .2 .55];
cET = [.3 .7 .6];

figure('Position',[300 300 1200 500])
subplot(1,4,1)
hold on
shadedErrorPlot(TRes(xRange)',qKtET(xRange,mySig,3),qKtET(xRange,mySig,lErr),qKtET(xRange,mySig,hErr),cET);
shadedErrorPlot(TRes(xRange)',qKtCE(xRange,mySig,3),qKtCE(xRange,mySig,lErr),qKtCE(xRange,mySig,hErr),cCE);
hold off
ylim([-60 40])
xlabel('Temporal Resolution [s]')
ylabel('Percent Error')
title('Ktrans')
legend('','ETM','','CERRM')

subplot(1,4,2)
hold on
shadedErrorPlot(TRes(xRange)',qKepET(xRange,mySig,3),qKepET(xRange,mySig,lErr),qKepET(xRange,mySig,hErr),cET);
shadedErrorPlot(TRes(xRange)',qKepCE(xRange,mySig,3),qKepCE(xRange,mySig,lErr),qKepCE(xRange,mySig,hErr),cCE);
hold off
ylim([-60 40])
xlabel('Temporal Resolution [s]')
ylabel('Percent Error')
title('kep')

subplot(1,4,3)
hold on
shadedErrorPlot(TRes(xRange)',qVeET(xRange,mySig,3),qVeET(xRange,mySig,lErr),qVeET(xRange,mySig,hErr),cET);
shadedErrorPlot(TRes(xRange)',qVeCE(xRange,mySig,3),qVeCE(xRange,mySig,lErr),qVeCE(xRange,mySig,hErr),cCE);
hold off
ylim([-60 40])
xlabel('Temporal Resolution [s]')
ylabel('Percent Error')
title('ve')

subplot(1,4,4)
hold on
shadedErrorPlot(TRes(xRange)',qVpET(xRange,mySig,3),qVpET(xRange,mySig,lErr),qVpET(xRange,mySig,hErr),cET);
shadedErrorPlot(TRes(xRange)',qVpCE(xRange,mySig,3),qVpCE(xRange,mySig,lErr),qVpCE(xRange,mySig,hErr),cCE);
hold off
ylim([-100 300])
xlabel('Temporal Resolution [s]')
ylabel('Percent Error')
title('vp')

%% Make the Error Maps
% These steps should've been moved to a new script, but were appended here
% The next steps load a specific dataset (noiseless, TemporalRes = 1s)
% and makes the percent error maps for kep (Fig 1 in  manuscript)
% and computes the error for RRM at vp=0.01 (mentioned in manuscript text)

% Pick which dataset to load
indTRes = 1; % Default (1) : Temporal resolutions = 1 s
indSigma = 1; % Default (1) : Noiseless
repF = 100; % Default (100) | Replications for each noise

curFile = ['Downsample-Noise-' num2str(indSigma) '-TRes-' num2str(indTRes) '.mat'];
load(fullfile(outDir, 'ETM', curFile));
load(fullfile(outDir, 'CERRM', curFile));

% Apply vpCutOff - unnecessary
pkCERRM(~vpMask,:)=[];
pkERRM(~vpMask,:)=[];
pkLRRM(~vpMask,:)=[];
pkERTM(~vpMask,:)=[];
pkETM(~vpMask,:)=[];
pkRTM(~vpMask,:)=[];
pkTM(~vpMask,:)=[];

ktError = PercentError(ktRR*pkCERRM(:,1),trueKt(:));
veError = PercentError(veRR*pkCERRM(:,2),trueVe(:));
kepError = PercentError(pkCERRM(:,3),trueKep(:));
vpError = PercentError(ktRR*pkCERRM(:,4),trueVp(:));

ktErrorCE = reshape(ktError,[nX nY repF]);
kepErrorCE = reshape(kepError,[nX nY repF]);
veErrorCE = reshape(veError,[nX nY repF]);
vpErrorCE = reshape(vpError,[nX nY repF]);

pkERRM = real(pkERRM);
ktError = PercentError(ktRR*pkERRM(:,1),trueKt(:));
veError = PercentError(veRR*pkERRM(:,2),trueVe(:));
kepError = PercentError(pkERRM(:,3),trueKep(:));
vpError = PercentError(ktRR*pkERRM(:,4),trueVp(:));

ktErrorE = reshape(ktError,[nX nY repF]);
kepErrorE = reshape(kepError,[nX nY repF]);
veErrorE = reshape(veError,[nX nY repF]);
vpErrorE = reshape(vpError,[nX nY repF]);

ktError = PercentError(ktRR*pkLRRM(:,1),trueKt(:));
veError = PercentError(veRR*pkLRRM(:,2),trueVe(:));
kepError = PercentError(pkLRRM(:,3),trueKep(:));

ktErrorR = reshape(ktError,[nX nY repF]);
kepErrorR = reshape(kepError,[nX nY repF]);
veErrorR = reshape(veError,[nX nY repF]);

% Missing axes labels
figure('Position',[300 300 1200 500])

subplot(1,3,1)
showErrMap(logModulus(flipud(mean(kepErrorR,3))),3)
title('RRM - kep Percent Error')

subplot(1,3,2)
showErrMap(flipud(mean(kepErrorE,3)),0.1)
title('ERRM - kep Percent Error')

subplot(1,3,3)
showErrMap(flipud(mean(kepErrorCE,3)),0.1)
title('CERRM - kep Percent Error')
%% Show the median percent error and interquartile for RRM at chosen vp

% Recompute errors from previous loaded data
ktError = PercentError(ktRR*pkLRRM(:,1),trueKt(:));
veError = PercentError(veRR*pkLRRM(:,2),trueVe(:));
kepError = PercentError(pkLRRM(:,3),trueKep(:));

chosenVp = 0.01;
estsKt = ktError(trueVp(:)==chosenVp,:);
estsKep = kepError(trueVp(:)==chosenVp,:);
estsVe = veError(trueVp(:)==chosenVp,:);

disp(['Median percent error and interquartile range for RRM when vp = ' num2str(chosenVp)])
disp('For Ktrans:')
disp([median(estsKt(:)) iqr(estsKt(:))])
disp('For kep:')
disp([median(estsKep(:)) iqr(estsKep(:))])
disp('For ve:')
disp([median(estsVe(:)) iqr(estsVe(:))])