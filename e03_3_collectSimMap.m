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
TRes = 1:60;
listSigmaC = [0:0.01:0.05];

ktRR = 0.07;
kepRR = 0.5;
veRR = ktRR/kepRR;

repF = 100;
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
    % Use only real part when calculating percent error
    pkERRM = real(pkERRM);
    % For ERRM: exclude fits where tissue of interest and reference tissue
    % have same kep values, as the ERRM fails to fit those cases
    goodERRM = (trueKep~=kepRR);
    
    %% Apply vpCutOff
    % Does nothing if vpCutoff = [-inf inf]
    pkCERRM(~vpMask,:)=[];
    pkCLRRM(~vpMask,:)=[];
    pkERRM(~vpMask,:)=[];
    pkERTM(~vpMask,:)=[];
    pkETM(~vpMask,:)=[];
    pkLRRM(~vpMask,:)=[];
    pkRTM(~vpMask,:)=[];
    pkTM(~vpMask,:)=[];
    %% Collect mean, std, and quartiles of the percent errors
    
    % TM
    [ktError, meanKtErrT(i,j), stdKtErrT(i,j)] = PercentError(pkTM(:,1),trueKt(:));
    [veError, meanVeErrT(i,j), stdVeErrT(i,j)] = PercentError(pkTM(:,1)./pkTM(:,2),trueVe(:));
    [kepError, meanKepErrT(i,j), stdKepErrT(i,j)] = PercentError(pkTM(:,2),trueKep(:));
    
    qKtT(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
    qVeT(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
    qKepT(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);
    
    % ETM
    [ktError, meanKtErrET(i,j), stdKtErrET(i,j)] = PercentError(pkETM(:,1),trueKt(:));
    [veError, meanVeErrET(i,j), stdVeErrET(i,j)] = PercentError(pkETM(:,1)./pkETM(:,2),trueVe(:));
    [kepError, meanKepErrET(i,j), stdKepErrET(i,j)] = PercentError(pkETM(:,2),trueKep(:));
    [vpError, meanVpErrET(i,j), stdVpErrET(i,j)] = PercentError(pkETM(:,3),trueVp(:));
    
    qKtET(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
    qVeET(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
    qKepET(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);
    qVpET(i,j,:) = quantile(vpError,[.05 .25 .5 .75 .95]);
    
    % RTM-TM
    [ktError, meanKtErrRT(i,j), stdKtErrRT(i,j)] = PercentError(pkRTM(:,1),trueKt(:));
    [veError, meanVeErrRT(i,j), stdVeErrRT(i,j)] = PercentError(pkRTM(:,1)./pkRTM(:,2),trueVe(:));
    [kepError, meanKepErrRT(i,j), stdKepErrRT(i,j)] = PercentError(pkRTM(:,2),trueKep(:));
    
    qKtRT(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
    qVeRT(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
    qKepRT(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);
    
    % ERTM (or RTM-ETM)
    [ktError, meanKtErrERT(i,j), stdKtErrERT(i,j)] = PercentError(pkERTM(:,1),trueKt(:));
    [veError, meanVeErrERT(i,j), stdVeErrERT(i,j)] = PercentError(pkERTM(:,1)./pkERTM(:,2),trueVe(:));
    [kepError, meanKepErrERT(i,j), stdKepErrERT(i,j)] = PercentError(pkERTM(:,2),trueKep(:));
    [vpError, meanVpErrERT(i,j), stdVpErrERT(i,j)] = PercentError(pkERTM(:,3),trueVp(:));
    
    qKtERT(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
    qVeERT(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
    qKepERT(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);
    qVpERT(i,j,:) = quantile(vpError,[.05 .25 .5 .75 .95]);
    
    % CERRM
    [ktError, meanKtErrCE(i,j), stdKtErrCE(i,j), badValKtCE(i,j)] = PercentError(ktRR*pkCERRM(:,1),trueKt(:));
    [veError, meanVeErrCE(i,j), stdVeErrCE(i,j), badValVeCE(i,j)] = PercentError(veRR*pkCERRM(:,2),trueVe(:));
    [kepError, meanKepErrCE(i,j), stdKepErrCE(i,j), badValKepCE(i,j)] = PercentError(pkCERRM(:,3),trueKep(:));
    [vpError, meanVpErrCE(i,j), stdVpErrCE(i,j), badValVpCE(i,j)] = PercentError(ktRR*pkCERRM(:,4),trueVp(:));
    
    qKtCE(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
    qVeCE(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
    qKepCE(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);
    qVpCE(i,j,:) = quantile(vpError,[.05 .25 .5 .75 .95]);
    
    % ERRM
    [ktError, meanKtErrE(i,j), stdKtErrE(i,j), badValKtE(i,j)] = PercentError(ktRR*pkERRM(:,1),trueKt(:));
    [veError, meanVeErrE(i,j), stdVeErrE(i,j), badValVeE(i,j)] = PercentError(veRR*pkERRM(:,2),trueVe(:));
    [kepError, meanKepErrE(i,j), stdKepErrE(i,j), badValKepE(i,j)] = PercentError(pkERRM(:,3),trueKep(:));
    [vpError, meanVpErrE(i,j), stdVpErrE(i,j), badValVpE(i,j)] = PercentError(ktRR*pkERRM(:,4),trueVp(:));
    
    qKtE(i,j,:) = quantile(ktError(goodERRM),[.05 .25 .5 .75 .95]);
    qVeE(i,j,:) = quantile(veError(goodERRM),[.05 .25 .5 .75 .95]);
    qKepE(i,j,:) = quantile(kepError(goodERRM),[.05 .25 .5 .75 .95]);
    qVpE(i,j,:) = quantile(vpError,[.05 .25 .5 .75 .95]);
    
    % RRM
    [ktError, meanKtErrR(i,j), stdKtErrR(i,j)] = PercentError(ktRR*pkLRRM(:,1),trueKt(:));
    [veError, meanVeErrR(i,j), stdVeErrR(i,j)] = PercentError(veRR*pkLRRM(:,2),trueVe(:));
    [kepError, meanKepErrR(i,j), stdKepErrR(i,j)] = PercentError(pkLRRM(:,3),trueKep(:));
    
    qKtR(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
    qVeR(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
    qKepR(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);
    
    % CLRRM - not mentioned in manuscript
    [ktError, meanKtErrCR(i,j), stdKtErrCR(i,j)] = PercentError(ktRR*pkCLRRM(:,1),trueKt(:));
    [veError, meanVeErrCR(i,j), stdVeErrCR(i,j)] = PercentError(veRR*pkCLRRM(:,2),trueVe(:));
    [kepError, meanKepErrCR(i,j), stdKepErrCR(i,j)] = PercentError(pkCLRRM(:,3),trueKep(:));
    
    qKtCR(i,j,:) = quantile(ktError,[.05 .25 .5 .75 .95]);
    qVeCR(i,j,:) = quantile(veError,[.05 .25 .5 .75 .95]);
    qKepCR(i,j,:) = quantile(kepError,[.05 .25 .5 .75 .95]);
    
    %%
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
outFile = fullfile('./dataResults',['e03a-downsampleMapResultsMean.csv']);
hdr=['FitMethod,TemporalRes,sigmaC,errKt,errVe,errKep,errVp,stdKt,stdVe,stdKep,stdVp'];
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
fclose(outID)
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
fclose(outID)
toc

%% Plots

%% Median +/- interquartile range vs StdDev of Noise

myX = listSigmaC;
myTres = 1;
lErr = 2;
hErr = 4;

figure('Position',[300 300 300 500])
errorbar(myX,qKtERT(myTres,:,3),abs(qKtERT(myTres,:,lErr)-qKtERT(myTres,:,3)),abs(qKtERT(myTres,:,hErr)-qKtERT(myTres,:,3)))
hold on
errorbar(myX,qKtCE(myTres,:,3),abs(qKtCE(myTres,:,lErr)-qKtCE(myTres,:,3)),abs(qKtCE(myTres,:,hErr)-qKtCE(myTres,:,3)))
hold on
errorbar(myX,qKtE(myTres,:,3),abs(qKtE(myTres,:,lErr)-qKtE(myTres,:,3)),abs(qKtE(myTres,:,hErr)-qKtE(myTres,:,3)))
hold off
ylim([-20 20])
xlabel('StdDev [mM]')
ylabel('Percent Error')
title('Ktrans')
legend('ERTM','CERRM','ERRM')

figure('Position',[300 300 300 500])
errorbar(myX,qKepERT(myTres,:,3),abs(qKepERT(myTres,:,lErr)-qKepERT(myTres,:,3)),abs(qKepERT(myTres,:,hErr)-qKepERT(myTres,:,3)))
hold on
errorbar(myX,qKepCE(myTres,:,3),abs(qKepCE(myTres,:,lErr)-qKepCE(myTres,:,3)),abs(qKepCE(myTres,:,hErr)-qKepCE(myTres,:,3)))
hold on
errorbar(myX,qKepE(myTres,:,3),abs(qKepE(myTres,:,lErr)-qKepE(myTres,:,3)),abs(qKepE(myTres,:,hErr)-qKepE(myTres,:,3)))
hold off
ylim([-20 20])
xlabel('StdDev [mM]')
ylabel('Percent Error')
title('kep')
legend('ERTM','CERRM','ERRM')

figure('Position',[300 300 300 500])
errorbar(myX,qVeERT(myTres,:,3),abs(qVeERT(myTres,:,lErr)-qVeERT(myTres,:,3)),abs(qVeERT(myTres,:,hErr)-qVeERT(myTres,:,3)))
hold on
errorbar(myX,qVeCE(myTres,:,3),abs(qVeCE(myTres,:,lErr)-qVeCE(myTres,:,3)),abs(qVeCE(myTres,:,hErr)-qVeCE(myTres,:,3)))
hold on
errorbar(myX,qVeE(myTres,:,3),abs(qVeE(myTres,:,lErr)-qVeE(myTres,:,3)),abs(qVeE(myTres,:,hErr)-qVeE(myTres,:,3)))
hold off
ylim([-20 20])
xlabel('StdDev [mM]')
ylabel('Percent Error')
title('ve')
legend('ERTM','CERRM','ERRM')

figure('Position',[300 300 300 500])
errorbar(myX,qVpERT(myTres,:,3),abs(qVpERT(myTres,:,lErr)-qVpERT(myTres,:,3)),abs(qVpERT(myTres,:,hErr)-qVpERT(myTres,:,3)))
hold on
errorbar(myX,qVpCE(myTres,:,3),abs(qVpCE(myTres,:,lErr)-qVpCE(myTres,:,3)),abs(qVpCE(myTres,:,hErr)-qVpCE(myTres,:,3)))
hold on
errorbar(myX,qVpE(myTres,:,3),abs(qVpE(myTres,:,lErr)-qVpE(myTres,:,3)),abs(qVpE(myTres,:,hErr)-qVpE(myTres,:,3)))
hold off
ylim([-50 50])
xlabel('StdDev [mM]')
ylabel('Percent Error')
title('vp')
legend('ERTM','CERRM','ERRM')

%% Ratios of interquartile range between ERRM and CERRM
quartRatio_Kt=(qKtE(myTres,:,hErr)-qKtE(myTres,:,lErr))./(qKtCE(myTres,:,hErr)-qKtCE(myTres,:,lErr))
quartRatio_Kep=(qKepE(myTres,:,hErr)-qKepE(myTres,:,lErr))./(qKepCE(myTres,:,hErr)-qKepCE(myTres,:,lErr))
quartRatio_Ve=(qVeE(myTres,:,hErr)-qVeE(myTres,:,lErr))./(qVeCE(myTres,:,hErr)-qVeCE(myTres,:,lErr))
quartRatio_Vp=(qVpE(myTres,:,hErr)-qVpE(myTres,:,lErr))./(qVpCE(myTres,:,hErr)-qVpCE(myTres,:,lErr))

%% Median +/- interquartile range vs Temporal Resolution
mySig = 2;
lErr = 2;
hErr = 4;

cCE = [.5 .2 .55];
cET = [.3 .7 .6];

figure('Position',[300 300 300 500])
hold on
shadedErrorPlot(TRes',qKtET(:,mySig,3),qKtET(:,mySig,lErr),qKtET(:,mySig,hErr),cET);
shadedErrorPlot(TRes',qKtCE(:,mySig,3),qKtCE(:,mySig,lErr),qKtCE(:,mySig,hErr),cCE);
hold off
ylim([-100 100])
xlabel('Temporal Resolution [s]')
ylabel('Percent Error')
title('Ktrans')

figure('Position',[300 300 300 500])
hold on
shadedErrorPlot(TRes',qKepET(:,mySig,3),qKepET(:,mySig,lErr),qKepET(:,mySig,hErr),cET);
shadedErrorPlot(TRes',qKepCE(:,mySig,3),qKepCE(:,mySig,lErr),qKepCE(:,mySig,hErr),cCE);
hold off
ylim([-100 100])
xlabel('Temporal Resolution [s]')
ylabel('Percent Error')
title('kep')

figure('Position',[300 300 300 500])
hold on
shadedErrorPlot(TRes',qVeET(:,mySig,3),qVeET(:,mySig,lErr),qVeET(:,mySig,hErr),cET);
shadedErrorPlot(TRes',qVeCE(:,mySig,3),qVeCE(:,mySig,lErr),qVeCE(:,mySig,hErr),cCE);
hold off
ylim([-100 100])
xlabel('Temporal Resolution [s]')
ylabel('Percent Error')
title('ve')

figure('Position',[300 300 300 500])
hold on
shadedErrorPlot(TRes',qVpET(:,mySig,3),qVpET(:,mySig,lErr),qVpET(:,mySig,hErr),cET);
shadedErrorPlot(TRes',qVpCE(:,mySig,3),qVpCE(:,mySig,lErr),qVpCE(:,mySig,hErr),cCE);
hold off
ylim([-100 1000])
xlabel('Temporal Resolution [s]')
ylabel('Percent Error')
title('vp')

%% Make the Error Maps
indTRes = 1;
indSigma = 1;
repF = 100;

curFile = ['Downsample-Noise-' num2str(indSigma) '-TRes-' num2str(indTRes) '.mat'];
load(fullfile(outDir, 'ETM', curFile));
load(fullfile(outDir, 'CERRM', curFile));

% Apply vpCutOff - unnecessary
pkCERRM(~vpMask,:)=[];
pkERRM(~vpMask,:)=[];
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

%% Plot the error maps
% Missing axes labels
figure('Position',[300 300 300 500])
showErrMap(flipud(mean(kepErrorCE,3)),0.1)
title('CERRM-kep')

figure('Position',[300 300 300 500])
showErrMap(flipud(mean(kepErrorE,3)),0.1)
title('ERRM-kep')
%%
