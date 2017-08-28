% Extra script that converts percent-error maps to json format.
% The json files are used by the interactive figures. 

% Minimal comments

addpath('./mfiles')
clearvars

inDir = './data/mapDownsample';
outDir = './data/mapDownsampleJSon';

if ~exist(outDir,'dir')
    mkdir(outDir);
end

%%
load('./data/simMap.mat');
TRes = 1:60;
listSigmaC = [0:0.01:0.05];

oCp = Cp;
oT = t;
ktRR = 0.07;
kepRR = 0.5;
veRR = ktRR/kepRR;

repF = 100;
trueKt = repmat(trueKt(:),[repF 1]);
trueKep = repmat(trueKep(:),[repF 1]);
trueVe = repmat(trueVe(:),[repF 1]);
trueVp = repmat(trueVp(:),[repF 1]);

%%
sX=nX; sY=nY;
veFilterQuantiles = [0.01 0.99];
mapDim = [sX sY repF];

for j=1:length(listSigmaC)
for i=1:length(TRes)
    curFile = ['Downsample-Noise-' num2str(j) '-TRes-' num2str(i) '.mat'];
    load(fullfile(inDir, 'ETM', curFile));
    load(fullfile(inDir, 'CERRM', curFile));
    
    numImagE(i,j) = sum( sum(abs(imag(pkERRM')))>0 );
    pkERRM = real(pkERRM);
    goodERRM = (trueKep~=kepRR);
    
    % ETM
    ktError = PercentError(pkETM(:,1),trueKt(:));
    veError = PercentError(pkETM(:,1)./pkETM(:,2),trueVe(:));
    kepError = PercentError(pkETM(:,2),trueKep(:));
    vpError = PercentError(pkETM(:,3),trueVp(:));
    
    [ktMean ktStd]  = summarizeErrorMap(ktError,[sX sY repF]);
    [veMean veStd] = summarizeErrorMap(veError,[sX sY repF]);
    [kepMean kepStd] = summarizeErrorMap(kepError,[sX sY repF]);
    [vpMean vpStd] = summarizeErrorMap(vpError,[sX sY repF]);
    
    jsonKt = x2json(ktMean','KTrans');
    jsonVe = x2json(veMean','ve');
    jsonKep = x2json(kepMean','kep');
    jsonVp = x2json(vpMean','vp');
    
    outFileName = ['errMap-ETM-Noise-' num2str(j) '-TRes-' num2str(i) '.json']; 

    fid = fopen(fullfile(outDir,outFileName),'w');
    fprintf(fid,'{\n%s,\n%s,\n%s,\n%s\n}',jsonKt,jsonVe,jsonKep,jsonVp);
    fclose(fid);

    % TM
    ktError = PercentError(pkTM(:,1),trueKt(:));
    veError = PercentError(pkTM(:,1)./pkTM(:,2),trueVe(:));
    kepError = PercentError(pkTM(:,2),trueKep(:));
    
    [ktMean ktStd]  = summarizeErrorMap(ktError,[sX sY repF]);
    [veMean veStd] = summarizeErrorMap(veError,[sX sY repF]);
    [kepMean kepStd] = summarizeErrorMap(kepError,[sX sY repF]);
    
    jsonKt = x2json(ktMean','KTrans');
    jsonVe = x2json(veMean','ve');
    jsonKep = x2json(kepMean','kep');
    
    outFileName = ['errMap-TM-Noise-' num2str(j) '-TRes-' num2str(i) '.json'];  
    fid = fopen(fullfile(outDir,outFileName),'w');
    fprintf(fid,'{\n%s,\n%s,\n%s,\n%s\n}',jsonKt,jsonVe,jsonKep,jsonVp);
    fclose(fid);

    % ERTM
    ktError = PercentError(pkERTM(:,1),trueKt(:));
    veError = PercentError(pkERTM(:,1)./pkERTM(:,2),trueVe(:));
    kepError = PercentError(pkERTM(:,2),trueKep(:));
    vpError = PercentError(pkERTM(:,3),trueVp(:));
    
    [ktMean ktStd]  = summarizeErrorMap(ktError,[sX sY repF]);
    [veMean veStd] = summarizeErrorMap(veError,[sX sY repF]);
    [kepMean kepStd] = summarizeErrorMap(kepError,[sX sY repF]);
    [vpMean vpStd] = summarizeErrorMap(vpError,[sX sY repF]);
    
    jsonKt = x2json(ktMean','KTrans');
    jsonVe = x2json(veMean','ve');
    jsonKep = x2json(kepMean','kep');
    jsonVp = x2json(vpMean','vp');
    
    outFileName = ['errMap-ERTM-Noise-' num2str(j) '-TRes-' num2str(i) '.json'];  
    fid = fopen(fullfile(outDir,outFileName),'w');
    fprintf(fid,'{\n%s,\n%s,\n%s,\n%s\n}',jsonKt,jsonVe,jsonKep,jsonVp);
    fclose(fid);
    
    % RTM
    ktError = PercentError(pkRTM(:,1),trueKt(:));
    veError = PercentError(pkRTM(:,1)./pkRTM(:,2),trueVe(:));
    kepError = PercentError(pkRTM(:,2),trueKep(:));
    
    [ktMean ktStd]  = summarizeErrorMap(ktError,[sX sY repF]);
    [veMean veStd] = summarizeErrorMap(veError,[sX sY repF]);
    [kepMean kepStd] = summarizeErrorMap(kepError,[sX sY repF]);
    
    jsonKt = x2json(ktMean','KTrans');
    jsonVe = x2json(veMean','ve');
    jsonKep = x2json(kepMean','kep');
    
    outFileName = ['errMap-RTM-Noise-' num2str(j) '-TRes-' num2str(i) '.json'];  
    fid = fopen(fullfile(outDir,outFileName),'w');
    fprintf(fid,'{\n%s,\n%s,\n%s,\n%s\n}',jsonKt,jsonVe,jsonKep,jsonVp);
    fclose(fid);
    
    % CERRM
    ktError = PercentError(ktRR*pkCERRM(:,1),trueKt(:));
    veError = PercentError(veRR*pkCERRM(:,2),trueVe(:));
    kepError = PercentError(pkCERRM(:,3),trueKep(:));
    vpError = PercentError(ktRR*pkCERRM(:,4),trueVp(:));
    
    [ktMean ktStd]  = summarizeErrorMap(ktError,[sX sY repF]);
    [veMean veStd] = summarizeErrorMap(veError,[sX sY repF]);
    [kepMean kepStd] = summarizeErrorMap(kepError,[sX sY repF]);
    [vpMean vpStd] = summarizeErrorMap(vpError,[sX sY repF]);
    
    jsonKt = x2json(ktMean','KTrans');
    jsonVe = x2json(veMean','ve');
    jsonKep = x2json(kepMean','kep');
    jsonVp = x2json(vpMean','vp');
    
    outFileName = ['errMap-CERRM-Noise-' num2str(j) '-TRes-' num2str(i) '.json'];  
    fid = fopen(fullfile(outDir,outFileName),'w');
    fprintf(fid,'{\n%s,\n%s,\n%s,\n%s\n}',jsonKt,jsonVe,jsonKep,jsonVp);
    fclose(fid);
    
    % ERRM
    ktError = PercentError(ktRR*pkERRM(:,1),trueKt(:));
    veError = PercentError(veRR*pkERRM(:,2),trueVe(:));
    kepError = PercentError(pkERRM(:,3),trueKep(:));
    vpError = PercentError(ktRR*pkERRM(:,4),trueVp(:));
    
    [ktMean ktStd]  = summarizeErrorMap(ktError,[sX sY repF]);
    [veMean veStd] = summarizeErrorMap(veError,[sX sY repF]);
    [kepMean kepStd] = summarizeErrorMap(kepError,[sX sY repF]);
    [vpMean vpStd] = summarizeErrorMap(vpError,[sX sY repF]);
    
    jsonKt = x2json(ktMean','KTrans');
    jsonVe = x2json(veMean','ve');
    jsonKep = x2json(kepMean','kep');
    jsonVp = x2json(vpMean','vp');
    
    outFileName = ['errMap-ERRM-Noise-' num2str(j) '-TRes-' num2str(i) '.json'];  
    fid = fopen(fullfile(outDir,outFileName),'w');
    fprintf(fid,'{\n%s,\n%s,\n%s,\n%s\n}',jsonKt,jsonVe,jsonKep,jsonVp);
    fclose(fid);
    
    % RRM
    ktError = PercentError(ktRR*pkLRRM(:,1),trueKt(:));
    veError = PercentError(veRR*pkLRRM(:,2),trueVe(:));
    kepError = PercentError(pkLRRM(:,3),trueKep(:));
    
    [ktMean ktStd]  = summarizeErrorMap(ktError,[sX sY repF]);
    [veMean veStd] = summarizeErrorMap(veError,[sX sY repF]);
    [kepMean kepStd] = summarizeErrorMap(kepError,[sX sY repF]);
    
    jsonKt = x2json(ktMean','KTrans');
    jsonVe = x2json(veMean','ve');
    jsonKep = x2json(kepMean','kep');
    
    outFileName = ['errMap-RRM-Noise-' num2str(j) '-TRes-' num2str(i) '.json'];  
    fid = fopen(fullfile(outDir,outFileName),'w');
    fprintf(fid,'{\n%s,\n%s,\n%s,\n%s\n}',jsonKt,jsonVe,jsonKep,jsonVp);
    fclose(fid);
    
    % CLRRM
    ktError = PercentError(ktRR*pkCLRRM(:,1),trueKt(:));
    veError = PercentError(veRR*pkCLRRM(:,2),trueVe(:));
    kepError = PercentError(pkCLRRM(:,3),trueKep(:));
    
    [ktMean ktStd]  = summarizeErrorMap(ktError,[sX sY repF]);
    [veMean veStd] = summarizeErrorMap(veError,[sX sY repF]);
    [kepMean kepStd] = summarizeErrorMap(kepError,[sX sY repF]);
    
    jsonKt = x2json(ktMean','KTrans');
    jsonVe = x2json(veMean','ve');
    jsonKep = x2json(kepMean','kep');
    
    outFileName = ['errMap-CLRRM-Noise-' num2str(j) '-TRes-' num2str(i) '.json'];  
    fid = fopen(fullfile(outDir,outFileName),'w');
    fprintf(fid,'{\n%s,\n%s,\n%s,\n%s\n}',jsonKt,jsonVe,jsonKep,jsonVp);
    fclose(fid);
end
end
return
%%
