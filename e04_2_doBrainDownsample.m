% Downsample GBM data to compare stability of CERRM and ETM

% Estimated runtime: 65 seconds

%% Pre-setup

clearvars
addpath('./mfiles')
rng(12345)

%% Load data
load('./data/gbm/gbmData.mat');
% The loaded data contains:
% t   [70x1] - time (in minutes) at each dce frame
% Cp  [70x1] - concentration in blood plasma
% Crr [70x1] - concentration in reference region
% Ct  [70x10813] - concentration in tumour voxels (10813 voxels)

tic

%% Fit to high temporal resolution (5.4 s - original) data

% CERRM
[tempPkCE, ~, estKepRR(1)]=CERRM(Ct,Crr,t);
% ETM
tempPkET=LLSQ(Ct,Cp,t,1);
% RRM
tempPkRR = LRRM(Ct,Crr,t);

% 'true' parameters are fits at original temporal resolution
trueKtCE = tempPkCE(:,1);
trueVeCE = tempPkCE(:,2);
trueKepCE = tempPkCE(:,3);
trueVpCE = tempPkCE(:,4);

trueKtET = tempPkET(:,1);
trueVeET = tempPkET(:,1)./tempPkET(:,2);
trueKepET = tempPkET(:,2);
trueVpET = tempPkET(:,3);

trueKtRR = tempPkRR(:,1);
trueVeRR = tempPkRR(:,2);
trueKepRR = tempPkRR(:,3);

TRes(1) = 60*(t(2)-t(1));

%% Downsample the data and refit
nVox = size(Ct,2);
for i=2:12  
    % Randomize the first frame when downsampling
    phaseValues = randi([0 i-1],nVox,1);
    clearvars tempPkCE tempPkET pkERRM
    for j=1:nVox
        curCt = downsample(Ct(:,j),i,phaseValues(j));
        curCrr = downsample(Crr,i,phaseValues(j));
        curCp = downsample(Cp,i,phaseValues(j));
        curT = downsample(t,i,phaseValues(j));
        
        pkERRM(j,:)=ERRM(curCt,curCrr,curT,0,0,0); % ERRM
        tempPkET(j,:)=LLSQ(curCt,curCp,curT,1); % ETM
        tempPkRR(j,:)=LRRM(curCt,curCrr,curT); % RRM
    end
    
    % Get estimate for kepRR - from ERRM
    rawKepRR = pkERRM(:,5);
    % Find voxels where all estimates are real and positive
    goodVals = pkERRM(:,1)>0 & pkERRM(:,2)>0 & pkERRM(:,3)>0 & pkERRM(:,4)>0 & pkERRM(:,5)>0 & imag(pkERRM(:,5))==0;
    estKepRR(i) = iqrMean(rawKepRR(goodVals));
    
    % Now fit CERRM
    for j=1:nVox
        curCt = downsample(Ct(:,j),i,phaseValues(j));
        curCrr = downsample(Crr,i,phaseValues(j));
        curCp = downsample(Cp,i,phaseValues(j));
        curT = downsample(t,i,phaseValues(j));
        
        tempPkCE(j,:)=CERRM(curCt,curCrr,curT,estKepRR(i));
    end
    
    % Calculate percent change
    tempErrKtCE=PercentError(tempPkCE(:,1),trueKtCE);
    tempErrVeCE=PercentError(tempPkCE(:,2),trueVeCE);
    tempErrKepCE=PercentError(tempPkCE(:,3),trueKepCE);
    tempErrVpCE=PercentError(tempPkCE(:,4),trueVpCE);
        
    tempErrKtET=PercentError(tempPkET(:,1),trueKtET);
    tempErrVeET=PercentError(tempPkET(:,1)./tempPkET(:,2),trueVeET);
    tempErrKepET=PercentError(tempPkET(:,2),trueKepET);
    tempErrVpET=PercentError(tempPkET(:,3),trueVpET);
    
    tempErrKtRR=PercentError(tempPkRR(:,1),trueKtRR);
    tempErrVeRR=PercentError(tempPkRR(:,2),trueVeRR);
    tempErrKepRR=PercentError(tempPkRR(:,3),trueKepRR);
    
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
    
    qKtRR(i,:)=quantile(tempErrKtCE(:),[.25 .5 .75]);
    qVeRR(i,:)=quantile(tempErrVeCE(:),[.25 .5 .75]);
    qKepRR(i,:)=quantile(tempErrKepCE(:),[.25 .5 .75]);
    
    TRes(i) = 60*(curT(2)-curT(1));
end
toc

%% Export to CSV
outFile = fullfile('./dataResults',['e04-downsampleGbmResults.csv']);
hdr=['FitMethod,TemporalRes,errKt,errVe,errKep,errVp,qt25Kt,qt75Kt,qt25Ve,qt75Ve,qt25Kep,qt75Kep,qt25Vp,qt75Vp'];
outID = fopen(outFile, 'w+');
fprintf(outID, '%s\n', hdr); % Print header into csv file
for i=1:12
   outLine = {'ETM',TRes(i),qKtET(i,2),qVeET(i,2),qKepET(i,2),qVpET(i,2),...
       qKtET(i,1),qKtET(i,3),qVeET(i,1),qVeET(i,3),qKepET(i,1),qKepET(i,3),...
       qVpET(i,1),qVpET(i,3)};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'RRM',TRes(i),qKtRR(i,2),qVeRR(i,2),qKepRR(i,2),NaN,...
       qKtRR(i,1),qKtRR(i,3),qVeRR(i,1),qVeRR(i,3),qKepRR(i,1),qKepRR(i,3),...
       NaN,NaN};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
%    outLine = {'ERTM',TRes(i),qKtERT(i,2),qVeERT(i,2),qKepERT(i,2),qVpERT(i,2),...
%        qKtERT(i,1),qKtERT(i,3),qVeERT(i,1),qVeERT(i,3),qKepERT(i,1),qKepERT(i,3),...
%        qVpERT(i,1),qVpERT(i,3)};
%    fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
   
   outLine = {'CERRM',TRes(i),qKtCE(i,2),qVeCE(i,2),qKepCE(i,2),qVpCE(i,2),...
       qKtCE(i,1),qKtCE(i,3),qVeCE(i,1),qVeCE(i,3),qKepCE(i,1),qKepCE(i,3),...
       qVpCE(i,1),qVpCE(i,3)};
   fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
end
fclose(outID)
%%
toc

%% Median +/- interquartile range vs Temporal Resolution
lErr = 1;
hErr = 3;

cCE = [.5 .2 .55];
cET = [.3 .7 .6];

figure('Position',[300 300 320 500])
hold on
shadedErrorPlot(TRes',qKtET(:,2),qKtET(:,lErr),qKtET(:,hErr),cET,1);
shadedErrorPlot(TRes',qKtCE(:,2),qKtCE(:,lErr),qKtCE(:,hErr),cCE,1);
hold off
ylim([-100 100])
xlabel('Temporal Resolution [s]')
ylabel('Percent Error')
title('Ktrans')

figure('Position',[300 300 320 500])
hold on
shadedErrorPlot(TRes',qKepET(:,2),qKepET(:,lErr),qKepET(:,hErr),cET,1);
shadedErrorPlot(TRes',qKepCE(:,2),qKepCE(:,lErr),qKepCE(:,hErr),cCE,1);
hold off
ylim([-100 100])
xlabel('Temporal Resolution [s]')
ylabel('Percent Error')
title('kep')

figure('Position',[300 300 320 500])
hold on
shadedErrorPlot(TRes',qVeET(:,2),qVeET(:,lErr),qVeET(:,hErr),cET,1);
shadedErrorPlot(TRes',qVeCE(:,2),qVeCE(:,lErr),qVeCE(:,hErr),cCE,1);
hold off
ylim([-100 100])
xlabel('Temporal Resolution [s]')
ylabel('Percent Error')
title('ve')

figure('Position',[300 300 320 500])
hold on
shadedErrorPlot(TRes',qVpET(:,2),qVpET(:,lErr),qVpET(:,hErr),cET,1);
shadedErrorPlot(TRes',qVpCE(:,2),qVpCE(:,lErr),qVpCE(:,hErr),cCE,1);
hold off
ylim([-100 300])
xlabel('Temporal Resolution [s]')
ylabel('Percent Error')
title('vp')