clearvars
addpath('./mfiles')

indTRes = [5, 15, 30]; % Hist order: top-to-down
indNoise = [2, 4, 6]; % Hist order: left-to-right
doPlot = true;

vMax = 1;
vMin = 0;

if doPlot
    figure('Position',[300 300 900 500])
end

for i=1:length(indNoise)
    for j=1:length(indTRes)
        curFile = ['Downsample-Noise-' num2str(indNoise(i)) '-TRes-' num2str(indTRes(j)) '.mat'];
        load(fullfile('./data/mapDownsample/CERRM',curFile));
        
        estKepRR = pkERRM(:,5); % All of the kepRR estimate from ERRM
        
        % Get number of imaginary estimates
        isImag = imag(estKepRR) ~= 0;
        numImag(i,j) = sum(isImag);
        
        % Get number of negative estimates
        estKepRR = estKepRR(~isImag);
        numNegative(i,j) = sum(estKepRR<0);
        
        % Get the interquartile mean of positive real kepRR estimates
        kepRR(i,j) = iqrMean(estKepRR(estKepRR>0));
        
        if doPlot
            vMask = estKepRR>vMin & estKepRR<vMax;
            subplot(length(indTRes),length(indNoise),i+length(indNoise)*(j-1))
            histogram(estKepRR,40,'BinLimits',[vMin vMax])
        end
        
    end
end

% Display the percent of fits which were negative or complex
% These were annotated onto the figure afterwards
disp('Percent of fits which were negative')
disp(100*numNegative./length(isImag))
disp('Percent of fits which were complex')
disp(100*numImag./length(isImag))