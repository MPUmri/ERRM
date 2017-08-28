function [ meanMap, stdMap ] = summarizeErrorMap( rawError, arrayDim, repDim )

    if nargin<3
        repDim = 3;
    end

    rawError = reshape(rawError,arrayDim);

    meanMap = squeeze(mean(rawError,repDim));
    stdMap = squeeze(std(rawError,0,repDim));
    
    meanMap = round(meanMap,3);
    stdMap = round(stdMap,3);
end

