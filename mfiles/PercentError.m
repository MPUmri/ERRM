function [percentErr, meanPercentErr, stdPercentErr, numBadEstimates] = PercentError(estVals, trueVals)

    if isrow(trueVals)
        trueVals = trueVals';
    end
    
    if size(estVals) ~= size(trueVals)
        [nX nY] = size(estVals);
        if nX == length(trueVals)
            trueVals = repmat(trueVals,[1 nY]);
        else
            trueVals = repmat(trueVals',[nX 1]);
        end
    end
    percentErr = 100*(estVals-trueVals)./trueVals;
    finiteVals = isfinite(percentErr);
    meanPercentErr = nanmean(percentErr(finiteVals(:)));
    stdPercentErr = nanstd(percentErr(finiteVals(:)));
    numBadEstimates = sum(~finiteVals);
end

