function [ outX, outMask ] = quantileFilter( inX, qtRange )
    % Removes values that are outside the quantiles in qtRange

    qtVals = quantile(inX,qtRange);
    outMask = inX>qtVals(1) & inX<qtVals(2);
    outX = inX(outMask);

end

