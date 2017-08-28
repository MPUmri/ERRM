function [ outX, outMask ] = quantileFilter( inX, qtRange )

    qtVals = quantile(inX,qtRange);
    outMask = inX>qtVals(1) & inX<qtVals(2);
    outX = inX(outMask);

end

