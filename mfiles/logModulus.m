function [ outVals ] = logModulus( inVals )

    outVals = sign(inVals) .* log10(abs(inVals)+1);

end

