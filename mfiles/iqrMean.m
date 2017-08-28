function [ iqrM, iqrStd, qtRange ] = iqrMean( X )
    % x can either be a vector of length M (column, preferrably)
    %       or a MxN matrix
    % output is a single value (if input is vector) or vector of length N
    
    if isrow(X)
        X = X';
    end
    [nS nT] = size(X);
    
    iqrM = zeros(nT,1);
    iqrStd = iqrM;
    for i=1:nT
        x = X(:,i);
        qtRange = quantile(x,[.25 .75]);
        xMask = x>qtRange(1) & x<qtRange(2);
        iqrM(i) = mean(x(xMask));
        iqrStd(i) = std(x(xMask));
    end
end

