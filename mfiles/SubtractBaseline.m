function [ outData ] = SubtractBaseline(imgData, numBaselines)

    if nargin<2
        numBaselines=1;
    end

    S = size(imgData);
    
    switch length(S)
        case 4
            baseLine = mean(imgData(:,:,:,1:numBaselines),4);
            outData = imgData - repmat(baseLine,[1 1 1 S(end)]);
        case 3
            baseLine = mean(imgData(:,:,1:numBaselines),3);
            outData = imgData - repmat(baseLine,[1 1 S(end)]);
        case 2
            baseLine = mean(imgData(:,1:numBaselines),2);
            outData = imgData - repmat(baseLine,[1 S(end)]);
        case 1
            baseLine = mean(imgData(1:numBaselines));
            outData = imgData - baseLine;
        otherwise
            outData = NaN;
    end

end

