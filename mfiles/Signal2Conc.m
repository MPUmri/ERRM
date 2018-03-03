function [ Ct ] = Signal2Conc(dataSignal, dataT1, dataM0, TR, FlipAngle, r1)
    % dataT1 in ms
    % TR in ms
    % FlipAngle in degrees
    % r1 in 1/(mM s)

    %%  
    imgSize = size(dataSignal);
    sT = imgSize(end);
    
    alpha = FlipAngle * pi /180;
    
    
    eR10 = repmat(exp(-TR ./ dataT1),[1 1 1 sT]);
    B = (1-eR10) ./ (1-cos(alpha) .* eR10);
    A =  SubtractBaseline(dataSignal) ./ (repmat(dataM0,[1 1 1 sT]) .* sin(alpha));

    num = 1 - (A+B);
    den = 1 - cos(alpha) .* (A+B);
    R1 = real(-1/(TR/1000) .* log(num ./ den));
    Ct = (R1 - repmat(1000./dataT1,[1 1 1 sT])) ./ r1;
    
    Ct(isnan(Ct))=0;
end

