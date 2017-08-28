function [ estPk, estCp ] = doRTM(Ct,Crr,t,pkParamsRR, doExt, doNonNeg)
% Reference tissue method

    if nargin<6
        doNonNeg = false;
    end

    if nargin < 5
        doExt = 1;
    end

    ktRR = pkParamsRR(1);
    veRR = pkParamsRR(2);
    
    estCp = (1/ktRR) * gradient(Crr,t(2)-t(1)) + (1/veRR) * Crr;
    
    estPk = LLSQ(Ct,estCp,t,doExt,doNonNeg);
end

