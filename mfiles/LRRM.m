function [pkParams, fittedCt, M] = LRRM(Ct, Crr, t, doPure)
  % Linear Reference Region fit
%%
    if nargin<4
        doPure = false;
    end
%%
    stepSize = t(2)-t(1); % Assuming constant step size throughout acquisition
    
    % Initialize matrices
    [sT, sX] = size(Ct);
    
    M1 = zeros(sT,1);
    M2 = M1;

    M1 = Crr;
    M2 =  stepSize*cumtrapz(Crr);

    fittedCt = zeros(sT,sX);
    pkParams = zeros(sX,3);

    % Solve the linear problem
    % Disabling warnings for speed (possible non-tissue regions in image give warnings)
    warning off

    for i=1:sX
            curCt = squeeze(Ct(:,i));
            M3 = -stepSize*cumtrapz(curCt);
            M = [M1, M2, M3];
            pkParams(i,:) = lscov(M,curCt);
            fittedCt(:,i) = M*pkParams(i,:)';
    end
    warning on
    
    % pkParams = [kTrans_TOI/kTrans_RR, kTrans_TOI/ve_RR, kTrans_TOI/ve_TOI]
    % desire: pkParams = [kTrans_TOI/ktrans_RR, ve_TOI/ve_RR, kTrans_TOI/ve_TOI]
    if doPure==false
        pkParams(:,2) = pkParams(:,2)./pkParams(:,3);
    end
end
