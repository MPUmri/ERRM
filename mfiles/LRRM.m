function [pkParams, fittedCt, originalFittedParams] = LRRM(Ct, Crr, t)
  % Linear Reference Region Model fitting, based on: 
  % - Li, et al. (2009) Medical Physics, 36, 3786–3794. http://doi.org/10.1118/1.3152113
  % - Cárdenas-Rodríguez, et al. (2013) MRI 31(4), 497–507. http://doi.org/10.1016/j.mri.2012.10.008
  % %%%%%%%%
  % Inputs:
  % Ct [MxN] - Concentration in tissue of interest for N voxels
  % Crr [Mx1] - Concentration in reference region
  % t [Mx1] - Time (in minutes) at each timepoint
  % %% 
  % Optional inputs: (default value in square brackets)
  % doPure [Bool-false] - keep fitting parameters without transforming them
  % %%%%%%%%
  % Outputs:
  % pkParams [Mx3] - [kt/ktRR, ve/veRR, kep]
  % fittedCt [MxN] - fitted curve, fits to cumtrapz(t,Ct)
  % originalFittedParams [Mx3] - [kt/ktRR, kt/veRR, kep]
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
    
    originalFittedParams = pkParams;
    % pkParams = [kTrans_TOI/kTrans_RR, kTrans_TOI/ve_RR, kTrans_TOI/ve_TOI]
    % desire: pkParams = [kTrans_TOI/ktrans_RR, ve_TOI/ve_RR, kTrans_TOI/ve_TOI]
    pkParams(:,2) = pkParams(:,2)./pkParams(:,3);
end
