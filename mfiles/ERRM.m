function [pkParams, fittedCt] = ERRM(Ct, Crr, t, doNonNeg, doPure, doReal)
  % Extended Linear Reference Region fit - Doing double integrals
  % %%%%%%%%
  % Inputs:
  % Ct [MxN] - Concentration in tissue of interest for N voxels
  % Crr [Mx1] - Concentration in reference region
  % t [Mx1] - Time (in minutes) at each timepoint
  % %% 
  % Optional inputs: (default value in square brackets)
  % doNonNeg [Bool-false] - use non-negative linear least squares
  % doPure [Bool-false] - keep fitting parameters without transforming them
  % doReal [Bool-true] - convert complex estimates to real
  % %%%%%%%%
  % Outputs:
  % pkParams [Mx5] - [kt/ktRR, ve/veRR, kep, vp/ktRR, kepRR]
  % fittedCt [MxN] - fitted curve, fits to cumtrapz(t,Ct)
  
  %% Potential changes to code
  % - Get rid of 'doNonNeg' since NNLS does not provide improvement
  % - Get rid of 'doPure' and return original fitted params as additional output
  % - Instead of 'doReal' being boolean, make it a switch where one option
  % allows complex fits to be discarded
  
  %% Set optional parameters to default values if they are not supplied

    % By default, do not use NonNegative Linear Least Squares
    if nargin<4
        doNonNeg = false;
    end
    % By default, we re-arrange the fitted parameters into more useful forms
    if nargin<5
        doPure = false;
    end
    % By default, only output the real fitting parameters
    if nargin<6
        doReal = true;
    end
    
    %% Do the fitting
    stepSize = t(2)-t(1); % Assuming constant step size throughout acquisition
    
    % Initialize matrices
    [sT, sX] = size(Ct);
    
    M1 = stepSize*cumtrapz(Crr);
    M2 = stepSize*stepSize*cumtrapz(cumtrapz(Crr));
    M4 = Crr;

    fittedCt = zeros(sT,sX);
    pkParams = zeros(sX,4);

    % Solve the linear problem
    % Disabling warnings for speed (possible non-tissue regions in image give warnings)
    warning off
    
    if doNonNeg
        % Use non-negative linear least squares
        for i=1:sX
                curCt = stepSize*cumtrapz(squeeze(Ct(:,i)));
                M3 = -stepSize*cumtrapz(curCt);
                M = [M1, M2, M3, M4];
                pkParams(i,:) = lsqnonneg(M,curCt);
                fittedCt(:,i) = M*pkParams(i,:)';
        end
    else
        % Use simple linear least squares
        for i=1:sX
                curCt = stepSize*cumtrapz(squeeze(Ct(:,i)));
                M3 = -stepSize*cumtrapz(curCt);
                M = [M1, M2, M3, M4];
                pkParams(i,:) = lscov(M,curCt);
                fittedCt(:,i) = M*pkParams(i,:)';
        end
    end
    warning on
    
    % Transform the fitted parameters into a useful form
    % pkParams = [MESSY, MESSY, kep, vp/ktRR]
    % desired = [kt/ktRR, ve/veRR, kep, vp/ktRR, kepRR]
    if doPure==false
        A = pkParams(:,1)./pkParams(:,4);
        B = pkParams(:,2)./pkParams(:,4);
        kepRR = ( A - sqrt(A.^2 - 4*B) ) / 2;
        ktvp = A - pkParams(:,3) - kepRR;
        ktRel = pkParams(:,4) .* ktvp;
        veRel = ktRel .* kepRR ./ pkParams(:,3);
        pkParams = [ktRel, veRel, pkParams(:,3), pkParams(:,4), kepRR];
    end
    
    if doReal
        pkParams = real(pkParams);
    end
    
end
