function [pkParams, fittedCt, kepRR, pkERRM, fittedCtERRM] = CERRM(Ct, Crr, t, kepRR, doNonNeg, doPure)
  % Constrained Extended Linear Reference Region fit - Doing double integrals
  % %%%%%%%%
  % Inputs:
  % Ct [MxN] - Concentration in tissue of interest for N voxels
  % Crr [Mx1] - Concentration in reference region
  % t [Mx1] - Time (in minutes) at each timepoint
  % %% 
  % Optional inputs: (default value in square brackets)
  % kepRR [1x1] - Default: [] - Value of kepRR, if known a-priori.
  % Otherwise, leave as [] to estimate it from ERRM.
  % doNonNeg [Bool-false] - use non-negative linear least squares
  % doPure [Bool-false] - keep fitting parameters without transforming them
  % %%%%%%%%
  % Outputs:
  % pkParams [Mx5] - [kt/ktRR, ve/veRR, kep, vp/ktRR, kepRR]
  % fittedCt [MxN] - fitted curve, fits to cumtrapz(t,Ct)
  % kepRR [1x1] - kepRR value used with CERRM (user provided or estimated
  % from ERRM)
  % Optional outputs: if input kepRR = [], then:
  % pkERRM [Mx5] - [kt/ktRR, ve/veRR, kep, vp/ktRR, kepRR] from ERRM
  % fittedCrERRM [MxN] - fitted curve, fits to cumtrapz(t,Ct)
  
  
  %% Potential changes to code
  % - Get rid of 'doNonNeg' since NNLS does not provide improvement
  % - Get rid of 'doPure' and return original fitted params as additional output
  % - Instead of 'doReal' being boolean, make it a switch where one option
  % allows complex fits to be discarded

  %% Do the fitting

    if nargin < 6
        doPure = false;
    end
    
    if nargin < 5
        doNonNeg = false;
    end

    if nargin < 4
        kepRR = [];
    end
    
    if isempty(kepRR)
       [pkERRM, fittedCtERRM] = ERRM(Ct,Crr,t,doNonNeg);
       rawKepRR = pkERRM(:,5);
       if std(rawKepRR)<1e-3
           % If the estimated kepRR from ERRM is closely grouped, then use median
           % This situation is very unlikely to happen in clinical data
           % Then practical purpose of this is that when simulating
           % noiseless data, the interquartile mean can't be used because
           % there is no fluctuation in the estimated kepRR from ERRM
           kepRR = nanmedian(rawKepRR);
       else
           % Find voxels where all estimates are real and positive
           goodVals = pkERRM(:,1)>0 & pkERRM(:,2)>0 & pkERRM(:,3)>0 & pkERRM(:,4)>0 & pkERRM(:,5)>0 & imag(pkERRM(:,5))==0;
           kepRR = iqrMean(rawKepRR(goodVals));
       end
    else
        pkERRM = NaN;
        fittedCtERRM = NaN;
        rawKepRR = NaN;
    end
    
    %%
    stepSize = t(2)-t(1); % Assuming constant step size throughout acquisition
    
    % Initialize matrices
    [sT, sX] = size(Ct);
    CrrInt1 = stepSize*cumtrapz(Crr);
    CrrInt2 = stepSize*cumtrapz(CrrInt1);
       
    M1 = CrrInt1+kepRR*CrrInt2;
    M2 = Crr + kepRR*CrrInt1;

    pkParams = zeros(sX,3);
    fittedCt = zeros(sT,sX);

    % Solve the linear problem
    % Disabling warnings for speed (possible non-tissue regions in image give warnings)
    warning off
    
    if doNonNeg
        for i=1:sX
                curCt = squeeze(Ct(:,i));
                y = stepSize*cumtrapz(curCt);
                M3 = -stepSize*cumtrapz(y);
                M = [M1, M2, M3];
                pkParams(i,:) = lsqnonneg(M,y);
                fittedCt(:,i) = M*pkParams(i,:)';
        end
    else
        for i=1:sX
                curCt = squeeze(Ct(:,i));
                y = stepSize*cumtrapz(curCt);
                M3 = -stepSize*cumtrapz(y);
                M = [M1, M2, M3];
                pkParams(i,:) = lscov(M,y);
                fittedCt(:,i) = M*pkParams(i,:)';
        end
    end
    warning on
    
    % pkParams = [kt/ktRR + vp*kep/ktRR, vp/ktRR, kep]
    % desired = [kt/ktRR, ve/veRR, kep, vp/ktRR, kepRR]
    if doPure==false
        vpKtRR = pkParams(:,2);
        kep = pkParams(:,3);
        ktRel = pkParams(:,1) - kep .* vpKtRR;
        veRel = ktRel .* kepRR ./ kep;
        pkParams = [ktRel, veRel, kep, vpKtRR];
    end
    pkParams(:,5) = rawKepRR;
    
end
