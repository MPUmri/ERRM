function [pkParams, resid, fitData] = LLSQ(Ct, theAif, t, modType, doNonNeg)
  % Linear Fit, based on Murase (2004) MRM, 51(4), 858–862. http://doi.org/10.1002/mrm.20022
  % Inputs:
  % Ct [MxN] - Concentration in tissue of interest for N voxels
  % theAif [Mx1] - Concentration in blood plasma
  % t [Mx1] - Time (in minutes) at each timepoint
  % modType [bool] - Tofts Model (false or 0), Extended Tofts Model(true or 1)
  % doNonNeg [bool] - do non-negative least squares
  % Outputs:
  % pkParams = [kTrans, kep] for Tofts model
  %            [kTrans, kep, vp] for Extended Tofts Model
  % resid = residual of the fit
  % fitData = fitted curves
  
  %% Things to change
  % - The name of the function should be changed from 'LLSQ' to something
  % more appropriate (e.g. 'LinearTofts').
  % - 'theAif' should be renamed to 'Cp'
  % - Code is poorly commented
  %%
  
  if nargin<5
      doNonNeg = false;
  end

  if nargin<4
    modType = 0;
  end
  
  stepSize = t(2)-t(1); % Assuming constant step size throughout acquisition

  % Initialize matrices
  datSize=size(Ct);
  [sT sX] = size(Ct);
  fitData = zeros(sT,sX);

  if modType == 0
    % Tofts Model
    resid = zeros(sX,1);
    pkParams = zeros(sX,2);
    
    trapzCp = stepSize * cumtrapz(theAif);
    
    % Solve the linear problem
    % Disabling warnings for speed (possible non-tissue regions in image give warnings)
    warning off
    if doNonNeg
        for i=1:sX
            curCt = squeeze(Ct(:,i));
            M = zeros(sT,2);
            M(:,1) = trapzCp;
            M(:,2) = -stepSize * cumtrapz(curCt);
            pkParams(i,:) = lsqnonneg(M,curCt);
            fitData(:,i) = M*pkParams(i,:)';
            resid(i) = norm(curCt-M*pkParams(i,:)');
        end
    else
        for i=1:sX
            curCt = squeeze(Ct(:,i));
            M = zeros(sT,2);
            M(:,1) = trapzCp;
            M(:,2) = -stepSize * cumtrapz(curCt);
            pkParams(i,:) = mldivide(M,curCt);
            fitData(:,i) = M*pkParams(i,:)';
            resid(i) = norm(curCt-M*pkParams(i,:)');
        end
    end
    warning on
    
  else
    % Extended Tofts Model
    resid = zeros(sX,1);
    pkParams = zeros(sX,3);
    
    trapzCp = stepSize * cumtrapz(theAif);
    
    warning off
    % Build the matrix M
    if doNonNeg
        for i=1:sX
            curCt = squeeze(Ct(:,i));
            M = zeros(sT,2);
            M(:,1) = trapzCp;
            M(:,2) = -stepSize * cumtrapz(curCt);
            M(:,3) = theAif;
            pkParams(i,:) = lsqnonneg(M,curCt);
            fitData(:,i) = M*pkParams(i,:)';
            resid(i) = norm(curCt-M*pkParams(i,:)');
        end
    else
        for i=1:sX
            curCt = squeeze(Ct(:,i));
            M = zeros(sT,2);
            M(:,1) = trapzCp;
            M(:,2) = -stepSize * cumtrapz(curCt);
            M(:,3) = theAif;
            pkParams(i,:) = mldivide(M,curCt);
            fitData(:,i) = M*pkParams(i,:)';
            resid(i) = norm(curCt-M*pkParams(i,:)');
        end
    end
    warning on

    % Form of pkParams is [kTrans + kep*vp, kEp, vp]
    % Have to convert to [Ktrans, kep, vp]
    pkParams(:,1) = pkParams(:,1) - pkParams(:,2) .* pkParams(:,3);
  end
end
