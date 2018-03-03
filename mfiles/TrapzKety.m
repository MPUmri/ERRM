function ct = TrapzKety(theAif, pkParams, t, doExt)
    % Alternative to ToftsKety()

	if nargin < 4
		doExt = 0;
	end

    kTrans = pkParams(1);
    kep = pkParams(2);
    stepSize = t(2)-t(1);
    
    ct = zeros(size(theAif));
    
    for k=1:length(t)
        ct(k) = kTrans*trapz(theAif(1:k).*exp(kep*(t(1:k)-t(k))));
    end
    
    ct = stepSize * ct';
    
    if doExt == 1
    	ct = ct + pkParams(3)*theAif';
    end
end

