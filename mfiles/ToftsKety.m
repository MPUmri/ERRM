function ct = ToftsKety(Cp, pkParams, t, doExt)

	if nargin < 4
		doExt = 0;
    end

    if isrow(Cp)
        Cp = Cp';
    end
    
    kTrans = pkParams(1);
    kep = pkParams(2);
    
    ct = kTrans * expConv(Cp,kep,t);
    
    if doExt == 1
    	ct = ct + pkParams(3)*Cp;
    end
end

