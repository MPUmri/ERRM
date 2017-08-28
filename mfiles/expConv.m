function [f] = expConv(A,B,t)
    % Does conv(A,exp(-B*t))
    % Based on Flouri et al. (2016) MRM 76(3), doi: 10.1002/mrm.25991
        
    n = length(t);
	
    x = B * ( t(2:n) - t(1:n-1) );
	dA = ( A(2:n) - A(1:n-1) ) ./ x;

	E = exp(-x);
	E0 = 1 - E;
	E1 = x - E0;

	iterAdd = A(1:n-1).*E0 + dA.*E1;
    %%
    f = zeros(n,1);
    for i=1:n-1
       f(i+1) = E(i)*f(i) + iterAdd(i); 
    end
    f = f ./ B;
end

%% Old version

