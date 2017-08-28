function [ ro,p ] = getCorrCoef(A,B)

    [R,P] = corrcoef(A,B);
    
    ro=R(1,2);
    p=P(1,2);
end

