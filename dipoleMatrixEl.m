function val = dipoleMatrixEl(F,Fp,mF,mFp,q)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
val = -1^(Fp-1+mF)*sqrt(2*F+1)*Wigner3j([Fp 1 F],[mFp, q, -mF]);
    
end
