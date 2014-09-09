function Gamma = DecayAmplitude(F,Fp,mF,mFp)

Gamma =0;
for q = -1:1
    Gamma = Gamma+dipoleMatrixEl(F,Fp,mF,mFp,q).^2;
end