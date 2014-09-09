function C = getCouplingForAllPolarizations(Elf, F,Fp, mF,mFp)
C = 0;
for q = -1:1
    C = C+ Elf(q+2)*dipoleMatrixEl(F,Fp,mF,mFp,q);
end
