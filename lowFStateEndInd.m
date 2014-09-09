function eInd = lowFStateEndInd(lowFStates,lFInd)

NStates = 2*lowFStates+1;
eInd = 0;
for ii = 1:lFInd
        eInd = eInd+NStates(ii);
end

