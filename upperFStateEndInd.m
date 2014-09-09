function eInd = upperFStateEndInd(upperFStates,NlowStates,uFInd)

NStates = 2*upperFStates+1;
eInd = NlowStates;
for ii = 1:uFInd
        eInd = eInd+NStates(ii);
end

