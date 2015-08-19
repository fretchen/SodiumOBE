function sInd = lowFStateStartInd(lowFStates,lFInd)

NStates = 2*lowFStates+1;
if lFInd == 1;
    sInd = 1;
else
    sInd = 1;
    for ii = 2:lFInd
        sInd = sInd+NStates(ii-1);
    end
end

