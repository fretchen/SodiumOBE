function sInd = upperFStateStartInd(upperFStates,NlowStates,uFInd)

NStates = 2*upperFStates+1;
sInd = NlowStates+1;
if uFInd > 1;
    for ii = 2:uFInd
        sInd = sInd+NStates(ii-1);
    end
end

