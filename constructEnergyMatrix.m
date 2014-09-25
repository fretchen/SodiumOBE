function En = constructEnergyMatrix(lowFStates,upperFStates,deltaLowFStates, deltaUpperFStates)
%% constructs the energies of the Hamiltonian

NlowF = length(lowFStates);
NupperF = length(upperFStates);

NlowTotal = lowFStateEndInd(lowFStates,NlowF);
NupperTotal = upperFStateEndInd(upperFStates,NlowTotal,NupperF)-NlowTotal;

Nlevel = NlowTotal+NupperTotal;

En = zeros(Nlevel);
for lFInd=1:NlowF
    detuning = deltaLowFStates(lFInd);
    sInd = lowFStateStartInd(lowFStates,lFInd);
    eInd = lowFStateEndInd(lowFStates,lFInd);
    for ii=sInd:eInd
        En(ii,ii) = detuning;
    end
end

for uFInd=1:NupperF
    detuning = deltaUpperFStates(uFInd);
    sInd = upperFStateStartInd(upperFStates,NlowTotal,uFInd);
    eInd = upperFStateEndInd(upperFStates,NlowTotal,uFInd);
    for ii=sInd:eInd
        En(ii,ii) = detuning;
    end
end