function EvolutionMatrix = constructDensityEvolution(t,Hamfun,args,DecayMatrix)

H = Hamfun(t,args{:});

numofelements = size(H);
if numofelements(1) ~= numofelements(2);
    error('Oh no!');
else
    numofelements = numofelements(1);
end

EvolutionMatrix = zeros(numofelements.^2);

% Build in the arbitrary Hamiltonian:
for ii=1:numofelements
    for jj=1:numofelements
        for kk=1:numofelements
            % First part of the commutator:
            EvolutionMatrix(returnDensityEvInd(ii,jj,numofelements),...
                            returnDensityEvInd(ii,kk,numofelements)) = ...
                EvolutionMatrix(returnDensityEvInd(ii,jj,numofelements),...
                                returnDensityEvInd(ii,kk,numofelements)) ...
                + 1i*H(kk,jj);
            % Second part of the commutator:
            EvolutionMatrix(returnDensityEvInd(ii,jj,numofelements),...
                            returnDensityEvInd(kk,jj,numofelements)) = ...
                EvolutionMatrix(returnDensityEvInd(ii,jj,numofelements),...
                                returnDensityEvInd(kk,jj,numofelements)) ...
                 - 1i*H(ii,kk);
        end
    end
end

% Looks like we still need to build in the decay matrix:
EvolutionMatrix = EvolutionMatrix + DecayMatrix;

end