clear all;

% the electric field and the coupling strength

Elf = 20*[1 0 1];

% the lower manifold
lowFStates = [1 2];
deltaLowFStates = [0 180];

% the upper manifold
upperFStates = [1 2];
deltaUpperFStates = [0 18];

tmax = 5;

NlowF = length(lowFStates);
NupperF = length(upperFStates);

NlowTotal = lowFStateEndInd(lowFStates,NlowF);
NupperTotal = upperFStateEndInd(upperFStates,NlowTotal,NupperF)-NlowTotal;

Nlevel = NlowTotal+NupperTotal;

NManifolds = NlowF +NupperF;

% the dipole matrix element we normalize too
normEl = 1;

%% construct the coupling, energy and loss matrix
Coupling = constructCouplingMatrix(Elf,lowFStates,upperFStates, normEl);
En = constructEnergyMatrix(lowFStates,upperFStates,deltaLowFStates, deltaUpperFStates);
Gamma = constructLossMatrix(lowFStates,upperFStates,normEl);

H = Coupling+En;

%% run it

args = {};
A = constructDensityEvolution(0,@(t)H,args,Gamma);

[V, D] = eig(H);

int_vector = zeros(1,size(A,1));
int_vector(returnDensityEvInd(1,1,Nlevel)) = 1;

[t, y] = ode45(@(t,y)(A*y),[0 tmax],int_vector);

%% test that the populations are conserved
for jj=1:Nlevel
    pop(jj,:) = y(:,returnDensityEvInd(jj,jj,Nlevel));
end

Ntot = sum(pop,1);

%% plot it
figure(1)
plot(t,real(Ntot),'ro')

figure(2)

clf;

for lFInd =1:NlowF
    subplot(NManifolds, 1, lFInd)
    F = lowFStates(lFInd);
    sInd = lowFStateStartInd(lowFStates,lFInd);
    eInd = lowFStateEndInd(lowFStates,lFInd);
    
    ll =1;
    for jj=sInd:eInd
        ve(ll) = returnDensityEvInd(jj,jj,Nlevel);
        mF = ll-(F+1);
        if mF<0
            M(ll,:) = ['-' num2str(abs(mF))];
        else
            M(ll,:) = ['+' num2str(mF)];
        end
        ll=ll+1;
    end
    plot(t,real(y(:,ve)))
    ylabel('Populations')
    legend(M)
  %  ylim([0 1]);
    grid on
    clear M ve;
end

for uFInd =1:NupperF
    subplot(NManifolds, 1, uFInd+NlowF);
    Fp = upperFStates(uFInd);
    sFInd = upperFStateStartInd(upperFStates,NlowTotal,uFInd);
    eFInd = upperFStateEndInd(upperFStates,NlowTotal,uFInd);
    ll =1;
    for jj=sFInd:eFInd
        ve(ll) = returnDensityEvInd(jj,jj,Nlevel);
        mF = ll-(Fp+1);
        if mF<0
            M(ll,:) = ['-' num2str(abs(mF))];
        else
            M(ll,:) = ['+' num2str(mF)];
        end
        ll=ll+1;
    end
    plot(t,real(y(:,ve)))
    ylabel('Populations')
    legend(M)
%    ylim([0 1]);
    grid on
    clear M;
end