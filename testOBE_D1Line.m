clear all;

%% Parameters for the simulation
% the electric field and the coupling strength
Elf = 4*[0 1 0];%s+ pi and s-
%Elf = 4*[1 0 1];%s+ pi and s-

% how long do we look at it
tmax = 20;

% we defined the zero by setting the F=1 and the F'=1 state to zero energy.
Detun = -18;% energies are measured in linewidths 

% F state we start out with
FstartInd = 1;
mFstart = -1;

debug = 0;
%% define the energy structure

% the lower manifold
lowFStates = [1 2];
deltaLowFStates = [0 180];

%
upperFStates = [1 2];
deltaUpperFStates = [0 18] + Detun;

J = 1/2;
Jp =1/2;

% the dipole matrix element we normalize too
normEl = 1;

%% construct the coupling matrix

Coupling = constructCouplingMatrix(Elf,lowFStates,upperFStates, normEl, J, Jp);
En = constructEnergyMatrix(lowFStates,upperFStates,deltaLowFStates, deltaUpperFStates);
Gamma = constructLossMatrix(lowFStates,upperFStates,normEl, J, Jp);

H = Coupling+En;

%% initialize into the right state
NlowF = length(lowFStates);
NupperF = length(upperFStates);

NlowTotal = lowFStateEndInd(lowFStates,NlowF);
NupperTotal = upperFStateEndInd(upperFStates,NlowTotal,NupperF)-NlowTotal;

Nlevel = NlowTotal+NupperTotal;

NManifolds = NlowF +NupperF;

startInd = lowFStateStartInd(lowFStates,FstartInd);
Fstart = lowFStates(FstartInd);
jjStart = mFstart+(startInd-1)+(Fstart+1);

%% run it

args = {};
A = constructDensityEvolution(0,@(t)H,args,Gamma);

int_vector = zeros(1,size(A,1));
int_vector(returnDensityEvInd(jjStart,jjStart,Nlevel)) = 1;

[V, D] = eig(H);

[t, y] = ode45(@(t,y)(A*y),[0 tmax],int_vector);

%% test that the populations are conserved
if debug
    for jj=1:Nlevel
        pop(jj,:) = y(:,returnDensityEvInd(jj,jj,Nlevel));
    end
    
    Ntot = sum(pop,1);
    figure(1)
    clf;
    plot(t,real(Ntot),'ro')
end

%% plot it

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
    ylim([0 1]);
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
    ylim([0 1]);
    grid on
    clear M;
end