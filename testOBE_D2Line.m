clear all;

%% Parameters for the simulation
% the electric field and the coupling strength
Elf = 4*[0 0 1];%s- pi and s+

% how long do we look at it
tmax = 20;

% we defined the zero by setting the F=1 and the F'=1 state to zero energy.
Detun = +170;% energies are measured in linewidths 

FstartInd = 2;% F state we start out with
mFstart = -2;% mF state we start out with

debug = 0;% if you want to debug things

%% define the energy structure

% the lower manifold
lowFStates = [1 2];
deltaLowFStates = [0 180];

% upper mF states
upperFStates = [0 1 2 3];
deltaUpperFStates = [-1.6 0 3.4 10] + Detun;

J = 1/2;
Jp = 3/2;

NlowF = length(lowFStates);
NupperF = length(upperFStates);

NlowTotal = lowFStateEndInd(lowFStates,NlowF);
NupperTotal = upperFStateEndInd(upperFStates,NlowTotal,NupperF)-NlowTotal;

Nlevel = NlowTotal+NupperTotal;

NManifolds = NlowF +NupperF;

% the dipole matrix element we normalize too
normEl = 1/sqrt(2);


%% construct the coupling matrix

Coupling = constructCouplingMatrix(Elf,lowFStates,upperFStates, normEl, J, Jp);
En = constructEnergyMatrix(lowFStates,upperFStates,deltaLowFStates, deltaUpperFStates);
Gamma = constructLossMatrix(lowFStates,upperFStates,normEl, J, Jp, 'debug', debug);
H = Coupling+En;

%% initialize into the right state
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
    if ~debug
        ylim([0 1]);
    end
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
    if ~debug
        ylim([0 1]);
    end
    grid on
    clear M;
end