clear all;

% the electric field and the coupling strength

Elf = 20*[0 0 1];

% the lower manifold
lowFStates = [1 2];
deltaLowFStates = [0 180];

% the upper manifold on the D2 line
%{
upperFStates = [0 1 2 3];
deltaUpperFStates = [0 1.5 5 11];
%}

% the upper manifold on the D1 line
upperFStates = [1 2];
deltaUpperFStates = [0 18];

tmax = 5;

NlowF = length(lowFStates);
NupperF = length(upperFStates);

NlowTotal = lowFStateEndInd(lowFStates,NlowF);
NupperTotal = upperFStateEndInd(upperFStates,NlowTotal,NupperF)-NlowTotal;

Nlevel = NlowTotal+NupperTotal;

NManifolds = NlowF +NupperF;

%% construct the coupling matrix
Coupling = zeros(Nlevel);

for lFInd=1:NlowF
    F = lowFStates(lFInd);
    sInd = lowFStateStartInd(lowFStates,lFInd);
    eInd = lowFStateEndInd(lowFStates,lFInd);
    for uFInd = 1:NupperF
        Fp = upperFStates(uFInd);
        sFInd = upperFStateStartInd(upperFStates,NlowTotal,uFInd);
        eFInd = upperFStateEndInd(upperFStates,NlowTotal,uFInd);
        for ii=sInd:eInd
            mF = ii-(sInd-1)-(F+1);
            for jj = sFInd:eFInd
                mFp = jj-(sFInd-1)-(Fp+1);
                C = getCouplingForAllPolarizations(Elf, F,Fp, mF,mFp);
                Coupling(ii,jj) = C;
                Coupling(jj,ii) = C;
            end
        end
    end
end

%% Construct the diagonal matrix
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

H = Coupling+En;

%% construct the loss matrix
Gamma = zeros(Nlevel^2);

% upper lifetime
invLT = 1;

for lFInd=1:NlowF
    F = lowFStates(lFInd);
    sInd = lowFStateStartInd(lowFStates,lFInd);
    eInd = lowFStateEndInd(lowFStates,lFInd);
    for uFInd = 1:NupperF
        Fp = upperFStates(uFInd);
        sFInd = upperFStateStartInd(upperFStates,NlowTotal,uFInd);
        eFInd = upperFStateEndInd(upperFStates,NlowTotal,uFInd);
        for ii=sInd:eInd
            mF = ii-(sInd-1)-(F+1);
            for jj = sFInd:eFInd
                mFp = jj-(sFInd-1)-(Fp+1);               
                lGind = returnDensityEvInd(ii,ii,Nlevel);
                uGind = returnDensityEvInd(jj,jj,Nlevel);
                Gamma(lGind,uGind) = DecayAmplitude(F,Fp,mF,mFp);
                Gamma(uGind,uGind) = Gamma(uGind,uGind) - DecayAmplitude(F,Fp,mF,mFp);
            end
        end
    end
end
sum(sum(Gamma))

% normalize the lifetime
for uFInd = 1:NupperF
    sFInd = upperFStateStartInd(upperFStates,NlowTotal,uFInd);
    eFInd = upperFStateEndInd(upperFStates,NlowTotal,uFInd);
    uGind = returnDensityEvInd(sFInd,sFInd,Nlevel);
    uFlt = Gamma(uGind,uGind);
    if uFlt
        for jj=sFInd:eFInd
            uGind = returnDensityEvInd(jj,jj,Nlevel);
            Gamma(uGind,uGind) = Gamma(uGind,uGind)/abs(uFlt)*invLT;
            for lFInd=1:NlowF
                sInd = lowFStateStartInd(lowFStates,lFInd);
                eInd = lowFStateEndInd(lowFStates,lFInd);
                for ii=sInd:eInd
                    lGind = returnDensityEvInd(ii,ii,Nlevel);
                    Gamma(lGind,uGind) = Gamma(lGind,uGind)/abs(uFlt)*invLT;
                end
            end
        end
    end
end

sum(sum(Gamma))

% now the lifetime of the coherences between excited and lower state
for lFInd=1:NlowF
    F = lowFStates(lFInd);
    sInd = lowFStateStartInd(lowFStates,lFInd);
    eInd = lowFStateEndInd(lowFStates,lFInd);
    for uFInd = 1:NupperF
        Fp = upperFStates(uFInd);
        sFInd = upperFStateStartInd(upperFStates,NlowTotal,uFInd);
        eFInd = upperFStateEndInd(upperFStates,NlowTotal,uFInd);
        for ii=sInd:eInd
            for jj = sFInd:eFInd
                Gamma(returnDensityEvInd(ii,jj,Nlevel),returnDensityEvInd(ii,jj,Nlevel)) =-invLT/2;
                
            end
        end
    end
end

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
figure(1)
plot(t,Ntot,'ro')

figure(2)
clf;
%% plot it
figure(1)
plot(t,Ntot,'ro')
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
    plot(t,y(:,ve))
    ylabel('Populations')
    legend(M)
    ylim([0 1]);
    grid on
    clear M;
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
    plot(t,y(:,ve))
    ylabel('Populations')
    legend(M)
    ylim([0 1]);
    grid on
    clear M;
end