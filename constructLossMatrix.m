function Gamma = constructLossMatrix(varargin)
%% constructs the loss matrix for the OBE
p = inputParser;

p.addRequired('lowFStates', @(x)isnumeric(x));
p.addRequired('upperFStates', @(x)isnumeric(x));
p.addRequired('normEl', @(x)isnumeric(x));

p.addOptional('J', 1/2, @(x)isnumeric(x));
p.addOptional('Jp', 1/2, @(x)isnumeric(x));
p.addParamValue('I', 3/2, @(x)isnumeric(x));

parse(p,varargin{:});

lowFStates = p.Results.lowFStates;
upperFStates = p.Results.upperFStates;
normEl = p.Results.normEl;

J = p.Results.J;
Jp = p.Results.Jp;
I = p.Results.I;

NlowF = length(lowFStates);
NupperF = length(upperFStates);

NlowTotal = lowFStateEndInd(lowFStates,NlowF);
NupperTotal = upperFStateEndInd(upperFStates,NlowTotal,NupperF)-NlowTotal;

Nlevel = NlowTotal+NupperTotal;

Gamma = zeros(Nlevel^2);

% upper lifetime
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
                amp = DecayAmplitude(F,Fp,mF,mFp,normEl,...
                    'J', J, 'Jp', Jp, 'I', I);
                if amp
                    Gamma(lGind,uGind) = amp;
                    Gamma(uGind,uGind) = Gamma(uGind,uGind) - amp;
                end
            end
        end
    end
end

sum(sum(Gamma))

% now the lifetime of the coherences between excited and lower state
for lFInd=1:NlowF
    sInd = lowFStateStartInd(lowFStates,lFInd);
    eInd = lowFStateEndInd(lowFStates,lFInd);
    for uFInd = 1:NupperF
        sFInd = upperFStateStartInd(upperFStates,NlowTotal,uFInd);
        eFInd = upperFStateEndInd(upperFStates,NlowTotal,uFInd);
        for ii=sInd:eInd
            for jj = sFInd:eFInd
                invLT = -Gamma(returnDensityEvInd(jj,jj,Nlevel),returnDensityEvInd(jj,jj,Nlevel));
                Gamma(returnDensityEvInd(ii,jj,Nlevel),returnDensityEvInd(ii,jj,Nlevel)) =-invLT/2;
                Gamma(returnDensityEvInd(jj,ii,Nlevel),returnDensityEvInd(jj,ii,Nlevel)) =-invLT/2;
                
            end
        end
    end
end