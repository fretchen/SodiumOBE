function Gamma = constructLossMatrix(varargin)
%% constructs the loss matrix for the OBE
%% !!!!! be careful the lifetime of the upper states is supposed to be one
p = inputParser;

p.addRequired('lowFStates', @(x)isnumeric(x));
p.addRequired('upperFStates', @(x)isnumeric(x));
p.addRequired('normEl', @(x)isnumeric(x));

p.addOptional('J', 1/2, @(x)isnumeric(x));
p.addOptional('Jp', 1/2, @(x)isnumeric(x));
p.addParamValue('I', 3/2, @(x)isnumeric(x));
p.addParamValue('debug', 0, @(x)isnumeric(x));

parse(p,varargin{:});

lowFStates = p.Results.lowFStates;
upperFStates = p.Results.upperFStates;
normEl = p.Results.normEl;

J = p.Results.J;
Jp = p.Results.Jp;
I = p.Results.I;
debug = p.Results.debug;

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

% now the lifetime of the coherences
for ii=1:Nlevel
    for jj = 1:Nlevel
        if ii-jj
            ind = returnDensityEvInd(ii,jj,Nlevel);
            Gamma(ind,ind) =-1/2;
        end
    end
end

for uFInd = 1:NupperF
    sFInd = upperFStateStartInd(upperFStates,NlowTotal,uFInd);
    
    invLT = -Gamma(returnDensityEvInd(sFInd,sFInd,Nlevel),returnDensityEvInd(sFInd,sFInd,Nlevel));
    if abs(invLT-1)>1e-12
        disp(invLT-1)
        error('The upper lifetime has to be equal too one for the moment.');
    end
end

if debug
    for uFInd = 1:NupperF
        sFInd = upperFStateStartInd(upperFStates,NlowTotal,uFInd);
        
        invLT = -Gamma(returnDensityEvInd(sFInd,sFInd,Nlevel),returnDensityEvInd(sFInd,sFInd,Nlevel));
        disp([ 'invLT =' num2str(invLT)]);
    end
end
