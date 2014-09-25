function Coupling = constructCouplingMatrix(varargin)
p = inputParser;

p.addRequired('Elf',@(x)isnumeric(x));
p.addRequired('lowFStates', @(x)isnumeric(x));
p.addRequired('upperFStates', @(x)isnumeric(x));
p.addRequired('normEl', @(x)isnumeric(x));

p.addOptional('J', 1/2, @(x)isnumeric(x));
p.addOptional('Jp', 1/2, @(x)isnumeric(x));
p.addParamValue('I', 3/2, @(x)isnumeric(x));

parse(p,varargin{:});

Elf = p.Results.Elf;
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
                C = getCouplingForAllPolarizations(Elf, F,Fp, mF,mFp,normEl,...
                    'J', J, 'Jp', Jp, 'I', I);
                Coupling(ii,jj) = C;
                Coupling(jj,ii) = C;
            end
        end
    end
end