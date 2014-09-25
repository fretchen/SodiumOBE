function C = getCouplingForAllPolarizations(varargin)
%% calaculate the coupling strength of the electric field between the 
%% specific states
p = inputParser;

p.addRequired('Elf',@(x)isnumeric(x));
p.addRequired('F', @(x)isnumeric(x));
p.addRequired('Fp', @(x)isnumeric(x));
p.addRequired('mF', @(x)isnumeric(x));
p.addRequired('mFp', @(x)isnumeric(x));

p.addOptional('dMEnorm', 1, @(x)isnumeric(x));
p.addParamValue('J', 1/2, @(x)isnumeric(x));
p.addParamValue('Jp', 1/2, @(x)isnumeric(x));
p.addParamValue('I', 3/2, @(x)isnumeric(x));

parse(p,varargin{:});

C = 0;
Elf = p.Results.Elf;
F = p.Results.F;
Fp = p.Results.Fp;
mF = p.Results.mF;
mFp = p.Results.mFp;

dMEnorm = p.Results.dMEnorm;
J = p.Results.J;
Jp = p.Results.Jp;
I = p.Results.I;

for q = -1:1
    C = C+ Elf(q+2)*dipoleMatrixEl(F,Fp,mF,mFp,q, J, Jp, 'I', I)./dMEnorm;
end
