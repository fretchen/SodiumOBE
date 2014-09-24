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

parse(p,varargin{:});

C = 0;
Elf = p.Results.Elf;
F = p.Results.F;
Fp = p.Results.Fp;
mF = p.Results.mF;
mFp = p.Results.mFp;

dMEnorm = p.Results.dMEnorm;

for q = -1:1
    C = C+ Elf(q+2)*dipoleMatrixEl(F,Fp,mF,mFp,q)./dMEnorm;
end
