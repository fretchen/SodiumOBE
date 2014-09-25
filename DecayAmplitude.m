function Gamma = DecayAmplitude(varargin)
p = inputParser;

p.addRequired('F', @(x)isnumeric(x));
p.addRequired('Fp', @(x)isnumeric(x));
p.addRequired('mF', @(x)isnumeric(x));
p.addRequired('mFp', @(x)isnumeric(x));

p.addOptional('dMEnorm', 1, @(x)isnumeric(x));
p.addParamValue('J', 1/2, @(x)isnumeric(x));
p.addParamValue('Jp', 1/2, @(x)isnumeric(x));
p.addParamValue('I', 3/2, @(x)isnumeric(x));

parse(p,varargin{:});

F = p.Results.F;
Fp = p.Results.Fp;
mF = p.Results.mF;
mFp = p.Results.mFp;

dMEnorm = p.Results.dMEnorm;
J = p.Results.J;
Jp = p.Results.Jp;
I = p.Results.I;

Gamma =0;
for q = -1:1
    Gamma = Gamma+dipoleMatrixEl(F,Fp,mF,mFp,q, J, Jp, 'I', I).^2./dMEnorm.^2;
end