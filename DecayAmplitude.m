function Gamma = DecayAmplitude(varargin)
p = inputParser;

p.addRequired('F', @(x)isnumeric(x));
p.addRequired('Fp', @(x)isnumeric(x));
p.addRequired('mF', @(x)isnumeric(x));
p.addRequired('mFp', @(x)isnumeric(x));

p.addOptional('dMEnorm', 1, @(x)isnumeric(x));

parse(p,varargin{:});

F = p.Results.F;
Fp = p.Results.Fp;
mF = p.Results.mF;
mFp = p.Results.mFp;

dMEnorm = p.Results.dMEnorm;

Gamma =0;
for q = -1:1
    Gamma = Gamma+dipoleMatrixEl(F,Fp,mF,mFp,q).^2./dMEnorm.^2;
end