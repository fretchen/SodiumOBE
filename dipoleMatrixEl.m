function val = dipoleMatrixEl(varargin)
%Calculates the dipole transition elements according to the conventions of
% the Steck pdf.
% By default it is assumed that the J = L+S value is 1/2
% Further we assume that the nuclear spin reads 3/2
p = inputParser;

p.addRequired('F', @(x)isnumeric(x));
p.addRequired('Fp', @(x)isnumeric(x));
p.addRequired('mF', @(x)isnumeric(x));
p.addRequired('mFp', @(x)isnumeric(x));
p.addRequired('q', @(x)isnumeric(x));
p.addOptional('J', 0.5, @(x)isnumeric(x));
p.addOptional('Jp', 0.5, @(x)isnumeric(x));
p.addParamValue('I', 1.5, @(x)isnumeric(x));
parse(p,varargin{:});

F = p.Results.F;
Fp = p.Results.Fp;
mF = p.Results.mF;
mFp = p.Results.mFp;
q = p.Results.q;
J = p.Results.J;
Jp = p.Results.Jp;
I = p.Results.I;

F_Fp_over = (-1)^(Fp+J+1+I)*sqrt((2*Fp+1)*(2*J+1))*Wigner6j(J, Jp, 1, Fp, F, I);
val = F_Fp_over*(-1)^(Fp-1+mF)*sqrt(2*F+1)*Wigner3j([Fp 1 F],[mFp, q, -mF]);
    
end
