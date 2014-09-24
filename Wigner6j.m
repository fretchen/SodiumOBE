function Wigner = Wigner6j(j1,j2,j3,J1,J2,J3)
% Wigner 6j-symbol calculator. Written by Amita B Deb, Clarendon Lab. 2007. 
% Improved by Richard A. Holt, Univ. of Western Ontario, 2009. 
% Further improved to deal with large input arguments by Lee Harper, 
% CRL Oxford, 2009. 

% Builds on the m-file Wigner6jcoeff posted by Richard A. Holt. 

% Calculates { j1, j2 ,j3}  using Racah formula. See: Sobelman: Atomic Spectra and Radiative Transitions.
%              J1  J2  J3

% Finding Triangular coefficients

tri1 = triangle_coeff(j1,j2,j3);
tri2 = triangle_coeff(j1,J2,J3);
tri3 = triangle_coeff(J1,j2,J3);
tri4 = triangle_coeff(J1,J2,j3);

if (tri1==0||tri2==0||tri3==0||tri4==0)
    Wigner=0;
    return
end

% Finding the range of summation in the Racah formula.

a(1) = j1 + j2 + j3;
a(2) = j1 + J2 + J3;
a(3) = J1 + j2 + J3;
a(4) = J1 + J2 + j3;

rangei = max(a);

k(1) = j1 + j2 + J1 + J2;
k(2) = j2 + j3 + J2 + J3;
k(3) = j3 + j1 + J3 + J1;

rangef = min(k);

Wigner = 0;

   for t=rangei:rangef
   
      Wigner = Wigner + ((-1)^t)*exp(gammaln(t+2) - fung(t,j1,j2,j3,J1,J2,J3));
      
   end
   
   Wigner = (tri1*tri2*tri3*tri4)^(0.5)*Wigner;
   
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Subfunctions:

function r = fung(t, j1,j2,j3,J1,J2,J3)
% Calculating the logarithm of the denominator in Racah Formula, using
% the gamma function in place of the factorial.
r = sum(gammaln([(t-j1-j2-j3);(t-j1-J2-J3);(t-J1-j2-J3);(t-J1-J2-j3);...
                 (j1+j2+J1+J2-t);(j2+j3+J2+J3-t);(j3+j1+J3+J1-t)] + 1));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tri = triangle_coeff(a,b,c)
% Calculates triangle coefficients for angular momenta.
% This version returns 0 if the triangle inequalities are violated.  (RAH)

if (a<0 || b<0 || c<0) 
    tri=0;
    return
end

for xa = abs(a-b):1:a+b
    if c==xa
        tri = factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/(factorial(a+b+c+1));
        return
    end
end

tri=0;

end

   