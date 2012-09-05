function y = saturatingExponential(coef, x)
% exponentialSaturating     
%       simple exponential fit: 1-e^-(x/tau)
%
%   exponentialSaturating(STARTINGK,X)
%
% See also NLINFIT
%
% Help added by TA 09052012
x(x<0) = 0;
y = 1 - exp(-x / coef(1));

