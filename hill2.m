function y = hill2(beta, x)
% see also hill

x(x<0) = 0;
y = 1.0 ./ (1 + (beta(1) ./ x).^beta(2));
