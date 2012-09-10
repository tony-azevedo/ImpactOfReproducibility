%% Nonlinearity Module

% Static non-linearities

% --------------------------
% saturatingExponential (out = 0 for in <0)
tau = 0.42323;
nlfunc = @(x)saturatingExponential(tau,x);
alpha = nlfunc(r_in);
% --------------------------

% --------------------------
% hill (out = 0 for in <0)
K = 0.22323;
h = 1.41001;
nlfunc = @(x)hill2([K,h],x);
alpha = nlfunc(r_in);
% --------------------------

% --------------------------
% Cum Gaussian 
mu_nl = 0.5;
sig_nl = 0.07;
nlfunc = @(x)cumulative_gauss_with_mean([sig_nl,mu_nl],x);
alpha = nlfunc(r_in);
% --------------------------


% --------------------------
% FieldRieke w(A) = P_s(A) / (P_s(A) + P_n(A))
P_s_of_A = zeros(size(r_in));
n = 1;
while sum(q==n)
    P_s_of_A = P_s_of_A + sum(q==n)/length(q) * normpdf(r_in,n*a_bar,sqrt(sig_d^2 + n*sig_i^2));
    n = n+1;
end
P_n_of_A = sum(q==0)/length(q) * normpdf(r_in,0*a_bar,sqrt(sig_d^2 + sig_i^2));
alpha = P_s_of_A./(P_s_of_A+P_n_of_A);

% --------------------------

r_out = alpha .* r_in;

a_out = r_out;


% Adaptive threshold - lets less through upon use
