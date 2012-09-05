%% Response Module

%% Responses are just a vector of amplitudes drawn from the Baylor Model
lambda = .2;
a_bar = 1.2;
sig_d = .2;
sig_i = 0.5572;

N = 100000;
r = poissrnd(lambda,1,N);
r = normrnd(r,sqrt(sig_d^2 + r.*sig_i^2));

% if
n = hist(r,optimalBinWidth(r));
plot(optimalBinWidth(r),n,'color', [1 1 1] * sig_i);

%% Responses are simulated

%% Responses are drawn from a distribution (loaded matrix)