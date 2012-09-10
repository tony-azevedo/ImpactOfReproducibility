%% Response Module

%% Responses are just a vector of amplitudes drawn from the Baylor Model

r = poissrnd(lambda,1,N);
q = r;
r = normrnd(r*a_bar,sqrt(sig_d^2 + r.*sig_i^2));
a = r;

%% Responses are simulated

%% Responses are drawn from a distribution (loaded matrix)