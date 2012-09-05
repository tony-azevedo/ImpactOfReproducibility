lambda = .2;
a_bar = 1.2;
sig_d = .2;
sig_i = .1;

N = 100000;
a = poissrnd(lambda,1,N);
a = normrnd(a,sqrt(sig_d^2 + a.*sig_i^2));

hist(a,30)