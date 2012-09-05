%% Response Module - generate responses

% --------------------------
% Responses are just a vector of amplitudes drawn from the Baylor Model
lambda = .2;
a_bar = 1.2;
sig_d = 0.12;
sig_i = 0.2;

N = 100000;
r = poissrnd(lambda,1,N);
q = r;
r = normrnd(r,sqrt(sig_d^2 + r.*sig_i^2));
a = r;
% --------------------------

% --------------------------
% Responses are simulated

% --------------------------

% --------------------------
% Responses are drawn from a distribution

% --------------------------

r_in = r;

% Filtering Module

% If r is a matrix of responses, apply filter

% Nonlinearity Module

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


r_out = alpha .* r_in;

a_out = r_out;


% Adaptive threshold - lets less through upon use

% Evaluation module

% ROC 
t = q >= 1;
a_l = min([a,a_out]);
a_r = max([a,a_out]);

a_norm = (a - a_l)/(a_r-a_l);
a_out_norm = (a_out - a_l)/(a_r-a_l);
[tpr_a,fpr_a,th_a] = roc(t,a_norm);
[tpr_a_out,fpr_a_out,th_a_out] = roc(t,a_out_norm);


% Plotting Module

% amplitude transform
figure(1), clf
[bins,deltax] = optimalBinWidth(a);
[n,bins] = hist(a,bins);
n = n/length(a)/deltax;
stairs(bins-deltax/2,n,'Color',[1 1 1] * sig_i);

hold on;

xlims = get(gca,'xlim');
ylims = get(gca,'ylim');

[bins,deltax] = optimalBinWidth(a_out);
[n,bins] = hist(a_out,bins);
n = n/length(a_out)/deltax;
stairs(bins-deltax/2,n,'Color',[1 1 1] * sig_i + [0 0 1]*(1-sig_i));

line(bins,nlfunc(bins),...
	'parent',gca,...
	'linestyle','--',...
	'color',[0 0 0],...
	'Tag','tag',...
	'DisplayName','legendName');

set(gca,'xlim',xlims);
set(gca,'ylim',ylims);
set(gca,'yScale','log')

%% ROC
figure(2),clf
line(fpr_a_out,tpr_a_out,...
	'linestyle','-',...
	'linewidth',2,...
	'color',[0 0 1],...
	'parent',gca,...
	'Tag','tag',...
	'DisplayName','legendName');

line(fpr_a_out,tpr_a_out,...
	'linestyle','-',...
	'color',[1 0 0],...
	'parent',gca,...
	'Tag','tag',...
	'DisplayName','legendName');


%% Pooling - i.e. Use output from above to add to the next layer