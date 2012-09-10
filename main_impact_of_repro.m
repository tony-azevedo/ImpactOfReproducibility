%% Response Module - generate responses

% --------------------------
% Responses are just a vector of amplitudes drawn from the Baylor Model
lambda = .2;
sig_d = 0.2;
a_bar = 1;
CV = 0.24;
sig_i = CV*a_bar;
N = 100000;


Rods = 2;
resp_mod = cell(Rods,4);

for i = 1:Rods;
    responseModule_Amp
    resp_mod(i,:) = {r_in,r,q,a};
end

a_ex = a;
r_ex = r_in;
q_ex = q;

% --------------------------

% --------------------------
% Responses are simulated

% --------------------------

% --------------------------
% Responses are drawn from data

% --------------------------


%% Filtering Module

% If r is a matrix of responses, apply filter

%% Nonlinearity Module

nl_mod = cell(Rods,2);

for i = 1:Rods;
    [r_in,r,q,a] = resp_mod{i,:};
    nonlinearityModule
    nl_mod(i,:) = {r_out,a_out};
end


%% Pooling - i.e. Use output from above to add to the next layer

r_out = zeros(size(resp_mod{1,1}));
a = r_out;
a_out = r_out;
q_out = r_out;

for i = 1:Rods;
    % linear
    [r_l,r,q,a_l] = resp_mod{i,:};
    a = a + a_l;    
    q_out = q_out + q;
    
    % nonlinear
    [r_nl,a_nl] = nl_mod{i,:};
    r_out = r_out + r_nl;
    a_out = a_out + a_nl;
end


%% Evaluation module

% ROC 
t = q_out >= 1;
a_l = min([a,a_out]);
a_r = max([a,a_out]);

a_norm = (a - a_l)/(a_r-a_l);
a_out_norm = (a_out - a_l)/(a_r-a_l);
[tpr_a,fpr_a,th_a] = roc(t,a_norm);
[tpr_a_out,fpr_a_out,th_a_out] = roc(t,a_out_norm);


%% Plotting Module

% amplitude transform
figure(1), clf

subplot(2,1,1)
[bins,deltax] = optimalBinWidth(a_ex);
[n,bins] = hist(a_ex,bins);
n = n/length(a_ex)/deltax;
stairs(bins-deltax/2,n,'Color',[1 1 1] * sig_i);

hold on;

% show the convolution of distributions
for i = 1:Rods;
    [n,bins] = hist(a_ex,bins);
end

xlims = get(gca,'xlim');
ylims = get(gca,'ylim');

%[bins,deltax] = optimalBinWidth(a_out);
[n,bins] = hist(a_out,bins);
n = n/length(a_out)/deltax;
stairs(bins-deltax/2,n,'Color',[1 1 1] * sig_i + [0 0 1]*(1-sig_i));

line(bins,nlfunc(bins),...
	'parent',gca,...
	'linestyle','--',...
	'color',[0 0 0],...
	'Tag','tag',...
	'DisplayName','legendName');

set(gca,'xlim',[-1.5,6]);
set(gca,'ylim',ylims);
%set(gca,'yScale','log')


subplot(2,1,2)
[bins,deltax] = optimalBinWidth(a);
[n,bins] = hist(a,bins);
n = n/length(a)/deltax;
stairs(bins-deltax/2,n,'Color',[1 1 1] * sig_i);

hold on;

xlims = get(gca,'xlim');
ylims = get(gca,'ylim');

%[bins,deltax] = optimalBinWidth(a_out);
[n,bins] = hist(a_out,bins);
n = n/length(a_out)/deltax;
stairs(bins-deltax/2,n,'Color',[1 1 1] * sig_i + [0 0 1]*(1-sig_i));

line(bins,nlfunc(bins),...
	'parent',gca,...
	'linestyle','--',...
	'color',[0 0 0],...
	'Tag','tag',...
	'DisplayName','legendName');

set(gca,'xlim',[-1.5,6]);
set(gca,'ylim',ylims);
%set(gca,'yScale','log')


% subplot(2,1,2)
% line(fpr_a_out,tpr_a_out,...
% 	'linestyle','-',...
% 	'linewidth',2,...
% 	'color',[0 0 1],...
% 	'parent',gca,...
% 	'Tag','tag',...
% 	'DisplayName','legendName');
% 
% line(fpr_a_out,tpr_a_out,...
% 	'linestyle','-',...
% 	'color',[1 0 0],...
% 	'parent',gca,...
% 	'Tag','tag',...
% 	'DisplayName','legendName');







