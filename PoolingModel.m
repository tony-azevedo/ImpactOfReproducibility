%%
% Model including non-discrete noise integration time and post-rod noise
% source

th = 2; % minimum number of quanta required for passing the nonlinearity
Ac = 512; % Effective rod convergence, number of rods, [rod] - power of 2

%model bg light and noise

Ibg = 0.0; % additional bg light, [R*/rod/sec]
ID =  0.003; %equivalent dark noise, [R*/rod/sec]
It = 0.08; % integration time, [s]
SynapticNoiseSD = 0.6;

% Frequency of seeing parameters
FOSTh = 10;          % threshold
FOSAc = 16;          % in units of Ac - i.e. FOS based on this many nonlinear subunits - power of 2
FOSIf = [0 1 3 10 30 100 300] * 1e-5;    % flash strengths

If = [2.06784E-4 4.35568E-4 8.99896E-4 0.00184 0.00375 0.00765]; % flash strengths,  [R*/rod]

% PART 1: create distribution of amplitudes of thermal and background events
% assuming exponential weighting

% time domain weights
tme = 0:0.001:1;
TimeDomainFilter = exp(-tme/It);
AmpStep = 0.025;

ScFact = ceil(max(If)*Ac*10);
DistX = -ScFact:AmpStep:ScFact;
SingleRodNoiseDist = hist(TimeDomainFilter, DistX);
SingleRodNoiseDist = SingleRodNoiseDist / (sum(SingleRodNoiseDist) * AmpStep);
Indices = find(DistX >= 0);
SingleRodNoiseDist = SingleRodNoiseDist * exp(-((Ibg+ID))) * (Ibg+ID);
SingleRodNoiseDist(Indices(1)) = SingleRodNoiseDist(Indices(1)) + exp(-(Ibg+ID)) / AmpStep;
SingleRodNoiseDist = SingleRodNoiseDist / (sum(SingleRodNoiseDist) * AmpStep);
PooledRodNoiseDist = SingleRodNoiseDist;

for rod = 1:round(log(Ac) / log(2))
    PooledRodNoiseDist = conv(PooledRodNoiseDist, PooledRodNoiseDist, 'same');
    PooledRodNoiseDist = PooledRodNoiseDist / (sum(PooledRodNoiseDist) * AmpStep);
end

% add post-pooling noise
SynapticNoiseDist = normpdf(DistX, 0, SynapticNoiseSD);
PooledRodNoiseDist = conv(PooledRodNoiseDist, SynapticNoiseDist, 'same');
PooledRodNoiseDist = PooledRodNoiseDist / (sum(PooledRodNoiseDist) * AmpStep);

VarLin = sum(PooledRodNoiseDist .* (DistX.^2) * AmpStep) - (sum(PooledRodNoiseDist .* DistX) * AmpStep).^2;
TempDist = PooledRodNoiseDist;
TempDist(1:Indices(1) + th/AmpStep) = 0;
TempDist(Indices(1)) = sum(PooledRodNoiseDist(1:Indices(1) + th/AmpStep));
VarNL = sum(TempDist .* (DistX.^2) * AmpStep) - (sum(TempDist .* DistX) * AmpStep).^2;

fprintf(1, 'noise reduction by nonlinearity = %d\n', VarLin/VarNL);

% PART 2: add noise distribution to flash distribution and compute mean
% responses
% of linear and nonlinear models
clear ResponseDist linear lin;
for flash = 1:length(If)
    SignalDist = zeros(1, length(DistX));
    for photon = 1:If(flash)*Ac*10+2
        SignalDist(Indices(1) + (photon-1)/AmpStep) = poisspdf(photon-1, If(flash)*Ac);
    end
    ResponseDist(flash, :) = conv(SignalDist, PooledRodNoiseDist, 'same');
    ResponseDist(flash, :) = ResponseDist(flash, :) / (sum(ResponseDist(flash, :)) * AmpStep);
    linear(flash) = (sum(DistX .* ResponseDist(flash, :)) - sum(DistX .* PooledRodNoiseDist)) * AmpStep / If(flash);
    lin(flash) = (sum(DistX(Indices(1) + th/AmpStep:length(DistX)) .* ResponseDist(flash, Indices(1) + th/AmpStep:length(DistX))) - sum(DistX(Indices(1) + th/AmpStep:length(DistX)) .* PooledRodNoiseDist(Indices(1) + th/AmpStep:length(DistX)))) * AmpStep / If(flash); 
end

lin = lin / max(linear);
linear = linear / max(linear);

FS = [0.28386 0.23836 0.50977 0.81547 0.99997 0.82412]; %flash sensitivity as measured in darkness, example cell

%083110 example cell, 0.008 R*/rod/sec bg light:

If2 = [0.0004 0.0009 0.0018 0.0038 0.0077]; % intensities from the exp
FS2 = [0.93808 1.00002 0.91115 0.93691 0.74348]; %measured flash sensitivity

% plot the models and the data, flash sensitivity

figure(1); clf; subplot(2, 2, 1);
plot(If, FS,'ko'), hold on                 % plot data in darkness, example cell
plot(If2, FS2,'go'), hold on               % plot data in background light, example cell   
plot(If, linear,'g', 'linewidth', 1)       % plot model prediction in darkness
plot(If, lin,'k', 'linewidth', 1)          % plot linear prediction        
axis tight
% let's plot the data and models in log scales

subplot(2, 2, 2);
loglog(If,FS.*If,'ko'), hold
loglog(If2,FS2.*If2,'go')
plot(If, linear.*If,'g', 'linewidth', 1)
plot(If, lin.*If,'k', 'linewidth', 1)
axis tight

% PART 3: FOS (relies on above for noise dist)

% compute full subunit response distribution
for flash = 1:length(FOSIf)
    SignalDist = zeros(1, length(DistX));
    for photon = 1:FOSIf(flash)*Ac*10+2
        SignalDist(Indices(1) + (photon-1)/AmpStep) = poisspdf(photon-1, FOSIf(flash)*Ac);
    end
    ResponseDist(flash, :) = conv(SignalDist, PooledRodNoiseDist, 'same');
    ResponseDist(flash, :) = ResponseDist(flash, :) / (sum(ResponseDist(flash, :)) * AmpStep);
end

% linear or nonlinear full RF dist
clear LinResponseDist NLResponseDist linFOS nlFOS;
for flash = 1:length(FOSIf)
    LinResponseDist(flash, :) = ResponseDist(flash, :);
    for subunit = 1:round(log(FOSAc) / log(2))
        LinResponseDist(flash, :) = conv(LinResponseDist(flash, :), LinResponseDist(flash, :), 'same');
    end
    NLResponseDist(flash, :) = ResponseDist(flash, :);
    NLResponseDist(flash, 1:Indices(1) + th/AmpStep) = 0;
    NLResponseDist(flash, Indices(1)) = sum(ResponseDist(flash, 1:Indices(1) + th/AmpStep));
    for subunit = 1:round(log(FOSAc) / log(2))
        NLResponseDist(flash, :) = conv(NLResponseDist(flash, :), NLResponseDist(flash, :), 'same');
    end
    LinResponseDist(flash, :) = LinResponseDist(flash, :) / (sum(LinResponseDist(flash, :)) * AmpStep);
    NLResponseDist(flash, :) = NLResponseDist(flash, :) / (sum(NLResponseDist(flash, :)) * AmpStep);
    linFOS(flash) = sum(LinResponseDist(flash, Indices(1) + FOSTh/AmpStep:length(DistX))) * AmpStep;
    nlFOS(flash) = sum(NLResponseDist(flash, Indices(1) + FOSTh/AmpStep:length(DistX))) * AmpStep;    
end

subplot(2, 2, 3);
semilogx(FOSIf, linFOS, 'g'); hold on
semilogx(FOSIf, nlFOS, 'k');
axis tight

fprintf(1, 'FP rate %d %d\n', linFOS(1), nlFOS(1));

% 2AFC
subplot(2, 2, 4)
clear LinPCorrect NLPCorrect LinPCorrectReject NLPCorrectReject
for flash = 1:length(FOSIf)
    LRatio = LinResponseDist(flash, :) ./ LinResponseDist(1, :);
    TempIndices = find(LRatio > 1);
    TempEqualIndices = find(LRatio == 1);
    LinPCorrect(flash) = (sum(LinResponseDist(flash, TempIndices)) + sum(LinResponseDist(flash, TempEqualIndices))/2) * AmpStep;
    TempIndices = find(LRatio < 1);
    LinPCorrectReject(flash) = (sum(LinResponseDist(1, TempIndices)) + sum(LinResponseDist(1, TempEqualIndices))/2) * AmpStep;
    LRatio = NLResponseDist(flash, :) ./ NLResponseDist(1, :);
    TempIndices = find(LRatio > 1);
    TempEqualIndices = find(LRatio == 1);
    NLPCorrect(flash) = (sum(NLResponseDist(flash, TempIndices)) + sum(NLResponseDist(flash, TempEqualIndices))/2) * AmpStep;
    TempIndices = find(LRatio < 1);
    NLPCorrectReject(flash) = (sum(NLResponseDist(1, TempIndices)) + sum(NLResponseDist(1, TempEqualIndices))/2) * AmpStep;
end
plot(FOSIf, (LinPCorrect+LinPCorrectReject)/2, 'og', FOSIf, (NLPCorrect+NLPCorrectReject)/2, 'ok');
axis tight

%%
% SIMPLE MODEL FOR NONLINEARITY

% entirely discrete model - assumes binary window for noise integration

th = 2; % minimum number of quanta required for passing the nonlinearity
Ac = 500; % Effective rod convergence, number of rods, [rod]

%model bg light 

Ibg = 0.0; % additional bg light, [R*/rod/sec]
ID =  0.003; %equivalent dark noise, [R*/rod/sec]
It = 0.110; % integration time, [sec]

IBg = Ac.*It.*(ID + Ibg); % background intensity, mean rate, [R*]

If = [2.06784E-4 4.35568E-4 8.99896E-4 0.00184 0.00375 0.00765]; % flash strengths,  [R*/rod]


M = Ac.* If; % ints in R* per Ac, [R*] - R* per RGC
% MBg = ones(1,length(M)).* IBg; % mean R* from BG light
Mtot = M + IBg; % mean , [R*] 

linear = ones(1, length(Mtot)) - poisscdf(0, Mtot); % the expected loss for a linear system, ie. the fraction of failures based on Poisson distribution
lin = ones(1, length(Mtot)) - poisscdf(th-1, Mtot); % the expected loss for the nonlinear system with a threshold, e.g. if the=2, all doubles pass

    
FS = [0.28386 0.23836 0.50977 0.81547 0.99997 0.82412]; %flash sensitivity as measured in darkness, example cell

%083110 exampe cell, 0.008 R*/rod/sec bg light:

If2 = [0.0004 0.0009 0.0018 0.0038 0.0077]; % intensities from the exp
FS2 = [0.93808 1.00002 0.91115 0.93691 0.74348]; %measured flash sensitivity

% let's plot the models and the data, flash sensitivity

figure(1); clf
plot(If, FS,'ko'), hold on                 % plot data in darkness, example cell
plot(If2, FS2,'go'), hold on               % plot data in background light, example cell   
plot(If, linear,'g', 'linewidth', 1)       % plot model prediction in darkness
plot(If, lin,'k', 'linewidth', 1)          % plot linear prediction        

% let's plot the data and models in log scales

figure(2); clf
loglog(If,FS.*If,'ko'), hold
loglog(If2,FS2.*If2,'go')
plot(If, linear.*If,'g', 'linewidth', 1)
plot(If, lin.*If,'k', 'linewidth', 1)




