clearvars;
% close all;
restoredefaultpath;

q = tf('q');
F = (0.7157 + 1.4315*q^-1 + 0.7157*q^-2)/(1 + 1.3490*q^-1 + 0.5140*q^-2);

%% Part 1: Understanding saturation and Butterworth filter
%%% 1.1 Determine saturation parameter M by experiment
% 100 second input sweep from -5 to 5
r1 = linspace(-5,5,100);

[u1,~] = assignment_sys_18(r1,'open loop');

M = max(u1); % == -1 * min(u)

% The combination of actuator saturation and Butterworth filter limits the
% amount of energy that can be applied to the system in the high frequency
% range. This will reduce the signal-to-noise ratio in that range.

%% Part 2: Nonparametric identification
%%% 2.1 Frequency behavior of G_0
% Reference design:
% - Within saturation bounds
% - Open loop FRF estimate: Schroeder phase multisine signal excites all 
%   frequencies with deterministic amplitude spectrum and good crest factor
% - Periodic signal to prevent leakage

N = 1500;   % Experiment length
n = 10;      % Number of periods
if mod(N,n) ~= 0
    error("n is not a clean divisor of N, reference signal cannot be periodic.")
end
p = N/n;    % Period length

r2 = repmat(M * multisine([0,0.49],1,p),1,n);

[u2,y2] = assignment_sys_18(r2,'open loop');

frfdata = iddata(y2,u2,1,Domain='Time',Period=p);

% Check for transient effects
u2last = u2((n-1)*p+1:n*p);
y2last = y2((n-1)*p+1:n*p);
figure(1);clf;
tiledlayout(TileSpacing='compact');
nexttile;
for i = 1:n
    u2plot = u2((i-1)*p+1:i*p) - u2last;
    plot((i-1)*p+1:i*p,u2plot);
    hold on;
end
ylabel('u');
nexttile;
for i = 1:n
    y2plot = y2((i-1)*p+1:i*p) - y2last; % Remove steady-state offset
    plot((i-1)*p+1:i*p,y2plot);
    hold on;
end
ylabel('y');
xlabel('t');

% Result: transient negligible, input affected for first ~20 samples

%%% 2.2 Identify system FRF
Ghat_spa = spa(frfdata);

f2 = figure(2);clf;
Ghat_bode = bodeplot(Ghat_spa);
Ghat_bode.PhaseWrappingEnabled = true;
showConfidence(Ghat_bode);
legend;
grid minor;

% Result:
% - Anti-resonance around 1.2 rad/s (two complex zeros)
% - Resonance around 2.03 rad/s (two complex poles)
% - Anti-resonance around 2.63 rad/s (two complex zeros)
% - DC gain -50 dB


%%% 2.3 Plot the magnitude of the estimated noise power spectrum
figure(3);clf;

error_spectrum = spectrumplot(Ghat_spa);
showConfidence(error_spectrum);
grid minor;

% The noise spectrum is concentrated around 2 rad/s

%% Part 3: Parametric identification and validation
%%% 3.1 Parametric model of G0
N_new = 3000; % New experiment length

% Generate a PRBS signal
% PRBS maximizes signal energy for a given amplitude bound
r3 = idinput(N_new, 'prbs', [0 1], [-M M]);

[u3, y3] = assignment_sys_18(r3, 'open loop');
data_param = iddata(y3, u3, 1, 'Domain', 'Time');

figure(4); clf;
tiledlayout(2,1, 'TileSpacing', 'compact');
nexttile;
plot(u3);
ylabel('u');
grid minor;
nexttile;
plot(y3);
ylabel('y');
xlabel('Samples');
grid minor;

% 50/50 split
ze = data_param(1:N_new/2);     % Estimation data
zv = data_param(N_new/2+1:end); % Validation data
ze = detrend(ze, 0);
zv = detrend(zv, 0);

%%% 3.2: Consistent parametric identification of G0
% Used some sweeps to find the best fit, phase still wrong though at 1.2rad/s:
% [nb nc nd nf nk]
% - nb=5 (4 zeros) captures the two anti-resonances.
% - nf=8 (8 poles) captures the resonance, the known Butterworth filter, 
%   and high-frequency roll-off, increased it to improve the fit (done by the sweep).
% - nc=5, nd=2 captures the colored noise spectrum concentrated around 2 rad/s,
%   determined by the sweep.
% - nk=0 accounts for the direct feedthrough observed in cross-correlation.

optimal_orders = [5, 5, 2, 8, 0];
sys_bj = bj(ze, optimal_orders);

% For a consistent estimate, Reu (cross-correlation) must be within intervals
figure(5); clf;
resid(zv, sys_bj);

% Time-domain cross validatipon
figure(6); clf;
compare(zv, sys_bj);
grid minor;

% Comparison with 2
figure(7); clf;
bode_spabj = bodeplot(Ghat_spa, 'b', sys_bj, 'r');
bode_spabj.PhaseWrappingEnabled = true;
showConfidence(bode_spabj);
legend('Nonparametric (SPA)', 'Parametric (BJ)');
grid minor;
% Chat specified that:
% If the BJ model passes the R_eu residual test but fails R_e, you have a consistent 
%   estimate of G0, but the noise model (nc, nd) might need tweaking.

%%% 3.3 Minimum variance estimate
% To prove whether the model achieves minimum variance, we can look at the 
% confidence regions of the estimated poles and zeros. 
% Over-parameterization leads to pole-zero cancellations and inflated variance.
present(sys_bj);

optimal_orders_reduced = [5, 5, 2, 4, 0]; % Changed nf from 8 to 4 (variance larger than parameter)
sys_bj_reduced = bj(ze, optimal_orders_reduced);
% Maintain fit
figure(8); clf;
compare(zv, sys_bj, sys_bj_reduced);
legend('Validation Data', 'BJ [5 5 2 8 0]', 'BJ Reduced [5 5 2 4 0]');
grid minor;
% Pole zero map with confidences
figure(9); clf;
hold on;
h_pz = iopzplot(sys_bj, sys_bj_reduced);
showConfidence(h_pz, 3);
legend('BJ [5 5 2 8 0]', 'BJ Reduced [5 5 2 4 0]');
grid minor;
axis equal;
% Bode to check, yay it also fixed the phase at 1.3rad/s
figure(10); clf;
bode_spabjreduced = bodeplot(Ghat_spa, 'b', sys_bj, 'r', sys_bj_reduced, 'g');
bode_spabjreduced.PhaseWrappingEnabled = true;
showConfidence(bode_spabjreduced);
legend('Nonparametric (SPA)', 'Parametric (BJ)', 'Parametric (BJ reduced)');
grid minor;
% Check parameter variance
present(sys_bj_reduced);