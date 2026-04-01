clearvars;
% close all;
restoredefaultpath;
rng(42); % Fix random seed for reproducibility

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
% BJ model: y = (B/F)u + (C/D)e — separates plant and noise model,
% so noise model misspecification does not bias the G0 estimate (consistency).
% We identify from u to y, so F(q) and S(·) are NOT part of G0.
%
% [nb nc nd nf nk]
% - nb=5 (4 zeros) captures the two anti-resonances seen in the FRF.
% - nf=8 (8 poles) captures the resonance, high-frequency roll-off, and
%   additional dynamics; started high to ensure a good fit, reduced in 3.3.
% - nc=5, nd=2 captures the colored noise spectrum concentrated around 2 rad/s.
% - nk=0 assumes direct feedthrough (no pure delay from u to y).
%   Verify with impulse response / cross-correlation:

% Verify nk choice: plot impulse response estimate to check for delay
figure(50); clf;
ir_est = impulseest(ze, 30);
impulseplot(ir_est);
title('Impulse response estimate (check if h(0) \approx 0 for nk=1)');
grid minor;
% If h(0) is significantly nonzero, nk=0 is appropriate (direct feedthrough).
% If h(0) ≈ 0, nk=1 (one-sample delay) may be more appropriate.

optimal_orders = [5, 5, 2, 8, 0];
sys_bj = bj(ze, optimal_orders);

% Also test nk=1 for comparison
optimal_orders_nk1 = [5, 5, 2, 8, 1];
sys_bj_nk1 = bj(ze, optimal_orders_nk1);
figure(51); clf;
compare(zv, sys_bj, sys_bj_nk1);
legend('Validation Data', 'BJ nk=0', 'BJ nk=1');
grid minor;
title('Delay comparison: nk=0 vs nk=1');
% Result: h(0) ≈ 0.005 (≈ 0), suggesting nk=1 is physically correct.
% However, nk=0 gives 83% vs 75% fit — the extra b_0 parameter improves
% the fit significantly. We keep nk=0 for the better fit; note in the
% report that the true delay is likely 1 sample.

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
% If the BJ model passes the R_eu residual test but fails R_e, G0 is
% consistently estimated, but the noise model (nc, nd) may need tweaking.

%%% 3.3 Minimum variance estimate
% To prove whether the model achieves minimum variance, we can look at the 
% confidence regions of the estimated poles and zeros. 
% Over-parameterization leads to pole-zero cancellations and inflated variance.
present(sys_bj);

% Reduce nf from 8 to 4: present(sys_bj) shows that several F coefficients
% have standard deviation larger than the parameter magnitude, indicating
% pole-zero near-cancellations and over-parameterization.
% The FRF suggests ~2 complex poles (resonance at 2.03 rad/s) → nf=4.
optimal_orders_reduced = [5, 5, 2, 4, 0];
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

%% Part 4: Experimental verification of variance estimates
%%% 4.1 Monte Carlo simulations
% Repeat Question 3 for 100 times: each run uses the same PRBS reference
% but a new noise realization (from assignment_sys_18) and estimates a BJ model.
n_mc = 100;
orders = optimal_orders_reduced; % [nb nc nd nf nk] = [5 5 2 4 0]
n_B = orders(1);
n_F = orders(4);

% Use getpvec to extract parameters — guarantees same ordering as getcov
% Parameter order: [B(n_B), C(n_C), D(n_D), F(n_F)]
n_total = sum(orders(1:4)); % total free parameters
params_all = zeros(n_mc, n_total);
cov_diags = zeros(n_mc, n_total); % store getcov diagonal from each run

% Fix the PRBS across all MC runs so that only the noise realization varies.
% This matches the getcov assumption of a fixed input design.
r_mc = idinput(N_new, 'prbs', [0 1], [-M M]);

for i = 1:n_mc
    % Same reference signal, new noise realization from assignment_sys_18
    [u_mc, y_mc] = assignment_sys_18(r_mc, 'open loop');
    data_mc = iddata(y_mc, u_mc, 1, 'Domain', 'Time');
    ze_mc = detrend(data_mc(1:N_new/2), 0);
    % Warm-start from Part 3 model to reduce local-minima effects
    sys_mc = bj(ze_mc, orders, sys_bj_reduced);
    params_all(i,:) = getpvec(sys_mc)';
    cov_diags(i,:) = diag(getcov(sys_mc))';
end

% Extract B and F columns from the parameter matrix
idx_B = 1:n_B;
idx_F = (n_total - n_F + 1):n_total;
params_B = params_all(:, idx_B);
params_F = params_all(:, idx_F);

% Plot parameter distributions
figure(11); clf;
subplot(2,1,1);
boxplot(params_B, 'Labels', compose('b_%d', 0:n_B-1));
ylabel('Value');
title('B(q) coefficients over 100 MC runs');
grid minor;
subplot(2,1,2);
boxplot(params_F, 'Labels', compose('f_%d', 1:n_F));
ylabel('Value');
title('F(q) coefficients over 100 MC runs');
grid minor;

% The parameters vary across runs because each experiment has a different
% noise realization e(t). Since the noise enters the output, each dataset
% yields a slightly different estimate. The variance depends on the
% signal-to-noise ratio and the experiment length.

%%% 4.2 Theoretical variance from one experiment
% getcov returns the covariance matrix with the same parameter ordering
% as getpvec: [B, C, D, F]
P_cov = getcov(sys_bj_reduced);
var_theo = diag(P_cov);
var_B_theo = var_theo(idx_B);
var_F_theo = var_theo(idx_F);

fprintf('\n--- Theoretical variances (from getcov) ---\n');
fprintf('B coefficients: '); fprintf('%.6e  ', var_B_theo); fprintf('\n');
fprintf('F coefficients: '); fprintf('%.6e  ', var_F_theo); fprintf('\n');

% The theoretical covariance is accurate when:
% 1. The model structure is correct (contains the true system)
% 2. The number of data points N is large (asymptotic result)
% 3. The noise model is correctly specified
% Under these conditions, the Cramér-Rao lower bound is achieved.

%%% 4.3 Compare Monte Carlo variance with theoretical variance
var_B_mc = var(params_B);
var_F_mc = var(params_F);
mean_B_mc = mean(params_B);
mean_F_mc = mean(params_F);

fprintf('\n--- Monte Carlo statistics (100 runs) ---\n');
fprintf('B mean:     '); fprintf('%.6f  ', mean_B_mc); fprintf('\n');
fprintf('B variance: '); fprintf('%.6e  ', var_B_mc); fprintf('\n');
fprintf('F mean:     '); fprintf('%.6f  ', mean_F_mc); fprintf('\n');
fprintf('F variance: '); fprintf('%.6e  ', var_F_mc); fprintf('\n');

fprintf('\n--- Comparison: MC variance / Theoretical variance (single experiment) ---\n');
fprintf('B ratios: '); fprintf('%.3f  ', var_B_mc ./ var_B_theo'); fprintf('\n');
fprintf('F ratios: '); fprintf('%.3f  ', var_F_mc ./ var_F_theo'); fprintf('\n');

% Diagnostic: average theoretical variance across all MC runs
% If this matches MC variance better, then any single getcov is just noisy.
% If it still doesn't match, the model structure causes systematic bias.
avg_cov_diag = mean(cov_diags);
var_B_theo_avg = avg_cov_diag(idx_B);
var_F_theo_avg = avg_cov_diag(idx_F);

fprintf('\n--- Comparison: MC variance / Average theoretical variance (over %d runs) ---\n', n_mc);
fprintf('B ratios: '); fprintf('%.3f  ', var_B_mc ./ var_B_theo_avg); fprintf('\n');
fprintf('F ratios: '); fprintf('%.3f  ', var_F_mc ./ var_F_theo_avg); fprintf('\n');

% The theoretical covariance (getcov) is an asymptotic (N->inf) Cramer-Rao
% lower bound based on a linearization of the prediction error.
% Remaining discrepancies between MC variance and theoretical variance
% are due to:
% 1. Finite sample effects (N=1500 estimation samples)
% 2. F parameters enter nonlinearly (as 1/F) in the BJ prediction error,
%    making the linearized Cramer-Rao bound a loose lower bound for finite N
% 3. Possible model misspecification if nf=4 does not fully capture G0's
%    denominator dynamics
% B parameters enter linearly, so their theoretical variance is tighter.

