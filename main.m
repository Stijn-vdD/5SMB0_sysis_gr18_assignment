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
% figure(50); clf;
% ir_est = impulseest(ze, 30);
% impulseplot(ir_est);
% title('Impulse response estimate (check if h(0) \approx 0 for nk=1)');
% grid minor;
% If h(0) is significantly nonzero, nk=0 is appropriate (direct feedthrough).
% If h(0) ≈ 0, nk=1 (one-sample delay) may be more appropriate.

% Sweep BJ orders under assignment constraints:
% nb >= 4, nf > 2, nd >= 2, nc >= 2, nk <= 1.
% Selection strategy:
% 1) Keep models with plant fit within 2% of the best fit.
% 2) Within that set, select the model with best residual diagnostics
%    (fewest R_eu and R_ee bound violations), then AIC.
% nb_range = 4:6;
% nc_range = 2:5;
% nd_range = 2:10;
% nf_range = 3:10;
% nk_range = 0:1;
% nb_range = 3:7;
% % nf_range = 4:10;
% nf_range = 5:9; %[2 4 6 8]
% nc_range = 1:2; %:2
% nd_range = 2:3; %:3
% nk_range = 0;

% nb_range = 5;
% % nf_range = 4:10;
% nf_range = 5:6; %[2 4 6 8]
% nc_range = 1:5; %:2
% nd_range = 1:5; %:3
% nk_range = 0;

% nb_range = 5;
% % nf_range = 4:10;
% nf_range = 4:5; %[2 4 6 8]
% nc_range = 1:2; %:2
% nd_range = 2:3; %:3
% nk_range = 0;

% nb_range = 5;
% % nf_range = 4:10;
% % nf_range = 4:5; %[2 4 6 8]
% nf_range = 4:5; %[2 4 6 8]
% nc_range = 2:10; %:2
% nd_range = 2:10; %:3
% nk_range = 0;

nb_range = 5;
% nf_range = 4:10;
% nf_range = 4:5; %[2 4 6 8]
nf_range = 4:5; %[2 4 6 8]
nc_range = 2; %:2
nd_range = 5; %:3
nk_range = 0;

n_comb = numel(nb_range) * numel(nc_range) * numel(nd_range) * ...
    numel(nf_range) * numel(nk_range);
orders_sweep = zeros(n_comb, 5);
fit_sweep = -inf(n_comb, 1);
aic_sweep = inf(n_comb, 1);
n_eu_out_sweep = inf(n_comb, 1);
n_ee_out_sweep = inf(n_comb, 1);
frac_ee_out_sweep = inf(n_comb, 1);
noise_order_sweep = inf(n_comb, 1); % nc + nd

N_sweep = size(zv.OutputData, 1);
max_lag_sweep = min(50, N_sweep - 1);
conf99_sweep = 2.576 / sqrt(N_sweep);

k = 0;
for nb = nb_range
    for nc = nc_range
        for nd = nd_range
            for nf = nf_range
                for nk = nk_range
                    k = k + 1;
                    orders_try = [nb, nc, nd, nf, nk];
                    orders_sweep(k, :) = orders_try;
                    try
                        sys_try = bj(ze, orders_try);
                        [~, fit_tmp] = compare(zv, sys_try);
                        fit_sweep(k) = mean(fit_tmp(:));
                        aic_sweep(k) = aic(sys_try);

                        % Independent error-model quality metrics from residuals
                        e_try_id = pe(zv, sys_try);
                        e_try = e_try_id.OutputData(:,1);
                        u_try = zv.InputData(:,1);
                        [Reu_try, ~] = xcorr(e_try, u_try, max_lag_sweep, 'coeff');
                        [Ree_try, lags_ee_try] = xcorr(e_try, e_try, max_lag_sweep, 'coeff');
                        n_eu_out_sweep(k) = sum(abs(Reu_try) > conf99_sweep);
                        n_ee_out_sweep(k) = sum(abs(Ree_try(lags_ee_try ~= 0)) > conf99_sweep);
                        frac_ee_out_sweep(k) = n_ee_out_sweep(k) / (numel(Ree_try) - 1);
                        noise_order_sweep(k) = orders_try(2) + orders_try(3);
                    catch
                        % Keep defaults for failed fits.
                    end
                end
            end
        end
    end
end

valid_idx = isfinite(aic_sweep) & isfinite(fit_sweep) & ...
    isfinite(n_eu_out_sweep) & isfinite(n_ee_out_sweep) & ...
    isfinite(noise_order_sweep);
if ~any(valid_idx)
    error('No valid BJ model found in the sweep ranges.');
end

whiteness_pass_idx = valid_idx & (n_ee_out_sweep == 0);

if any(whiteness_pass_idx)
    best_fit_pass = max(fit_sweep(whiteness_pass_idx));
    fit_eps = 1e-9;
    best_fit_rows = find(whiteness_pass_idx & (abs(fit_sweep - best_fit_pass) <= fit_eps));

    % Requested criterion among whiteness-passing models:
    % highest plant fit, then lowest noise-model order (nc+nd), then AIC.
    sel_table = [noise_order_sweep(best_fit_rows), aic_sweep(best_fit_rows)];
    [~, sel_sort_idx] = sortrows(sel_table, [1 2]);
    best_idx = best_fit_rows(sel_sort_idx(1));
    selection_mode = 'Whiteness-pass constrained';
else
    % Fallback when no model passes whiteness: maximize fit, then minimize
    % whiteness violations and noise-model order.
    best_fit_valid = max(fit_sweep(valid_idx));
    fit_eps = 1e-9;
    best_fit_rows = find(valid_idx & (abs(fit_sweep - best_fit_valid) <= fit_eps));
    sel_table = [n_ee_out_sweep(best_fit_rows), noise_order_sweep(best_fit_rows), aic_sweep(best_fit_rows)];
    [~, sel_sort_idx] = sortrows(sel_table, [1 2 3]);
    best_idx = best_fit_rows(sel_sort_idx(1));
    selection_mode = 'No whiteness-pass fallback';
end

optimal_orders = orders_sweep(best_idx, :);
sys_bj = bj(ze, optimal_orders);

fprintf('\n--- BJ order sweep (Section 3.2) ---\n');
fprintf('Best orders [nb nc nd nf nk] = [%d %d %d %d %d]\n', optimal_orders);
fprintf('Validation fit = %.2f%%, AIC = %.4g\n', fit_sweep(best_idx), aic_sweep(best_idx));
fprintf('R_eu outside 99%% bounds = %d, R_ee outside 99%% bounds = %d\n', ...
    n_eu_out_sweep(best_idx), n_ee_out_sweep(best_idx));
fprintf('Noise-model order (nc+nd) = %d\n', noise_order_sweep(best_idx));
fprintf('Whiteness-pass candidates: %d / %d\n', sum(whiteness_pass_idx), sum(valid_idx));
fprintf('Selection mode: %s\n', selection_mode);

score_table = [orders_sweep(valid_idx, :), fit_sweep(valid_idx), aic_sweep(valid_idx), ...
    n_eu_out_sweep(valid_idx), n_ee_out_sweep(valid_idx), frac_ee_out_sweep(valid_idx), ...
    noise_order_sweep(valid_idx), whiteness_pass_idx(valid_idx)];
[~, sort_idx] = sortrows(score_table, [-12, -6, 11, 7]);
disp(['All valid BJ candidates (sorted by passWhiteness desc, fit desc, ', ...
    'noiseOrder asc, AIC asc):']);
disp(array2table(score_table(sort_idx, :), ...
    'VariableNames', {'nb','nc','nd','nf','nk','fitPct','AIC', ...
    'nReuOut99','nReeOut99','fracReeOut99','noiseOrderNcNd','passWhiteness99'}));

% Compare chosen nk with the alternate nk (0 <-> 1), holding other orders fixed.
figure(51); clf;
alt_orders = optimal_orders;
alt_orders(5) = 1 - optimal_orders(5);
try
    sys_bj_altnk = bj(ze, alt_orders);
    compare(zv, sys_bj, sys_bj_altnk);
    legend('Validation Data', ...
        sprintf('BJ best [%d %d %d %d %d]', optimal_orders), ...
        sprintf('BJ alt nk [%d %d %d %d %d]', alt_orders));
catch
    compare(zv, sys_bj);
    legend('Validation Data', sprintf('BJ best [%d %d %d %d %d]', optimal_orders));
end
grid minor;
title('Delay comparison for best BJ structure');

% For a consistent estimate, Reu (cross-correlation) must be within intervals
figure(5); clf;
resid(zv, sys_bj);

% Explicit correlation diagnostics for report figure:
% - R_eu(tau): consistency check of plant model G0
% - R_ee(tau): whiteness test of residuals (noise model quality)
e_id = pe(zv, sys_bj);
e_val = e_id.OutputData;
u_val = zv.InputData;
N_res = size(e_val, 1);
max_lag = min(50, N_res - 1);
conf99 = 2.576 / sqrt(N_res);

[Reu, lags_eu] = xcorr(e_val(:,1), u_val(:,1), max_lag, 'coeff');
[Ree, lags_ee] = xcorr(e_val(:,1), e_val(:,1), max_lag, 'coeff');

figure(52); clf;
tiledlayout(2,1, 'TileSpacing', 'compact');
nexttile;
stem(lags_eu, Reu, 'filled');
hold on;
yline(conf99, 'r--', '99% bound');
yline(-conf99, 'r--');
yline(0, 'k-');
grid minor;
xlabel('Lag \tau [samples]');
ylabel('R_{\epsilon u}(\tau)');
title(sprintf('Residual-input cross-correlation (BJ [%d %d %d %d %d])', optimal_orders));

nexttile;
stem(lags_ee, Ree, 'filled');
hold on;
yline(conf99, 'r--', '99% bound');
yline(-conf99, 'r--');
yline(0, 'k-');
grid minor;
xlabel('Lag \tau [samples]');
ylabel('R_{\epsilon\epsilon}(\tau)');
title('Residual autocorrelation (whiteness test)');

n_eu_out = sum(abs(Reu) > conf99);
n_ee_out = sum(abs(Ree(lags_ee ~= 0)) > conf99);
fprintf('\n--- Residual correlation diagnostics (Part 3) ---\n');
fprintf('R_eu outside 99%% bounds: %d / %d lags\n', n_eu_out, numel(Reu));
fprintf('R_ee outside 99%% bounds (excluding tau=0): %d / %d lags\n', ...
    n_ee_out, numel(Ree)-1);
if n_eu_out == 0
    fprintf('Plant consistency test: PASS (no significant input-residual correlation).\n');
else
    fprintf('Plant consistency test: FAIL/WEAK (remaining input-residual correlation).\n');
end
if n_ee_out == 0
    fprintf('Residual whiteness test: PASS (noise model appears adequate).\n');
else
    fprintf('Residual whiteness test: FAIL/WEAK (residuals not fully white).\n');
end

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
% Figure 52 (R_eu and R_ee with 99% bounds) is used for the consistency
% and whiteness assessment in Part 3.3.
present(sys_bj);

%% Part 4: Experimental verification of variance estimates
%%% 4.1 Monte Carlo simulations
% Repeat Question 3 for 100 times: each run uses the same PRBS reference
% but a new noise realization (from assignment_sys_18) and estimates a BJ model.
n_mc = 100;
orders = optimal_orders; % [nb nc nd nf nk] selected in 3.2
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
    % Use the toolbox-compatible BJ syntax for each Monte Carlo run.
    sys_mc = bj(ze_mc, orders);
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
P_cov = getcov(sys_bj);
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

%% Part 5: MIMO identification
%%% 5.1 Parametric identification of a 2x2 MIMO system G0
N_mimo = 3000;

% Input design: two uncorrelated PRBS signals (one per input channel).
% Generate both channels in one idinput call for MATLAB-version compatibility.
% Independent channels are critical for separating the contribution of each input to each
% output. Uncorrelated inputs prevent ill-conditioning in the regressor
% matrix (analogous to the requirement that Phi_uu is full rank at all
% frequencies for MIMO identifiability).
% Band = [0 1] excites all frequencies up to Nyquist, amplitude within
% saturation bounds [-M, M].
r5 = idinput([N_mimo 2], 'prbs', [0 1], [-M M]); % size (N_mimo, 2)

[u5, y5] = assignment_sys_18(r5, 'MIMO');

data_mimo = iddata(y5, u5, 1, 'Domain', 'Time');

% Plot input/output signals
figure(12); clf;
tiledlayout(2, 2, 'TileSpacing', 'compact');
for ch = 1:2
    nexttile;
    plot(u5(:, ch));
    ylabel(sprintf('u_%d', ch));
    grid minor;
end
for ch = 1:2
    nexttile;
    plot(y5(:, ch));
    ylabel(sprintf('y_%d', ch));
    grid minor;
end
xlabel('Samples');

% 50/50 split for estimation and validation
ze_mimo = detrend(data_mimo(1:N_mimo/2), 0);
zv_mimo = detrend(data_mimo(N_mimo/2+1:end), 0);

%%% Model structure choice: state-space via subspace identification (n4sid)
% For MIMO systems, state-space models are more natural than polynomial
% (BJ/OE) models because:
% 1. A single state-space model captures all input-output channels
%    simultaneously with shared state dynamics.
% 2. Polynomial MIMO models require specifying separate orders for each
%    transfer function entry (4 B polynomials + 4 F polynomials for a 2x2),
%    which is cumbersome and prone to over-parameterization.
% 3. Subspace methods (n4sid) provide a non-iterative, numerically robust
%    initial estimate without local minima issues.

% Determine model order using singular value analysis
% figure(13); clf;
n4sid(ze_mimo);
title('Singular values for MIMO model order selection');
% Inspect the singular value plot: look for a clear gap indicating the
% appropriate model order. Select the order where the singular values
% drop significantly.

% From the SISO identification (Part 3), we expect:
% - 2 complex pole pairs (resonance + anti-resonance dynamics) per channel
%   → order ~4-8 for the plant dynamics
% Try a range of orders and compare validation fit
n_orders = 2:12;
sys_mimo_candidates = cell(length(n_orders), 1);
fit_scores = zeros(length(n_orders), 1);

for idx = 1:length(n_orders)
    sys_mimo_candidates{idx} = n4sid(ze_mimo, n_orders(idx), 'Focus', 'sim');
end

% Compare all candidates on validation data
figure(14); clf;
compare(zv_mimo, sys_mimo_candidates{:});
legend(['Validation Data', compose('n=%d', n_orders)]);
grid minor;
title('MIMO model order comparison');

% Select the best order based on validation fit (best compare %)
% Use a parsimony criterion: pick the lowest order that achieves a fit
% close to the best (within ~2-3% of the maximum)
for idx = 1:length(n_orders)
    [~, fit_tmp] = compare(zv_mimo, sys_mimo_candidates{idx});
    fit_scores(idx) = mean(fit_tmp); % average fit across output channels
end
fprintf('\n--- MIMO model order selection ---\n');
for idx = 1:length(n_orders)
    fprintf('  n=%2d: avg fit = %.1f%%\n', n_orders(idx), fit_scores(idx));
end

[best_fit, best_idx] = max(fit_scores);
% Pick lowest order within 3% of the best fit (parsimony)
threshold = best_fit - 3;
chosen_idx = find(fit_scores >= threshold, 1, 'first');
n_chosen = n_orders(chosen_idx);
fprintf('Best fit: n=%d (%.1f%%), chosen (parsimony): n=%d (%.1f%%)\n', ...
    n_orders(best_idx), best_fit, n_chosen, fit_scores(chosen_idx));

sys_mimo_n4sid = sys_mimo_candidates{chosen_idx};

%%% Refine with ssest (prediction error minimization)
% n4sid provides a consistent but not necessarily efficient estimate.
% ssest refines it by minimizing the prediction error (PEM), using the
% n4sid result as initial condition to avoid local minima.
sys_mimo = ssest(ze_mimo, sys_mimo_n4sid);

% Validation: compare on validation data
figure(15); clf;
compare(zv_mimo, sys_mimo_n4sid, sys_mimo);
legend('Validation Data', sprintf('n4sid (n=%d)', n_chosen), ...
    sprintf('ssest (n=%d)', n_chosen));
grid minor;
title('MIMO: n4sid vs ssest refinement');

% Residual analysis for consistency check
figure(16); clf;
resid(zv_mimo, sys_mimo);
title('MIMO residual analysis');

% Bode plot of the identified 2x2 transfer matrix
figure(17); clf;
bode_mimo = bodeplot(sys_mimo);
showConfidence(bode_mimo);
grid minor;
title('Identified MIMO system G_0(q)');

% Pole-zero map
figure(18); clf;
pzmap(sys_mimo);
grid minor;
title('MIMO pole-zero map');

fprintf('\n--- MIMO identification summary ---\n');
fprintf('Model order: n = %d states\n', n_chosen);
fprintf('Method: n4sid (subspace) + ssest (PEM refinement)\n');
present(sys_mimo);

