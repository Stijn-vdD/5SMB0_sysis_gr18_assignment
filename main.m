clearvars;
close all;
restoredefaultpath;
rng(42); % Fix random seed for reproducibility
plotting_header;

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


% Check for transient effects
u2last = u2((n-1)*p+1:n*p);
y2last = y2((n-1)*p+1:n*p);

f_transient = figure('units', 'centimeters', ...
    'Position', [5, 5, fig_setup.fig_wd, fig_setup.fig_hgt*0.7], ...
    'Name', 'Part 2: Transient Detection');clf;
tiledlayout(2,1,TileSpacing='tight');
nexttile;hold on;box on;
for i = 1:n-1
    u2plot = log(abs( u2((i-1)*p+1:i*p) - u2last ));
    plot((i-1)*p+1:i*p,u2plot);
end
xline(1*p, ':k');
yline(log(0.001),'-k','0.1\% error','Interpreter','latex','LineWidth',1.25,'LabelVerticalAlignment','bottom');
yline(log(eps),'-k','Precision floor','Interpreter','latex','LineWidth',1.25);
ylabel('$\log \left| u-u_n \right|$','Interpreter','latex');
xlim([0 N]);
ylim([-40 0]);
xticklabels([]);

nexttile;hold on;box on;
for i = 1:n-1
    y2plot = log(abs( y2((i-1)*p+1:i*p) - y2last ));
    plot((i-1)*p+1:i*p,y2plot,'Color',color(mod(i-1,7)+1));
    xline(i*p, ':k');
end
ylabel('$\log \left| y-y_n \right|$','Interpreter','latex');
xlabel('Time [samples]','Interpreter','latex');
xlim([0 N-p]);
ylim([-10 0]);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fig_setup.fntsize);

if fig_setup.export_figures
    filename = "output-figures/transient detection" + fig_setup.img_ext;
    exportgraphics(f_transient, filename, ...
        'ContentType', fig_setup.img_format, ...
        'BackgroundColor', 'none');
    fprintf('Figure exported: %s\n', filename);
end

% Result: transient negligible, input affected for first ~20 samples
frfdata = iddata(y2,u2,1,Domain='Time',Period=p);

%%% 2.2 Identify system FRF
w = 2*pi/p:2*pi/p:pi;
Ghat_spa = spa(frfdata,p,w);
Ghat_etfe = etfe(frfdata);


f_bode_etfe = figure('units', 'centimeters', ...
    'Position', [5, 5, fig_setup.fig_wd, fig_setup.fig_hgt*0.7], ...
    'Name', 'Part 2: Bode plot of ETFE model');

bodeopts = bodeoptions;
bodeopts.XLim = [1e-1, 1e1];
bodeopts.PhaseWrapping = 'on';
bodeopts.Title.String = '';
bodeopts.Grid = 'on';
bodeopts.TickLabel.FontSize = fig_setup.fntsize;
bodeopts.XLabel.Interpreter = 'latex';
bodeopts.YLabel.Interpreter = 'latex';
bodeopts.IOGrouping = 'all';

Ghat_plot = Ghat_etfe;
Ghat_plot.InputName = '';
Ghat_plot.OutputName = '';

bode_etfe = bodeplot(Ghat_plot, bodeopts);

ax_bode = findall(f_bode_etfe, 'Type', 'axes');
min_y = inf;
bottom_ax = [];
for i_ax = 1:numel(ax_bode)
    pos = get(ax_bode(i_ax), 'Position');
    if pos(2) < min_y
        min_y = pos(2);
        bottom_ax = ax_bode(i_ax);
    end
end

for i_ax = 1:numel(ax_bode)
    set(ax_bode(i_ax), 'TickLabelInterpreter', 'latex', 'FontSize', fig_setup.fntsize);
    set(ax_bode(i_ax), 'MinorGridLineStyle', 'none');

    ylab = get(ax_bode(i_ax), 'YLabel');
    if ~isempty(get(ylab, 'String'))
        ystr = char(get(ylab, 'String'));
        ystr = strrep(strrep(ystr, '(', '['), ')', ']');
        set(ylab, 'String', ystr, 'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
    end

    xlab = get(ax_bode(i_ax), 'XLabel');
    if ax_bode(i_ax) == bottom_ax
        set(xlab, 'String', 'Frequency [rad/s]', 'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
    else
        set(xlab, 'String', '', 'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
    end
end

line_bode = findall(f_bode_etfe, 'Type', 'line');
for i_ln = 1:numel(line_bode)
    line_bode(i_ln).Color = sscanf(char(color(1)), '#%2x%2x%2x', [1 3]) / 255;
    line_bode(i_ln).LineWidth = 1.25;
end

% h_etfe = line(NaN, NaN, 'Color', color(1), 'LineWidth', 1.25, 'Parent', bottom_ax);
% legend(bottom_ax, h_etfe, {'ETFE'}, 'Interpreter', 'latex', 'Location', 'northwest');

if fig_setup.export_figures
    filename = "output-figures/bode_etfe" + fig_setup.img_ext;
    exportgraphics(f_bode_etfe, filename, ...
        'ContentType', fig_setup.img_format, ...
        'BackgroundColor', 'none');
    fprintf('Figure exported: %s\n', filename);
end

% Result:
% - Anti-resonance around 1.2 rad/s (two complex zeros)
% - Resonance around 2.0 rad/s (two complex poles)
% - Anti-resonance around 2.7 rad/s (two complex zeros)
% - DC gain -50 dB


%%% 2.3 Plot the magnitude of the estimated noise power spectrum
f_error_spectrum = figure('units', 'centimeters', ...
    'Position', [5, 5, fig_setup.fig_wd, fig_setup.fig_hgt*0.5], ...
    'Name', 'Part 2: Estimated noise power spectrum');

specopts = spectrumoptions;
specopts.XLim = [1e-1, 1e1];
specopts.Title.String = '';
specopts.Grid = 'on';
specopts.TickLabel.FontSize = fig_setup.fntsize;
specopts.XLabel.Interpreter = 'latex';
specopts.YLabel.Interpreter = 'latex';

Ghat_spa_plot = Ghat_spa;
Ghat_spa_plot.OutputName = '';

error_spectrum = spectrumplot(Ghat_spa_plot, specopts);

ax_spec = findall(f_error_spectrum, 'Type', 'axes');
for i_ax = 1:numel(ax_spec)
    set(ax_spec(i_ax), 'TickLabelInterpreter', 'latex', 'FontSize', fig_setup.fntsize);
    set(ax_spec(i_ax), 'MinorGridLineStyle', 'none');

    ylab = get(ax_spec(i_ax), 'YLabel');
    if ~isempty(get(ylab, 'String'))
        ystr = char(get(ylab, 'String'));
        ystr = strrep(strrep(ystr, '(', '['), ')', ']');
        set(ylab, 'String', ystr, 'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
    end

    xlab = get(ax_spec(i_ax), 'XLabel');
    if i_ax == 1
        set(xlab, 'String', 'Frequency [rad/s]', 'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
    else
        set(xlab, 'String', '', 'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
    end
end

line_spec = findall(f_error_spectrum, 'Type', 'line');
for i_ln = 1:numel(line_spec)
    xdat = get(line_spec(i_ln), 'XData');
    % Keep vertical reference lines (e.g., Nyquist marker) at default color.
    if isnumeric(xdat) && numel(xdat) >= 2 && (max(xdat) - min(xdat) < 1e-12)
        line_spec(i_ln).LineWidth = 1.25;
        continue;
    end
    line_spec(i_ln).Color = sscanf(char(color(1)), '#%2x%2x%2x', [1 3]) / 255;
    line_spec(i_ln).LineWidth = 1.25;
end

if fig_setup.export_figures
    filename = "output-figures/error_spectrum_spa" + fig_setup.img_ext;
    exportgraphics(f_error_spectrum, filename, ...
        'ContentType', fig_setup.img_format, ...
        'BackgroundColor', 'none');
    fprintf('Figure exported: %s\n', filename);
end

% The noise spectrum is concentrated around 2 rad/s

%% Part 3: Parametric identification and validation
%%% 3.1 Parametric model of G0
N_new = 3000; % New experiment length

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
nb_range = 5;
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

% figure(51); clf;
% alt_orders = optimal_orders;
% alt_orders(5) = 1 - optimal_orders(5);
% try
%     sys_bj_altnk = bj(ze, alt_orders);
%     compare(zv, sys_bj, sys_bj_altnk);
%     legend('Validation Data', ...
%         sprintf('BJ best [%d %d %d %d %d]', optimal_orders), ...
%         sprintf('BJ alt nk [%d %d %d %d %d]', alt_orders));
% catch
%     compare(zv, sys_bj);
%     legend('Validation Data', sprintf('BJ best [%d %d %d %d %d]', optimal_orders));
% end
% grid minor;
% title('Delay comparison for best BJ structure');

% figure(5); clf;
% resid(zv, sys_bj);

e_id = pe(zv, sys_bj);
e_val = e_id.OutputData;
u_val = zv.InputData;
N_res = size(e_val, 1);
max_lag = min(50, N_res - 1);
conf99 = 2.576 / sqrt(N_res);

[Reu, lags_eu] = xcorr(e_val(:,1), u_val(:,1), max_lag, 'coeff');
[Ree, lags_ee] = xcorr(e_val(:,1), e_val(:,1), max_lag, 'coeff');

f_residual = figure('units', 'centimeters', ...
    'Position', [5, 5, fig_setup.fig_wd, fig_setup.fig_hgt*0.6], ...
    'Name', 'Part 3: Residual correlation diagnostics');

tiledlayout(2,1, 'TileSpacing', 'tight');
nexttile;
stem(lags_eu, Reu, 'filled', 'Color', color(1));
hold on;
yline(conf99, '--', 'LineWidth', 1.25, 'Color', color(2));
yline(-conf99, '--', 'LineWidth', 1.25, 'Color', color(2));
yline(0, 'k-', 'LineWidth', 1.25);
grid;
ylabel('$R_{\epsilon u}(\tau)$', 'Interpreter', 'latex');
xticklabels([]);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fig_setup.fntsize);
set(gca, 'MinorGridLineStyle', 'none');
grid on;

nexttile;
stem(lags_ee, Ree, 'filled', 'Color', color(1));
hold on;
yline(conf99, '--', 'LineWidth', 1.25, 'Color', color(2));
yline(-conf99, '--', 'LineWidth', 1.25, 'Color', color(2));
yline(0, 'k-', 'LineWidth', 1.25);
grid minor;
xlabel('Lag $\tau$ [-]', 'Interpreter', 'latex');
ylabel('$R_{\epsilon\epsilon}(\tau)$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fig_setup.fntsize);
set(gca, 'MinorGridLineStyle', 'none');
grid on;

if fig_setup.export_figures
    filename = "output-figures/residual_correlation" + fig_setup.img_ext;
    exportgraphics(f_residual, filename, ...
        'ContentType', fig_setup.img_format, ...
        'BackgroundColor', 'none');
    fprintf('Figure exported: %s\n', filename);
end

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

f_time_bj_validation = figure('units', 'centimeters', ...
    'Position', [5, 5, fig_setup.fig_wd, fig_setup.fig_hgt*0.45], ...
    'Name', 'Part 3: Time-domain validation of BJ model');
[y_bj, fit_bj, ~] = compare(zv, sys_bj);
plot(zv.SamplingInstants, zv.OutputData, 'LineWidth', 0.25, 'Color', color(6), 'DisplayName', 'Data');
hold on;
plot(y_bj.SamplingInstants, y_bj.OutputData, 'LineWidth', 0.25, 'Color', color(1), 'DisplayName', 'BJ');
grid on;
xlabel('Time [samples]', 'Interpreter', 'latex');
ylabel('Output $y(t)$', 'Interpreter', 'latex');
legend('Location', 'northwest', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fig_setup.fntsize);
set(gca, 'MinorGridLineStyle', 'none');
if fig_setup.export_figures
    filename = "output-figures/time_validation_bj" + fig_setup.img_ext;
    exportgraphics(f_time_bj_validation, filename, ...
        'ContentType', fig_setup.img_format, ...
        'BackgroundColor', 'none');
    fprintf('Figure exported: %s\n', filename);
end

% f_ps_prbs = figure('units', 'centimeters', ...
%                'Position', [5, 5, fig_setup.fig_wd, fig_setup.fig_hgt*0.45], ...
%                'Name', 'Part 3: Auto Power Spectrum of PRBS input');
% [Puu, w] = pwelch(u3, [], [], [], 2*pi, 'power');
% plot(w, 10*log10(Puu), 'LineWidth', 1.25, 'Color', color(1));
% grid on;
% xlabel('Frequency [rad/sample]', 'Interpreter', 'latex');
% ylabel('Power Spectrum [dB]', 'Interpreter', 'latex');
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fig_setup.fntsize);
% set(gca, 'MinorGridLineStyle', 'none');
% if fig_setup.export_figures
%     filename = "output-figures/ps_prbs" + fig_setup.img_ext;
%     exportgraphics(f_ps_prbs, filename, ...
%         'ContentType', fig_setup.img_format, ...
%         'BackgroundColor', 'none');
%     fprintf('Figure exported: %s\n', filename);
% end

f_bode_etfe_bj = figure('units', 'centimeters', ...
    'Position', [5, 5, fig_setup.fig_wd, fig_setup.fig_hgt*0.7], ...
    'Name', 'Part 3: Bode plots of etfe and BJ models');

bodeopts = bodeoptions;
bodeopts.XLim = [1e-1, 1e1];
bodeopts.PhaseWrapping = 'on';
bodeopts.Title.String = '';
bodeopts.Grid = 'on';
bodeopts.TickLabel.FontSize = fig_setup.fntsize;
bodeopts.XLabel.Interpreter = 'latex';
bodeopts.YLabel.Interpreter = 'latex';
bodeopts.IOGrouping = 'all';

% Create temporary copies to clear Input/Output names so the built-in labels disappear
Ghat_plot = Ghat_etfe;
Ghat_plot.InputName = '';
Ghat_plot.OutputName = '';
sys_bj_plot = sys_bj;
sys_bj_plot.InputName = '';
sys_bj_plot.OutputName = '';

bode_etfebj = bodeplot(Ghat_plot, '-', sys_bj_plot, '-', bodeopts);
showConfidence(bode_etfebj, 1);
% drawnow;

ax_bode = findall(f_bode_etfe_bj, 'Type', 'axes');
min_y = inf;
bottom_ax = [];
for i_ax = 1:numel(ax_bode)
    pos = get(ax_bode(i_ax), 'Position');
    if pos(2) < min_y
        min_y = pos(2);
        bottom_ax = ax_bode(i_ax);
    end
end

for i_ax = 1:numel(ax_bode)
    set(ax_bode(i_ax), 'TickLabelInterpreter', 'latex', 'FontSize', fig_setup.fntsize);
    set(ax_bode(i_ax), 'MinorGridLineStyle', 'none');

    ylab = get(ax_bode(i_ax), 'YLabel');
    if ~isempty(get(ylab, 'String'))
        ystr = char(get(ylab, 'String'));
        ystr = strrep(strrep(ystr, '(', '['), ')', ']');
        set(ylab, 'String', ystr, 'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
    end

    xlab = get(ax_bode(i_ax), 'XLabel');
    if ax_bode(i_ax) == bottom_ax
        % Setting arbitrary text can trigger bodeplot listeners to append (rad/s)
        % To prevent this effectively, we append a trailing space or intercept it.
        set(xlab, 'String', 'Frequency [rad/s]', 'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
    else
        set(xlab, 'String', '', 'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
    end
end

line_bode = findall(f_bode_etfe_bj, 'Type', 'line');
for i_ln = 1:numel(line_bode)
    ln_col = line_bode(i_ln).Color;
    if isnumeric(ln_col) && numel(ln_col) == 3
        if ln_col(3) >= ln_col(1)
            line_bode(i_ln).Color = sscanf(char(color(1)), '#%2x%2x%2x', [1 3]) / 255;
        else
            line_bode(i_ln).Color = sscanf(char(color(2)), '#%2x%2x%2x', [1 3]) / 255;
        end
    end
    line_bode(i_ln).LineWidth = 1.25;
end

% Create dummy line objects using 'line' to bypass LTI/Bodeplot strict axes checks
h_etfe = line(NaN, NaN, 'Color', color(1), 'LineWidth', 1.25, 'Parent', bottom_ax);
h_bj  = line(NaN, NaN, 'Color', color(2), 'LineWidth', 1.25, 'Parent', bottom_ax);
legend(bottom_ax, [h_etfe, h_bj], {'ETFE', 'BJ'}, 'Interpreter', 'latex', 'Location', 'northwest');

if fig_setup.export_figures
    filename = "output-figures/bode_etfe_bj" + fig_setup.img_ext;
    exportgraphics(f_bode_etfe_bj, filename, ...
        'ContentType', fig_setup.img_format, ...
        'BackgroundColor', 'none');
    fprintf('Figure exported: %s\n', filename);
end

%%% 3.3 Minimum variance estimate
present(sys_bj);

%% Part 4: Experimental verification of variance estimates
%%% 4.1 Monte Carlo simulations
n_mc = 100;
orders = optimal_orders; % [nb nc nd nf nk] selected in 3.2
n_B = orders(1);
n_F = orders(4);

n_total = sum(orders(1:4)); % total free parameters
params_all = zeros(n_mc, n_total);
cov_diags = zeros(n_mc, n_total); % store getcov diagonal from each run
sys_mc_all = cell(n_mc, 1); % store all MC models for envelope/pz visualization

r_mc = idinput(N_new, 'prbs', [0 1], [-M M]);

for i = 1:n_mc
    [u_mc, y_mc] = assignment_sys_18(r_mc, 'open loop');
    data_mc = iddata(y_mc, u_mc, 1, 'Domain', 'Time');
    ze_mc = detrend(data_mc(1:N_new/2), 0);
    sys_mc = bj(ze_mc, orders);
    sys_mc_all{i} = sys_mc;
    params_all(i,:) = getpvec(sys_mc)';
    cov_diags(i,:) = diag(getcov(sys_mc))';
end

% figure(19); clf;
% w_env = logspace(-3, log10(pi), 400);
% mag_mc_db = nan(n_mc, numel(w_env));
% for i = 1:n_mc
%     G_mc = squeeze(freqresp(sys_mc_all{i}, w_env));
%     mag_mc_db(i, :) = 20*log10(abs(G_mc));
% end
% mag_min = min(mag_mc_db, [], 1);
% mag_max = max(mag_mc_db, [], 1);
% mag_med = median(mag_mc_db, 1);
% G_frf = squeeze(freqresp(Ghat_etfe, w_env));
% mag_frf_db = 20*log10(abs(G_frf));

% fill([w_env, fliplr(w_env)], [mag_min, fliplr(mag_max)], ...
%     [0.85 0.90 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
% hold on;
% semilogx(w_env, mag_med, 'b-', 'LineWidth', 1.5);
% semilogx(w_env, mag_frf_db, 'r-', 'LineWidth', 1.5);
% grid minor;
% xlabel('\omega [rad/sample]');
% ylabel('Magnitude [dB]');
% title('Part 4: Bode magnitude envelope of MC BJ estimates');
% legend('MC envelope (min-max)', 'MC median', 'Nonparametric FRF (etfe)', ...
%     'Location', 'best');

% figure(20); clf;
% hold on;
% th = linspace(0, 2*pi, 500);
% plot(cos(th), sin(th), 'k--', 'LineWidth', 1.0);
% for i = 1:n_mc
%     p_i = pole(sys_mc_all{i});
%     z_i = zero(sys_mc_all{i});
%     plot(real(p_i), imag(p_i), 'rx', 'MarkerSize', 4);
%     plot(real(z_i), imag(z_i), 'bo', 'MarkerSize', 4);
% end
% p_best = pole(sys_bj);
% z_best = zero(sys_bj);
% plot(real(p_best), imag(p_best), 'r+', 'MarkerSize', 8, 'LineWidth', 1.5);
% plot(real(z_best), imag(z_best), 'bs', 'MarkerSize', 6, 'LineWidth', 1.2);
% grid minor;
% axis equal;
% xlabel('Real');
% ylabel('Imaginary');
% title('Part 4: Pole-zero map of MC BJ estimates');
% legend('Unit circle', 'MC poles', 'MC zeros', 'Selected BJ poles', 'Selected BJ zeros', ...
%     'Location', 'best');

idx_B = 1:n_B;
idx_F = (n_total - n_F + 1):n_total;
params_B = params_all(:, idx_B);
params_F = params_all(:, idx_F);

% figure(11); clf;
% subplot(2,1,1);
% boxplot(params_B, 'Labels', compose('b_%d', 0:n_B-1));
% ylabel('Value');
% title('B(q) coefficients over 100 MC runs');
% grid minor;
% subplot(2,1,2);
% boxplot(params_F, 'Labels', compose('f_%d', 1:n_F));
% ylabel('Value');
% title('F(q) coefficients over 100 MC runs');
% grid minor;

%%% 4.2 Theoretical variance from one experiment
P_cov = getcov(sys_bj);
var_theo = diag(P_cov);
var_B_theo = var_theo(idx_B);
var_F_theo = var_theo(idx_F);

fprintf('\n--- Theoretical variances (from getcov) ---\n');
fprintf('B coefficients: '); fprintf('%.6e  ', var_B_theo); fprintf('\n');
fprintf('F coefficients: '); fprintf('%.6e  ', var_F_theo); fprintf('\n');

p_ref = getpvec(sys_bj)';
std_theo = sqrt(var_theo)';
params_all_norm = bsxfun(@rdivide, bsxfun(@minus, params_all, p_ref), std_theo);
params_B_norm = params_all_norm(:, idx_B);
params_F_norm = params_all_norm(:, idx_F);

col1_rgb = sscanf(char(color(1)), '#%2x%2x%2x', [1 3]) / 255;
col2_rgb = sscanf(char(color(2)), '#%2x%2x%2x', [1 3]) / 255;

f_MC_violin = figure('units', 'centimeters', ...
    'Position', [5, 5, fig_setup.fig_wd, fig_setup.fig_hgt*0.6], ...
    'Name', 'Part 4: Normalized parameter distribution over MC runs (violin)');

tiledlayout(2,1, 'TileSpacing', 'compact');
nexttile;
hold on;
for j = 1:n_B
    v = params_B_norm(:, j);
    [f, xi] = ksdensity(v);
    f = 0.35 * f / max(f);
    patch([j - f, fliplr(j + f)], [xi, fliplr(xi)], col1_rgb, ...
        'FaceAlpha', 0.35, 'EdgeColor', col1_rgb);
    plot(j + 0*v, v, '.', 'Color', color(1), 'MarkerSize', 5);
end
yline(0, 'k-', 'LineWidth', 1.25);
yline(3, '--', 'Color', col2_rgb, 'LineWidth', 1.25);
yline(-3, '--', 'Color', col2_rgb, 'LineWidth', 1.25);
set(gca, 'XTick', 1:n_B, 'XTickLabel', compose('$b_{%d}$', 0:n_B-1));
xlim([0.5, n_B + 0.5]);
ylabel('z-score [-]', 'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fig_setup.fntsize);
set(gca, 'MinorGridLineStyle', 'none');
grid on;

nexttile;
hold on;
for j = 1:n_F
    v = params_F_norm(:, j);
    [f, xi] = ksdensity(v);
    f = 0.35 * f / max(f);
    patch([j - f, fliplr(j + f)], [xi, fliplr(xi)], col1_rgb, ...
        'FaceAlpha', 0.35, 'EdgeColor', col1_rgb);
    plot(j + 0*v, v, '.', 'Color', color(1), 'MarkerSize', 5);
end
yline(0, 'k-', 'LineWidth', 1.25);
yline(3, '--', 'Color', color(2), 'LineWidth', 1.25);
yline(-3, '--', 'Color', color(2), 'LineWidth', 1.25);
set(gca, 'XTick', 1:n_F, 'XTickLabel', compose('$f_{%d}$', 1:n_F));
xlim([0.5, n_F + 0.5]);
ylabel('z-score [-]', 'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fig_setup.fntsize);
set(gca, 'MinorGridLineStyle', 'none');
grid on;

if fig_setup.export_figures
    filename = "output-figures/MC_violin" + fig_setup.img_ext;
    exportgraphics(f_MC_violin, filename, ...
        'ContentType', fig_setup.img_format, ...
        'BackgroundColor', 'none');
    fprintf('Figure exported: %s\n', filename);
end

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

%% Part 5: MIMO identification
%%% 5.1 Parametric identification of a 2x2 MIMO system G0
N_mimo = 3000;

% Input design: two uncorrelated PRBS signals (one per input channel).
% Independent channels for separating the contribution of each input to each
% output. Uncorrelated inputs prevent ill-conditioning in the regressor
% matrix (analogous to the requirement that Phi_uu is full rank at all
% frequencies for MIMO identifiability).
% Band = [0 1] excites all frequencies up to Nyquist, amplitude within
% saturation bounds [-M, M].
r5 = idinput([N_mimo 2], 'prbs', [0 1], [-M M]); % size (N_mimo, 2)

[u5, y5] = assignment_sys_18(r5, 'MIMO');
nu = size(u5, 2);
ny = size(y5, 2);

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

%%% Kung algorithm (Ho-Kalman adaptation for noisy data)
% Step 0: estimate Markov parameters g_tau from data using impulseest.
ng = 80;      % number of impulse coefficients (must satisfy ng >= no + nc)
no = 20;      % number of block rows in the Hankel matrix
nc = 20;      % number of block columns in the Hankel matrix

if ng < no + nc
    error('Choose ng >= no + nc for Kung realization.');
end

sys_ir = impulseest(ze_mimo, ng);
t_ir = 0:ng;
g_ir = impulse(sys_ir, t_ir); % size: (ng+1) x ny x nu

g_seq = cell(ng + 1, 1);
for tau = 0:ng
    g_seq{tau + 1} = reshape(g_ir(tau + 1, :, :), [ny, nu]);
end

% Steps 1-2: build block Hankel matrices and perform SVD.
H = zeros(no * ny, nc * nu);
H_shift = zeros(no * ny, nc * nu);
for i = 1:no
    row_idx = (i - 1) * ny + (1:ny);
    for j = 1:nc
        col_idx = (j - 1) * nu + (1:nu);
        H(row_idx, col_idx) = g_seq{i + j};      % g_{i+j-1}
        H_shift(row_idx, col_idx) = g_seq{i + j + 1}; % g_{i+j}
    end
end

[U_h, S_h, V_h] = svd(H, 'econ');
sing_vals = diag(S_h);
n_max_plot = 20;

f_svd_mimo = figure('units', 'centimeters', ...
    'Position', [5, 5, fig_setup.fig_wd, fig_setup.fig_hgt*0.5], ...
    'Name', 'Part 5: Singular values of block Hankel matrix'); clf;
stem(1:numel(sing_vals(1:n_max_plot)), sing_vals(1:n_max_plot), 'o-', 'LineWidth', 1.25, 'Color', color(1));
grid on;
xlabel('Order index', 'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
ylabel('Singular value [-]', 'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fig_setup.fntsize);
set(gca, 'MinorGridLineStyle', 'none');

% Step 2.5: choose reduced rank n directly from the SVD drop-off.
% n_chosen = 6 based on being the last significant singular value.
n = 6;

hold on;
xline(n, '--', 'Color', color(2), 'LineWidth', 1.25);
hold off;

if fig_setup.export_figures
    filename = "output-figures/mimo_kung_svd" + fig_setup.img_ext;
    exportgraphics(f_svd_mimo, filename, ...
        'ContentType', fig_setup.img_format, ...
        'BackgroundColor', 'none');
    fprintf('Figure exported: %s\n', filename);
end

sv_ratio = sing_vals(1:end-1) ./ sing_vals(2:end);
drop_db = 20*log10(sv_ratio);

fprintf('\n--- Kung MIMO model-order selection (SVD-based) ---\n');
fprintf('Chosen order n = %d\n', n);
fprintf('Singular-value drop at n: %.2f dB\n', drop_db(n));

% Steps 3-4: realization for chosen reduced rank n.
D_kung = g_seq{1}; % direct feedthrough estimate (g_0)
U_n = U_h(:, 1:n);
S_n = S_h(1:n, 1:n);
V_n = V_h(:, 1:n);

sqrt_Sn = diag(sqrt(diag(S_n)));
inv_sqrt_Sn = diag(1 ./ sqrt(diag(S_n)));

O_n = U_n * sqrt_Sn;
R_n = sqrt_Sn * V_n';

C_n = O_n(1:ny, :);
B_n = R_n(:, 1:nu);
A_n = inv_sqrt_Sn * (U_n' * H_shift * V_n) * inv_sqrt_Sn;

sys_kung_ss = ss(A_n, B_n, C_n, D_kung, 1);
sys_mimo_kung = idss(sys_kung_ss);

% PEM-SS refinement initialized from the Kung realization.
opt_pem_mimo = ssestOptions;
opt_pem_mimo.Focus = 'prediction';
sys_mimo = ssest(ze_mimo, sys_mimo_kung, opt_pem_mimo);

[~, fit_kung] = compare(zv_mimo, sys_mimo_kung);
[~, fit_pem] = compare(zv_mimo, sys_mimo);
fit_kung_avg = mean(fit_kung(:));
fit_pem_avg = mean(fit_pem(:));

fprintf('\n--- MIMO model refinement (Kung -> PEM-SS) ---\n');
fprintf('Validation fit (Kung)   : %.2f%%\n', fit_kung_avg);
fprintf('Validation fit (PEM-SS) : %.2f%%\n', fit_pem_avg);
fprintf('Fit improvement         : %+0.2f%%\n', fit_pem_avg - fit_kung_avg);

% Validation plot: show only the final PEM-SS model (Kung and PEM are
% typically very close, so numerical fit comparison is reported above).
figure(15); clf;
compare(zv_mimo, sys_mimo);
legend('Validation Data', sprintf('PEM-SS (init: Kung, n=%d)', n));
grid minor;
title('MIMO: final PEM-SS model initialized from Kung realization');

% Residual analysis for consistency check (custom plot, SISO-style).
% Using explicit xcorr avoids linked y-axes behavior of resid() and allows
% zooming the cross-correlations to the confidence bounds.
e_mimo_id = pe(zv_mimo, sys_mimo);
e_mimo = e_mimo_id.OutputData;
u_mimo_val = zv_mimo.InputData;
N_res_mimo = size(e_mimo, 1);
max_lag_mimo = min(50, N_res_mimo - 1);
conf99_mimo = 2.576 / sqrt(N_res_mimo);

f_resid_mimo = figure('units', 'centimeters', ...
    'Position', [5, 5, 1.2*fig_setup.fig_wd, fig_setup.fig_hgt*0.7], ...
    'Name', 'Part 5: MIMO residual analysis'); clf;
tiledlayout(ny, nu+1, 'TileSpacing', 'tight', 'Padding', 'compact');

for out_idx = 1:ny
    [Ree_mimo, lags_ee_mimo] = xcorr(e_mimo(:, out_idx), e_mimo(:, out_idx), max_lag_mimo, 'coeff');

    nexttile((out_idx-1) * (nu+1) + 1);
    stem(lags_ee_mimo, Ree_mimo, 'filled', 'Color', color(1));
    hold on;
    yline(conf99_mimo, '--', 'LineWidth', 1.25, 'Color', color(2));
    yline(-conf99_mimo, '--', 'LineWidth', 1.25, 'Color', color(2));
    yline(0, 'k-', 'LineWidth', 1.25);
    grid on;
    xlim([-max_lag_mimo, max_lag_mimo]);
    ylim([-0.1, 1.1]);
    ylabel(sprintf('$R_{\\epsilon_%d\\epsilon_%d}(\\tau)$', out_idx, out_idx), ...
        'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fig_setup.fntsize);
    set(gca, 'MinorGridLineStyle', 'none');
    if out_idx < ny
        xticklabels([]);
    else
        xlabel('Lag $\tau$ [-]', 'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
    end

    for in_idx = 1:nu
        [Reu_mimo, lags_eu_mimo] = xcorr(e_mimo(:, out_idx), u_mimo_val(:, in_idx), max_lag_mimo, 'coeff');

        nexttile((out_idx-1) * (nu+1) + 1 + in_idx);
        stem(lags_eu_mimo, Reu_mimo, 'filled', 'Color', color(1));
        hold on;
        yline(conf99_mimo, '--', 'LineWidth', 1.25, 'Color', color(2));
        yline(-conf99_mimo, '--', 'LineWidth', 1.25, 'Color', color(2));
        yline(0, 'k-', 'LineWidth', 1.25);
        grid on;
        xlim([-max_lag_mimo, max_lag_mimo]);
        ylim([-0.1, 0.1]);
        ylabel(sprintf('$R_{\\epsilon_%d u_%d}(\\tau)$', out_idx, in_idx), ...
            'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fig_setup.fntsize);
        set(gca, 'MinorGridLineStyle', 'none');
        if out_idx < ny
            xticklabels([]);
        else
            xlabel('Lag $\tau$ [-]', 'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
        end
    end
end

if fig_setup.export_figures
    filename = "output-figures/mimo_residual" + fig_setup.img_ext;
    exportgraphics(f_resid_mimo, filename, ...
        'ContentType', fig_setup.img_format, ...
        'BackgroundColor', 'none');
    fprintf('Figure exported: %s\n', filename);
end

% Bode plots of the identified 2x2 transfer matrix, split per input.
w_bode_mimo = logspace(-1, log10(pi), 500);
for in_idx = 1:nu
    G_in = squeeze(freqresp(sys_mimo(:, in_idx), w_bode_mimo)); % ny x Nw
    if isvector(G_in)
        G_in = reshape(G_in, [1, numel(w_bode_mimo)]);
    end

    mag_db = 20*log10(abs(G_in));
    phase_deg = unwrap(angle(G_in), [], 2) * 180/pi;

    f_bode_mimo_in = figure('units', 'centimeters', ...
        'Position', [5, 5, fig_setup.fig_wd, fig_setup.fig_hgt*0.7], ...
        'Name', sprintf('Part 5: Bode plot input u_%d', in_idx)); clf;

    tiledlayout(2,1, 'TileSpacing', 'tight');

    nexttile;
    hold on;
    h_out = gobjects(ny, 1);
    for out_idx = 1:ny
        h_out(out_idx) = semilogx(w_bode_mimo, mag_db(out_idx, :), ...
            'LineWidth', 1.25, 'Color', color(out_idx));
    end
    grid on;
    xlim([1e-1, pi]);
    ylabel('Magnitude [dB]', 'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fig_setup.fntsize);
    set(gca, 'MinorGridLineStyle', 'none');
    set(gca, 'XScale', 'log');
    xticklabels([]);
    legend(h_out, compose('$y_{%d}$', 1:ny), 'Interpreter', 'latex', 'Location', 'best');

    nexttile;
    hold on;
    for out_idx = 1:ny
        semilogx(w_bode_mimo, phase_deg(out_idx, :), ...
            'LineWidth', 1.25, 'Color', color(out_idx));
    end
    grid on;
    xlim([1e-1, pi]);
    xlabel('Frequency [rad/s]', 'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
    ylabel('Phase [deg]', 'Interpreter', 'latex', 'FontSize', fig_setup.fntsize);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fig_setup.fntsize);
    set(gca, 'MinorGridLineStyle', 'none');
    set(gca, 'XScale', 'log');

    if fig_setup.export_figures
        filename = "output-figures/bode_mimo_kung_u" + string(in_idx) + fig_setup.img_ext;
        exportgraphics(f_bode_mimo_in, filename, ...
            'ContentType', fig_setup.img_format, ...
            'BackgroundColor', 'none');
        fprintf('Figure exported: %s\n', filename);
    end
end
