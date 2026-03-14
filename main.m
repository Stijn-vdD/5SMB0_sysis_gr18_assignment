clearvars;

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
n = 6;      % Number of periods
if mod(N,n) ~= 0
    error("n is not a clean divisor of N, reference signal cannot be periodic.")
end
p = N/n;    % Period length

r2 = repmat(M * multisine([0,0.49],1,p),1,N);

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
G0_etfe = etfe(frfdata);

f2 = figure(2);clf;
G0_bode = bodeplot(G0_etfe);
G0_bode.PhaseWrappingEnabled = true;
legend;