%% demo.m  -  TVPWY explosiveness test with high-frequency volatility
%
% Simulates intraday (5-minute) log prices with stochastic volatility and
% jumps, then tests for an explosive bubble in the daily price. Case 1 has no
% bubble; case 2 has an explosive regime starting on day 350.

clear;
n = 500;            % trading days
M = 78;             % 5-minute returns per 6.5-hour trading day
start = 350;        % bubble origination day (case 2)

R0 = simulate_hf(n, M, 0,    start, 1);   % no bubble
R1 = simulate_hf(n, M, 0.01, start, 2);   % explosive drift kappa = 0.01 per day

for k = 1:2
    if k == 1, R = R0; label = 'No bubble'; else, R = R1; label = 'Bubble from day 350'; end
    res = tvpwy_test(R);
    fprintf('\n%s\n', label);
    fprintf('  TVPWY = %6.3f   p-value = %.3f   5%% cv = %.3f   (PWY on raw prices = %.3f)\n', ...
            res.tvpwy, res.pval, res.cv(2), res.pwy);
    fprintf('  HF returns truncated as jumps: %.2f%%\n', 100*res.frac_trunc);
    if res.reject5          % date stamps are interpreted only after the test rejects
        fprintf('  explosive regime date-stamped from day %d\n', res.origination);
    end
end

% detector for the bubble case
figure;
plot(res.tvdf_index, res.tvdf); hold on;
yline(res.cvline, 'r--', 'cv line');  xline(start, 'k:', 'true start');
xlabel('day'); ylabel('TVDF_\tau'); title('Real-time TVDF detector'); hold off;

% If you already have daily truncated returns and truncated RV (e.g. built
% from 5-minute data elsewhere), use the second entry point directly:
%   res = tvpwy_from_tv(tr_ret, tv, 'ret', ret, 'dates', dates);

%% ---------------------------------------------------------------------------
function R = simulate_hf(n, M, kappa, start, seed)
% M-by-n matrix of 5-minute log returns. Daily variance follows a persistent
% log-AR(1) process, jumps arrive at rate 0.1 per day, and from day START the
% log price has explosive drift kappa * p (per day).
    rng(seed);
    dt   = 1/M;
    lv   = zeros(n,1);                         % log daily variance
    for t = 2:n, lv(t) = 0.98*lv(t-1) + 0.2*randn; end
    sig2 = 0.02^2 * exp(lv - mean(lv));        % daily variance, ~2% daily vol
    p = 0;  R = zeros(M, n);
    for t = 1:n
        for j = 1:M
            jump = (rand < 0.1*dt) * 0.03*randn;
            dp = (t >= start) * kappa * p * dt + sqrt(sig2(t)*dt)*randn + jump;
            R(j,t) = dp;  p = p + dp;
        end
    end
end
