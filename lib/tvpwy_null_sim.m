function sims = tvpwy_null_sim(n, t0, nsim, seed)
% Simulate the null distribution of the TVPWY statistic for a pseudo-sample
% of n intervals: under H0 the pseudo-sample is a random walk with unit-
% variance increments, so the statistic is the sup over tau in [t0,1] of the
% unit-variance DF statistic (test_lr_c) of x(1:i), i = floor(t0*(n+1))+1..n+1.
% Uses running sums (O(n) per path, same quantity as test_lrpwy_c).
% The caller's RNG state is restored on exit.
    old = rng;  cleanup = onCleanup(@() rng(old));
    rng(seed, 'twister');
    N  = n + 1;                       % length of x including x_0 = 0
    r0 = floor(t0 * N);
    sims = zeros(nsim, 1);
    chunk = 1000;
    for a = 1:chunk:nsim
        b = min(a + chunk - 1, nsim);  k = b - a + 1;
        x  = [zeros(1, k); cumsum(randn(n, k), 1)];   % N-by-k
        dx = diff(x, 1, 1);  lx = x(1:end-1, :);       % rows j = 1..n
        Sxy = cumsum(dx .* lx, 1);  Sd = cumsum(dx, 1);
        Sy  = cumsum(lx, 1);        Syy = cumsum(lx.^2, 1);
        m   = (1:n)';                                  % subsample x(1:m+1)
        num = Sxy - Sd .* Sy ./ m;
        den = sqrt(Syy - Sy.^2 ./ m);
        st  = num ./ den;                              % st(m,:) = TVDF of x(1:m+1)
        sims(a:b) = max(st(r0:end, :), [], 1)';        % i = m+1 >= r0+1
    end
end
