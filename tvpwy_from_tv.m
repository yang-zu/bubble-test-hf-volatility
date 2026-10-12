function res = tvpwy_from_tv(tr_ret, tv, varargin)
%TVPWY_FROM_TV  TVPWY explosiveness test from truncated returns and truncated RV.
%
%   RES = TVPWY_FROM_TV(TR_RET, TV) runs the TVPWY test of Boswijk, Yu and Zu
%   given, for each of n low-frequency intervals, the truncated return TR_RET
%   and the truncated realized variance TV (both n-by-1; see TV_TRUNCATE). The
%   low-frequency increments are devolatized to form the pseudo-sample
%
%       x_0 = 0,   x_i = sum_{s<=i} TR_RET(s) / sqrt(TV(s)),
%
%   and TVPWY is the supremum over tau in [t0, 1] of the Dickey-Fuller
%   statistic TVDF_tau of x over the first fraction tau of the sample, computed
%   with the error variance fixed at one. Under the null of no explosive
%   behaviour TVPWY has the same asymptotic distribution as the PWY statistic
%   (Phillips, Wu and Yu, 2011), free of the volatility process and jumps. The
%   test is right-tailed.
%
%   Use TVPWY_TEST to start from high-frequency returns instead.
%
%   Options (Name, Value):
%     't0'      minimum window fraction tau_0              (default 0.1)
%     'cvline'  critical value line for date-stamping       (default -0.08, the
%               5% right-tail critical value of the DF t distribution)
%     'nsim'    simulations for p-value and critical values (default 5000)
%     'seed'    RNG seed for the simulation (caller's RNG is restored)
%     'ret'     untruncated interval returns (n-by-1); if given, the classical
%               PWY statistic on the raw log price is also reported
%     'dates'   n-by-1 dates of the intervals, used to label the date stamps
%
%   RES is a struct:
%     tvpwy        the TVPWY statistic
%     pval         simulated p-value (upper tail)
%     cv           simulated critical values at the 10%, 5%, 1% levels
%     reject5      tvpwy > cv(2)
%     tvdf         TVDF detector sequence (the recursive statistics)
%     tvdf_index   interval index s of each detector value (subsample end)
%     origination  first interval with TVDF > cvline            (NaN if none)
%     conclusion   first interval >= origination + log(n) with TVDF < cvline
%                  (NaN if TVDF is still above the line at the end)
%     episodes     K-by-2 [origination conclusion] of every episode found by
%                  re-applying the same rule after each conclusion
%     cv_pwy       PWY (2011) asymptotic critical values used in the paper
%                  (only for t0 = 0.1)
%     pseudo       the pseudo-sample x (n+1-by-1, x(1) = 0)
%     pwy          classical PWY statistic (only if 'ret' is given)
%     *_date       date versions of tvdf_index/origination/conclusion/episodes
%                  (if 'dates' is given)
%
%   Reference:
%     Boswijk, H. P., Yu, J. and Zu, Y. Testing for explosiveness in financial
%     asset prices using high-frequency volatility. arXiv:2405.02087.

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, 'lib'));

    tr_ret = tr_ret(:);  tv = tv(:);
    n = numel(tr_ret);
    if numel(tv) ~= n, error('tvpwy_from_tv:size', 'tr_ret and tv must have the same length.'); end
    if any(tv <= 0), error('tvpwy_from_tv:tv', 'TV must be positive for every interval.'); end

    ip = inputParser;
    ip.addParameter('t0', 0.1, @(x) isscalar(x) && x > 0 && x < 1);
    ip.addParameter('cvline', -0.08, @isscalar);
    ip.addParameter('nsim', 5000, @(x) isscalar(x) && x >= 0);
    ip.addParameter('seed', 20240503, @isscalar);
    ip.addParameter('ret', [], @(x) isempty(x) || numel(x) == n);
    ip.addParameter('dates', [], @(x) isempty(x) || numel(x) == n);
    ip.parse(varargin{:});
    o = ip.Results;

    % pseudo-sample and recursive detector (same computation as the paper's code)
    x = cumsum([0; tr_ret ./ sqrt(tv)]);
    [tvpwy, tvdf] = test_lrpwy_c(x, o.t0);
    r0 = floor(o.t0 * (n + 1));
    s_idx = (r0:n)';                         % detector k uses x(1:r0+k) -> interval r0+k-1

    ep = tvpwy_datestamp(tvdf, s_idx, o.cvline, n);
    if isempty(ep), orig = NaN; conc = NaN; else, orig = ep(1,1); conc = ep(1,2); end

    res = struct('tvpwy', tvpwy, 'pval', NaN, 'cv', NaN(1,3), 'reject5', NaN, ...
                 'tvdf', tvdf, 'tvdf_index', s_idx, ...
                 'origination', orig, 'conclusion', conc, 'episodes', ep, 'pseudo', x, ...
                 't0', o.t0, 'cvline', o.cvline, 'n', n);

    if o.nsim > 0
        sims = tvpwy_null_sim(n, o.t0, round(o.nsim), o.seed);
        res.pval    = mean(sims >= tvpwy);
        ss = sort(sims);                                % base-MATLAB order statistics
        res.cv      = ss(ceil([0.90 0.95 0.99] * numel(ss)))';
        res.reject5 = tvpwy > res.cv(2);
    end
    if abs(o.t0 - 0.1) < 1e-12
        res.cv_pwy = [1.184 1.468 2.094];   % PWY (2011) asymptotic values used in the paper
    end

    if ~isempty(o.ret)
        res.pwy = test_pwy_c([0; cumsum(o.ret(:))], o.t0);
    end

    if ~isempty(o.dates)
        d = o.dates(:);
        res.tvdf_date = d(s_idx);
        res.origination_date = d([]);  res.conclusion_date = d([]);   % empty if no stamp
        if ~isnan(orig), res.origination_date = d(orig); end
        if ~isnan(conc), res.conclusion_date  = d(conc); end
        ed = reshape(d(ones(numel(ep),1)), size(ep));  % K-by-2, same class as dates
        ed(:) = missing;                                % NaT/NaN = still open
        ed(~isnan(ep)) = d(ep(~isnan(ep)));
        res.episodes_date = ed;
    end
end
