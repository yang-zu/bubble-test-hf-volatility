function res = tvpwy_test(R, varargin)
%TVPWY_TEST  Test for an explosive bubble using high-frequency volatility (TVPWY).
%
%   RES = TVPWY_TEST(R) tests for explosive behaviour in a low-frequency
%   (e.g. daily) log price, using the high-frequency (e.g. 5-minute) returns
%   observed within each low-frequency interval, following Boswijk, Yu and Zu.
%   R is an M-by-n matrix: column s holds the M intraday log returns of
%   interval s. R may also be a cell array of return vectors of varying length
%   (e.g. days with missing intraday data, or daily returns within months).
%   The test
%
%     1. removes jumps by truncating each high-frequency return at the local
%        threshold phi_s = c * sqrt(MedRV_s) * h^varpi;
%     2. devolatizes each interval's truncated return by the truncated realized
%        volatility sqrt(TV_s), giving a pseudo-sample that is a unit-variance
%        random walk under the null;
%     3. computes the sup of recursive Dickey-Fuller statistics (TVPWY) on the
%        pseudo-sample, whose null distribution is that of the PWY statistic and
%        does not depend on the volatility process or jumps; and
%     4. date-stamps the origination (and conclusion) of the explosive regime
%        in real time from the recursive statistics.
%
%   Options (Name, Value):
%     'c', 'varpi', 'h'   truncation settings (defaults 3, 0.49, 1/M)
%     't0'                minimum window fraction (default 0.1)
%     'cvline'            date-stamping critical value (default -0.08)
%     'nsim', 'seed'      simulation for p-value / critical values (default 5000)
%     'dates'             n-by-1 dates of the intervals
%
%   RES contains all fields of TVPWY_FROM_TV (tvpwy, pval, cv, reject5, tvdf,
%   origination, conclusion, pseudo, pwy, ...) plus
%     tr_ret, tv   truncated interval returns and truncated realized variance
%     ret          untruncated interval returns
%     frac_trunc   fraction of high-frequency returns truncated
%   RES.pwy is the classical PWY statistic on the raw log price, for comparison.
%
%   Example:
%     R   = randn(288, 500) * 0.01/sqrt(288);  % 500 days of 5-minute returns
%     res = tvpwy_test(R);
%
%   Reference:
%     Boswijk, H. P., Yu, J. and Zu, Y. Testing for explosiveness in financial
%     asset prices using high-frequency volatility. arXiv:2405.02087.

    ip = inputParser;
    ip.KeepUnmatched = true;
    ip.addParameter('c', 3);
    ip.addParameter('varpi', 0.49);
    ip.addParameter('h', []);
    ip.parse(varargin{:});
    trunc = {'c', ip.Results.c, 'varpi', ip.Results.varpi};
    if ~isempty(ip.Results.h), trunc = [trunc, {'h', ip.Results.h}]; end
    rest = namedargs2cell(ip.Unmatched);

    [tr_ret, tv, ret, ~, frac] = tv_truncate(R, trunc{:});
    res = tvpwy_from_tv(tr_ret, tv, 'ret', ret, rest{:});
    res.tr_ret = tr_ret;
    res.tv = tv;
    res.ret = ret;
    res.frac_trunc = frac;
end
