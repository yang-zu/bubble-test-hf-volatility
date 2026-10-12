function [tr_ret, tv, ret, phi, frac] = tv_truncate(R, varargin)
%TV_TRUNCATE  Jump-robust truncated returns and truncated realized variance.
%
%   [TR_RET, TV] = TV_TRUNCATE(R) takes an M-by-n matrix R of high-frequency
%   log returns, where column s holds the M intraday returns of low-frequency
%   interval s (e.g. the 288 five-minute returns of day s). Each return is
%   kept only if its magnitude is below the local, jump-robust threshold
%
%       phi_s = c * sqrt(MedRV_s) * h^varpi,
%
%   where MedRV_s is the median realized variance of interval s (Andersen,
%   Dobrev and Schaumburg, 2012). It returns, for every interval,
%     TR_RET  the truncated return   sum_j r_{s,j} 1{|r_{s,j}| < phi_s}   (n-by-1)
%     TV      the truncated realized variance (Mancini, 2009)
%             sum_j r_{s,j}^2 1{|r_{s,j}| < phi_s}                          (n-by-1)
%
%   [TR_RET, TV, RET, PHI, FRAC] also returns the untruncated interval
%   returns RET, the thresholds PHI (n-by-1) and the fraction of
%   high-frequency returns that were truncated.
%
%   R may also be a cell array {r_1, ..., r_n} of return vectors of varying
%   length (e.g. the daily returns within each month, or days with missing
%   intraday observations); then h defaults to 1/numel(r_s) for each interval.
%
%   Options (Name, Value):
%     'c'      threshold multiple                  (default 3)
%     'varpi'  threshold exponent                  (default 0.49)
%     'h'      intraday step entering h^varpi      (default 1/M)
%
%   Reference:
%     Boswijk, H. P., Yu, J. and Zu, Y. Testing for explosiveness in financial
%     asset prices using high-frequency volatility. arXiv:2405.02087.

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, 'lib'));

    if iscell(R)
        [tr_ret, tv, ret, phi, frac] = truncate_cell(R, varargin{:});
        return
    end
    if isvector(R)
        error('tv_truncate:shape', ...
            'R must be an M-by-n matrix (one column of M intraday returns per interval) or a cell array.');
    end
    M = size(R, 1);
    ip = inputParser;
    ip.addParameter('c', 3, @(x) isscalar(x) && x > 0);
    ip.addParameter('varpi', 0.49, @(x) isscalar(x) && x > 0 && x < 0.5);
    ip.addParameter('h', 1/M, @(x) isscalar(x) && x > 0);
    ip.parse(varargin{:});
    c = ip.Results.c;  varpi = ip.Results.varpi;  h = ip.Results.h;

    medrv = RV_med(R', 1)';                 % 1-by-n jump-robust local variance
    phi   = c * sqrt(medrv) * h^varpi;      % 1-by-n local thresholds
    keep  = abs(R) < phi;                   % implicit expansion over columns
    Rtr   = R .* keep;

    tr_ret = sum(Rtr, 1)';
    tv     = sum(Rtr.^2, 1)';
    ret    = sum(R, 1)';
    phi    = phi';
    frac   = mean(~keep(:));
end

function [tr_ret, tv, ret, phi, frac] = truncate_cell(R, varargin)
% Variable-length intervals: one return vector per cell.
    ip = inputParser;
    ip.addParameter('c', 3, @(x) isscalar(x) && x > 0);
    ip.addParameter('varpi', 0.49, @(x) isscalar(x) && x > 0 && x < 0.5);
    ip.addParameter('h', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    ip.parse(varargin{:});
    c = ip.Results.c;  varpi = ip.Results.varpi;
    n = numel(R);
    tr_ret = zeros(n,1); tv = zeros(n,1); ret = zeros(n,1); phi = zeros(n,1);
    ntr = 0; ntot = 0;
    for s = 1:n
        r = R{s}(:);
        if numel(r) < 3, error('tv_truncate:short', 'Interval %d has fewer than 3 returns.', s); end
        h = ip.Results.h;  if isempty(h), h = 1/numel(r); end
        phi(s) = c * sqrt(RV_med(r', 1)) * h^varpi;
        keep = abs(r) < phi(s);
        rt = r .* keep;
        tr_ret(s) = sum(rt);
        tv(s)     = sum(rt.^2);
        ret(s)    = sum(r);
        ntr = ntr + sum(~keep);  ntot = ntot + numel(r);
    end
    frac = ntr / ntot;
end
