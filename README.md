# Testing for explosive bubbles using high-frequency volatility — MATLAB

MATLAB code implementing the TVPWY test and its real-time date-stamping rule
proposed in:

> Boswijk, H. P., Yu, J. and Zu, Y. Testing for explosiveness in financial asset
> prices using high-frequency volatility. arXiv:2405.02087.
> [arxiv.org/abs/2405.02087](https://arxiv.org/abs/2405.02087)

The test detects explosive behaviour (bubbles) in a low-frequency price, e.g. a
daily price, by using the high-frequency (e.g. 5-minute) returns observed
within each day. It works in four steps:

1. **Removes jumps.** High-frequency returns larger than the local threshold
   `phi_s = c * sqrt(MedRV_s) * h^varpi` are discarded, where `MedRV_s` is the
   jump-robust median realized variance of day `s`.
2. **Devolatizes.** Each day's truncated return is divided by the square root of
   its truncated realized variance (TV), giving a pseudo-price that is a
   unit-variance random walk when there is no bubble.
3. **Tests.** It computes the sup of recursive Dickey–Fuller statistics (TVPWY) on
   the pseudo-price. Its null distribution is the PWY distribution of Phillips,
   Wu and Yu (2011), whatever the volatility process and jumps, so time-varying
   volatility and jumps no longer distort size or mask a bubble.
4. **Date-stamps.** The recursive statistics give the origination (and end) of
   an explosive regime in real time.

This is a ready-to-use implementation, not replication code.

## Requirements

Base MATLAB only (R2020b or later). No toolboxes needed.

## Usage

```matlab
% R : M-by-n matrix of high-frequency log returns
%     (column s = the M intraday returns of day s, e.g. M = 288 five-minute returns)
res = tvpwy_test(R);
res = tvpwy_test(R, 'dates', d);      % d: n-by-1 dates, to label the date stamps
```

If you already have each day's truncated return and truncated realized variance:

```matlab
res = tvpwy_from_tv(tr_ret, tv, 'ret', ret, 'dates', d);
[tr_ret, tv, ret] = tv_truncate(R);   % how they are built from R
```

Main fields of `res`:

| field | meaning |
|-------|---------|
| `tvpwy` | the TVPWY statistic (right-tailed test) |
| `pval`, `cv`, `reject5` | simulated p-value; critical values at 10%, 5%, 1%; rejection at 5% |
| `cv_pwy` | PWY (2011) asymptotic critical values used in the paper (for `t0 = 0.1`) |
| `tvdf`, `tvdf_index` | the recursive TVDF detector and the day each value refers to |
| `origination`, `conclusion` | date stamps of the first explosive episode (NaN if none or still open) |
| `episodes` | every episode `[origination conclusion]` found by re-applying the rule |
| `pseudo` | the devolatized pseudo-price |
| `pwy` | classical PWY statistic on the raw price, for comparison |
| `tr_ret`, `tv`, `frac_trunc` | truncated returns, truncated RV, share of returns truncated as jumps |

Options: truncation `'c'` (3), `'varpi'` (0.49), `'h'` (1/M); minimum window
`'t0'` (0.1); date-stamping line `'cvline'` (−0.08, the 5% right-tail critical
value of the Dickey–Fuller distribution); `'nsim'` (5000) and `'seed'` for the
simulated p-value.

## Notes

- **Date stamps.** The origination is the first day the detector exceeds
  `cvline`. The conclusion is the first day, at least `log(n)` days later, that
  it falls back below. Interpret date stamps only when the test rejects, since
  the detector can cross the line by chance when there is no bubble.
- **P-values.** The p-value and `cv` are simulated from the null distribution of
  the statistic at your sample size. The paper reports the tabulated asymptotic
  PWY values in `cv_pwy`; both lead to the same conclusions in its applications.
- **Validation.** This code reproduces the paper's results exactly. The TVPWY
  statistics and detectors for Bitcoin and Ethereum (3.521 and 3.555) match the
  replication files to machine precision, including the Bitcoin origination date
  of 4 April 2019 and the Ethereum episode starting 12 May 2020.

See [`demo.m`](demo.m) for a worked example on simulated 5-minute data with
stochastic volatility and jumps, with and without a bubble.

## How to cite

If you use this code, please cite the paper:

```bibtex
@unpublished{boswijk2026testing,
  title         = {Testing for explosiveness in financial asset prices using high-frequency volatility},
  author        = {Boswijk, H. Peter and Yu, Jun and Zu, Yang},
  year          = {2026},
  note          = {Working paper},
  eprint        = {2405.02087},
  archivePrefix = {arXiv},
  primaryClass  = {econ.EM}
}
```

## License

MIT — see [`LICENSE`](LICENSE).
