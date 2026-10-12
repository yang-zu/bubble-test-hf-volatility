function RVmed=RV_med(rtn,T)
% Median realized variance (MedRV) jump-robust IV estimator (Andersen,
% Dobrev & Schaumburg 2012; Mancini-style truncation building block).
% Ported from the PS-paper code (signhf_acf/Simulations1) for the TVPWY
% truncation-threshold robustness check.
%
% input:
%   rtn : n1 x M matrix of high-frequency returns (each row = one
%         low-frequency interval, M intraday returns)
%   T   : if 1, return per-interval MedRV (n1 x 1); else return the sum.
% output:
%   RVmed : n1 x 1 (T==1) MedRV per interval, an estimate of the interval
%           integrated variance, robust to jumps.
lag=1;
r10=rtn(:,1:end-lag-1);
r20=rtn(:,lag+1:end-lag);
r30=rtn(:,lag+2:end);

const=pi/(6-4*sqrt(3)+pi);
M = size(rtn,2);                 % intraday returns per interval
ssc = M/(M-2);                   % ADS(2012) finite-sample scaling (restored)

n1=size(rtn,1);
RV_med0=nan(n1,1);
for t=1:n1
    tmp_med = median([abs(r10(t,:));abs(r20(t,:));abs(r30(t,:))]);
    RV_med0(t)=const*ssc*sum(tmp_med.^2,'omitnan');
end


if T==1
    RVmed=RV_med0;
else
    RVmed=sum(RV_med0);
end

end
