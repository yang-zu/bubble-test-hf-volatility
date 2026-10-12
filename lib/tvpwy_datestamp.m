function ep = tvpwy_datestamp(tvdf, s_idx, cv, n)
% Real-time date stamps (Boswijk, Yu and Zu, Section 5), applied repeatedly:
%   origination: first interval with TVDF > cv;
%   conclusion : first interval at least log(n) intervals after the
%                origination with TVDF < cv (NaN if still above at the end).
% After a conclusion the search restarts, so every episode is returned.
% tvdf : detector sequence; s_idx : interval index of each detector value.
% ep   : K-by-2 matrix [origination, conclusion] of interval indices.
    tvdf = tvdf(:);  s_idx = s_idx(:);
    ep = zeros(0, 2);
    k0 = 1;
    while k0 <= numel(tvdf)
        ke = find(tvdf(k0:end) > cv, 1, 'first');
        if isempty(ke), break; end
        ke = ke + k0 - 1;
        orig = s_idx(ke);
        kf = find(s_idx >= orig + log(n) & tvdf < cv & (1:numel(tvdf))' > ke, 1, 'first');
        if isempty(kf)
            ep(end+1, :) = [orig, NaN]; %#ok<AGROW>
            break;
        end
        ep(end+1, :) = [orig, s_idx(kf)]; %#ok<AGROW>
        k0 = kf + 1;
    end
end
