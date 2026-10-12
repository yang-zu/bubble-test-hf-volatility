function dft = test_df_c(pri)
% t version DF statistic in OLS with a constant
% pri can be a matrix, where each column is a sample
dy = diff(pri,1);
ly = pri(1:end-1);
n = length(ly);
ly_dm = ly - mean(ly);
khat = sum(dy.*ly_dm)./sum(ly_dm.^2);
alpha = mean(dy) - khat*mean(ly);
u = dy - khat*ly - alpha;
dft = sum(dy.*ly_dm)./sqrt(sum(ly_dm.^2))./sqrt(sum(u.^2)/(n-2));


end