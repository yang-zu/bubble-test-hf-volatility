function dft = test_lr_c(pri)
% LR statistic in OLS with a constant
% pri can be a matrix, where each column is a sample
dy = diff(pri,1,1);
ly = pri(1:end-1,:);
n = size(ly,1);
ly_dm = ly - repmat(mean(ly,1),n,1);
dft = sum(dy.*ly_dm,1)./sqrt(sum(ly_dm.^2,1));
end