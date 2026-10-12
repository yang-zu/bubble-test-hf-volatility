function [supdft,dft,i_pwy] = test_lrpwy_c( y,t0 )
% sup DF t type statistic
n=length(y);    
r0=floor(t0*n);
dft=zeros(n,1);
for i=r0+1:n
    dft(i)=test_lr_c(y(1:i));
    %[~,~,dft(i-r0,:),~,~]=adftest(y(1:i,:),'model','ARD');
end
[supdft,i_pwy]=max(dft(r0+1:end));
dft = dft(r0+1:end);
end