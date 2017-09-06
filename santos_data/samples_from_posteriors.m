function [Q2] = samples_from_posteriors(Q1,V1,N)
%generates samples from posteriors devined by Q1 values and V1 weights, N
%is the number of samples

Qmax=prctile(Q1,90);
I=Q1<repmat(Qmax,size(Q1,1),1);
%size(I)
Qm=zeros(1,size(Q1,2));
Qs=zeros(1,size(Q1,2));

for i=1:size(Q1,2)
    Qi=Q1(I(:,i),i);
    Vi=V1(I(:,i));
    Vi=Vi/sum(Vi);
    [Qmi,Qsi]=mean_standard_error(Qi,Vi);
    Qm(i)=Qmi;
    Qs(i)=Qsi;
end
Q2=repmat(Qm,N,1)+randn(N,size(Q1,2)).*repmat(Qs,size(Q1,1),1);
end

