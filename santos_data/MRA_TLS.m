function [ r ] = MRA_TLS( R, P )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n=size(R,1);
r=[];
for i=1:n
    R1=getDataVec(R,P,i);
    [u v w]= svd(R1');
    r=[r ; (w(:,n)/(-w(i,n)))'];
end

end

