function [ r ] = MRA_LS( R, P )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n=size(R,1);
r=[];
for i=1:n
    R1=getDataVec(R,P,i);
    I=setdiff(1:n,i);
    Y=R1(i,:)';
    X=R1(I,:)';
    ri=pinv(X'*X)*X'*Y;
    ri1=zeros(1,3);ri1(i)=-1;ri1(I)=ri;
    r=[r ;ri1 ];
end

end

