function [Q1,I1] = remove_outliers(Q)
%Identify and remove outliers
I1=true(size(Q,1),1);
for i=1:size(Q,2)
    l=prctile(Q(:,i),1);
    h=prctile(Q(:,i),99);
    I=(Q(:,i)>l)&(Q(:,i)<h);
    I1=I1&I;
end

Q1=Q(I1,:);
end

