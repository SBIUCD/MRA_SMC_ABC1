function [ X ] = getDataVec(data,pert,k)
%provides the data from the perturbation experiments where the node in
%concern is not perturbed
s1=size(pert);
X=[];
for i=1:s1(2)
    if pert(k,i)==0 
        X=[X data(:,i)];
    end
end

end

