function [Ym,Ys] = mean_standard_error(Y,V)
%Calculates mean and standard error, Y= N X D matrix of data, 
%V=1 X N matrix of probabilities
V=V/sum(V); % normalize just to make sure
Ym=sum(Y.*repmat(V',1,size(Y,2)));
Ys= sqrt(sum(((Y-repmat(Ym,size(Y,1),1)).^2).*repmat(V',1,size(Y,2))))/sqrt(size(Y,1));
end

