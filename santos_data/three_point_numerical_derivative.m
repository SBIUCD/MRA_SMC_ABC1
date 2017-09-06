function [d] = three_point_numerical_derivative(Y,X)
%Three point numerical derivative estimator ay Xi
Yi_1=Y(1);
Yi=Y(2);
Yip1=Y(3);

Xi_1=X(1);
Xi=X(2);
Xip1=X(3);

del1=(Xi-Xip1)/((Xi_1-Xi)*(Xi_1-Xip1));
del2=(2*Xi-Xi_1-Xip1)/((Xi-Xi_1)*(Xi-Xip1));
del3=(Xi-Xi_1)/((Xip1-Xi_1)*(Xip1-Xi));
d=Yi_1*del1+Yi*del2+Yip1*del3;
end

