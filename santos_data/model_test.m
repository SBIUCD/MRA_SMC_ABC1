load results/Q1;
load results/V1;
[Qm,Vm]=mean_standard_error(exp(Q1),V1);
M=SimpleMAPKModel;

P=[Qm];
M.params=[P 1 0];
M.EGF=15.62;
M.NGF=1.88;
M.timespan=1:60;
M.y0=zeros(3,1);
y=simulate_model(M);
figure;plot(y');

M.params=[P 0 1];
y=simulate_model(M);
figure;plot(y');