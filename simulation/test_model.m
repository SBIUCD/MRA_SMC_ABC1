% kf4=O.params(1);
%             Kmf4=O.params(2);
%             Vm4=O.params(3);
%             Km4=O.params(4);
%             
%             kf5=O.params(5);
%             Kmf5=O.params(6);
%             Vm5=O.params(7);
%             Km5=O.params(8);
%             
%             kf6=O.params(9);
%             Kmf6=O.params(10);
%             Vm6=O.params(11);
%             Km6=O.params(12);
%             
%             kr4=O.params(13);
%             Kmr4=O.params(14);
%             kr5=O.params(15);
%             Kmr5=O.params(16);
M=SimpleMAPKModel;
M.EGF=2;
M.timespan=0:0.1:60;
M.tot=[100 500 1000];%[100 500 1000];
M.params=exp(log([10 100 20 100 10 100 10 100 1 100 50 100 10 100  10 100]*0.1));% model parameters
M.kf=0.3;
M.sigma=0;
y=simulate_model(M);
figure;plot(M.timespan,y');
r=simulate_local_responses_from_noisy_data(M);
r1=r_from_Jacobian_for_all_EGFs(M);
