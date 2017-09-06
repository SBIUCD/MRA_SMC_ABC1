%parpool(4);
Nq=27;
A=ABC_SMC_WEST;
A.Nq=Nq;
stimulants=[100 50];
A.stimulants=stimulants;
T=[0 5 10 15 30 60];
 E1=[0 0.76 0.33 0.10 0.03 0.02];E1=E1/max(E1);% ERK
 E2=[0 0.8 0.55 0.30 0.30 0.20];E2=E2/max(E2);

r1=[[-1 -0.10 -0.53];[1.11 -1 -0.35];[0.09 1.00 -1]];r1=r1(:);r1=r1([2:4 6:8]);% LRC EGF
r2=[[-1 -0.17 0.40];[6.18 -1 -3.73];[0.96 0.63 -1]];r2=r2(:);r2=r2([2:4 6:8]); % LRC NGF

% E1=[0 0.7 1 0.405 0.11 0.02 0.02 0.02 0.02]; % time course EGF
% E2=[0 0.2 0.9 1 0.68 0.52 0.58 0.51 0.51];% time course NGF

%load Q1_prior;load V1_prior;[Qm,Qs]=mean_standard_error(Q1,V1);
% MU=[2*ones(1,Nq-3) 20 100 500];
 MU=[2*ones(1,Nq-3) 10 50 250];
% % MU([10 16 22])=[50 100 200];
% % MU([2 5 11 17])=[50 50 100 200];
% % MU([4 8 14 20])=[50 50 100 200];
A.prior_params=struct('mu',log(MU),'sig',[2*ones(1,Nq-3) 2 2 2],'Indexes',1:Nq);%prior distribution of the parameter

% MU=Qm;
% A.prior_params=struct('mu',MU,'sig',Qs,'Indexes',1:Nq);%prior distribution of the parameter


A.N=100;
A.Nb=60;

Yobs=[r1' E1 0 0 0 r2' E2 0 0 0 ];

dr1=1/norm(r1);de1=1/norm([E1]);dr2=1/norm(r2);de2=1/norm([E2]);%sd=dr+de;dr=dr/sd;de=de/sd;
Yobs_params=diag([dr1*ones(6,1);de1*ones(9,1);dr2*ones(6,1);de2*ones(9,1)]);

A.D_obs=Yobs;
A.D_params=Yobs_params;
A.timespan=T;%[0 2 5 8 10 15 30 40 50];
load Q1;load V1;
[Q1,V1]=A.adaptive_weights_abc_smc_restart_from_saved_parameters(A,Q1,V1,[1.15:-0.05:1 0.975:-0.025:0.8]);