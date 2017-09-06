%parpool(4);
pctRunOnAll warning off;
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


MU=[2*ones(1,Nq-5) 0.2 0.2 25 100 400];
A.prior_params=struct('mu',log(MU),'sig',[2*ones(1,Nq-3) 2 2 2],'Indexes',1:Nq);%prior distribution of the parameter




A.N=100;
A.Nb=60;

Yobs=[r1' E1 0 0 0 r2' E2 0 0 0 ];

dr1=1/norm(r1);de1=1/norm([E1]);dr2=1/norm(r2);de2=1/norm([E2]);%sd=dr+de;dr=dr/sd;de=de/sd;
Yobs_params=diag([dr1*ones(6,1);de1*ones(9,1);dr2*ones(6,1);de2*ones(9,1)]);

A.D_obs=Yobs;
A.D_params=Yobs_params;
A.timespan=T;%[0 2 5 8 10 15 30 40 50];
[Q1,V1]=adaptive_weights_abc_smc(A);
%% Model fit
load Q1;
load V1;

    
stimulants=[100 50];
 T=[0 5 10 15 30 60];
 E1=[0 0.76 0.33 0.10 0.03 0.02];E1=E1/max(E1);% ERK
 E2=[0 0.8 0.55 0.30 0.30 0.20];E2=E2/max(E2);

 
r1=[[-1 -0.10 -0.53];[1.11 -1 -0.35];[0.09 1.00 -1]];r1=r1(:);r1=r1([2:4 6:8]);% LRC EGF
r2=[[-1 -0.17 0.40];[6.18 -1 -3.73];[0.96 0.63 -1]];r2=r2(:);r2=r2([2:4 6:8]); % LRC NGF 


pr1=zeros(size(Q1,1),6);
Y1=zeros(size(Q1,1),length(T));
pr2=zeros(size(Q1,1),6);
Y2=zeros(size(Q1,1),length(T));
flag1=zeros(size(Q1,1),1);
di=zeros(size(Q1,1),1);
parfor i=1:size(Q1,1)
    M=SimpleMAPKModel;
    M.timespan=T;%A.timespan;
    M.EGF=stimulants(1);M.NGF=stimulants(2);    
    Oi=M.simulate_observed_metric(M,exp(Q1(i,:)));
    if(length(Oi)>1)
        flag1(i)=1;
        pr1(i,:)=Oi(1:6);
        Y1(i,:)=Oi(7:12);
        pr2(i,:)=Oi(16:21);
        Y2(i,:)=Oi(22:27);
        size(Yobs)
        size(Oi)
        diff=(Yobs-Oi);        
        d(i)=sqrt(diff*Yobs_params*diff');
    end
end
disp('...........')
I1=flag1==1;
Y1_1=Y1(I1,:);
V1_1=V1(I1);
pr1=pr1(I1,:);
[Ym,Ys]=mean_standard_error(Y1_1,V1_1);
[rm,rs]=mean_standard_error(pr1,V1_1);

%I2=flag2==1;
Y1_2=Y2(I1,:);
V1_2=V1(I1);
pr2=pr2(I1,:);
[Ym1,Ys1]=mean_standard_error(Y1_2,V1_2);
[rm1,rs1]=mean_standard_error(pr2,V1_2);

fs=12;
Xlab={'RAF->MEK','RAF->ERK','MEK->RAF','MEK->ERK','ERK->RAF','ERK->MEK'};
figure('Position',[0 0 300 1000]);subplot(4,1,1);bar(r1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',1.5);hold on;errorbar(1:6,rm,rs,'r*','MarkerSize',8,'LineWidth',1.5);xlim([0.5 6.5]);ylim([-0.7 1.5]);set(gca,'XTick',1:6,'XTickLabel',Xlab,'FontSize',fs);rotateXLabels(gca,45);hold off;
subplot(4,1,2);bar(r2,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',1.5);hold on;errorbar(1:6,rm1,rs1,'r*','MarkerSize',8,'LineWidth',1.5);xlim([0.5 6.5]);ylim([-5.2 7]);set(gca,'XTick',1:6,'XTickLabel',Xlab,'FontSize',fs);rotateXLabels(gca,45);hold off;
subplot(4,1,3);plot(T,E1,'ko-','MarkerSize',8,'LineWidth',1.2);hold on;errorbar(T,Ym,Ys,'r*-','MarkerSize',8,'LineWidth',1.2);xlabel('time(minutes)','Fontsize',fs);ylabel('pERK','FontSize',fs);legend('Experimental data','Model fit');set(gca,'FontSize',fs);xlim([0 50]),ylim([0 1]);
subplot(4,1,4);plot(T,E2,'ko-','MarkerSize',8,'LineWidth',1.2);hold on;errorbar(T,Ym1,Ys1,'r*-','MarkerSize',8,'LineWidth',1.2);xlabel('time(minutes)','Fontsize',fs);ylabel('pERK','FontSize',fs);legend('Experimental data','Model fit');set(gca,'FontSize',fs);xlim([0 50]);ylim([0 1.4]);

export_fig(gcf,'model_fit','-jpg','-r300','-q100','-transparent');

%% RAF and MEK phosphorylation
load Q1;
load V1;

stimulants=[100 50];


T=[0 5 10 15 30 60];
%pr1=zeros(size(Q1,1),6);
%Y1=zeros(size(Q1,1),length(A.timespan));
Y1_r=zeros(size(Q1,1),length(T));
Y1_m=zeros(size(Q1,1),length(T));
Y2_r=zeros(size(Q1,1),length(T));
Y2_m=zeros(size(Q1,1),length(T));
flag1=zeros(size(Q1,1),1);
flag2=zeros(size(Q1,1),1);
parfor i=1:size(Q1,1)
    M=SimpleMAPKModel;
    M.timespan=T;
    M.EGF=stimulants(1);M.NGF=stimulants(2);
    
    % Simulate for EGF
    M.params=[exp(Q1(i,:)) 1 0] ;
    y=simulate_model(M);
    if size(y,2)==length(T)
    Y1_r(i,:)=y(1,:)/max(y(1,:));
    Y1_m(i,:)=y(2,:)/max(y(2,:));
    %pr1(i,:)=ri;
    flag1(i)=1;    
    end
    
    M.params=[exp(Q1(i,:)) 0 1];
    y_1=simulate_model(M);
    if size(y_1,2)==length(T)
    %y1=y_1(:,6);
    %
    %ri1=M.r_from_Jacobian_custom(M,y1,EGFt1,NGFt1);ri1=ri1(:);ri1=ri1([2:4 6:8]);    
    %y3=y_1(3,:)/max(y_1(3,:));
    %Y2(i,:)=y3;
    Y2_r(i,:)=y_1(1,:)/max(y_1(1,:));
    Y2_m(i,:)=y_1(2,:)/max(y_1(2,:));
    %pr2(i,:)=ri1;
    flag2(i)=1;
    end
    
end





T=[0 5 10 15 30 60];
E1=[0 0.76 0.33 0.10 0.03 0.02];E1=E1/max(E1);
R1=[0 0.79 0.47 0.34 0.15 0.09];R1=R1/max(R1);
M1=[0 0.76 0.4 0.3 0.14 0.09];M1=M1/max(M1);

E2=[0 0.8 0.55 0.30 0.30 0.20];E2=E2/max(E2);
R2=[0 0.78 0.31 0.39 0.26 0.35];R2=R2/max(R2);
M2=[0 0.75 0.34 0.42 0.35 0.30];M2=M2/max(M2);

I1=flag1==1;
Y1_r1=Y1_r(I1,:);
V1_1=V1(I1);
[Ymr,Ysr]=mean_standard_error(Y1_r1,V1_1);


Y1_m1=Y1_m(I1,:);
[Ymm,Ysm]=mean_standard_error(Y1_m1,V1_1);


I2=flag2==1;
Y2_r1=Y2_r(I2,:);
V2_1=V1(I2);
[Ymr1,Ysr1]=mean_standard_error(Y2_r1,V1_1);

Y2_m1=Y2_m(I2,:);
[Ymm1,Ysm1]=mean_standard_error(Y2_m1,V2_1);

figure('Position',[0 0 700 450]);
subplot(2,2,1);plot(T,R1,'o--','MarkerSize',8,'LineWidth',1.2);hold on;plot(T,M1,'s--','MarkerSize',8,'LineWidth',1.2);hold off; xlabel('time(minutes)','Fontsize',14);ylabel('Rel.Concn.(A.U.)','FontSize',14);legend('pRAF','pMEK');set(gca,'FontSize',14);xlim([0 50]);ylim([0 1]);
subplot(2,2,2);errorbar(T,Ymr,Ysr*sqrt(1000),'o-','MarkerSize',8,'LineWidth',1.2);hold on;errorbar(T,Ymm,Ysm,'s-','MarkerSize',8,'LineWidth',1.2);xlabel('time(minutes)','Fontsize',14);ylabel('Rel.Concn.(A.U.)','FontSize',14);legend('aRAF','aMEK');set(gca,'FontSize',14);xlim([0 50]);ylim([0 1.2]);
subplot(2,2,3);plot(T,R2,'o--','MarkerSize',8,'LineWidth',1.2);hold on;plot(T,M2,'s--','MarkerSize',8,'LineWidth',1.2);xlabel('time(minutes)','Fontsize',14);ylabel('Rel.Concn.(A.U.)','FontSize',14);legend('pRAF','pMEK');set(gca,'FontSize',14);xlim([0 50]),ylim([0 1]);
subplot(2,2,4);hold on;errorbar(T,Ymr1,Ysr1*sqrt(1000),'o-','MarkerSize',8,'LineWidth',1.2);errorbar(T,Ymm1,Ysm1,'s-','MarkerSize',8,'LineWidth',1.2);xlabel('time(minutes)','Fontsize',14);ylabel('Rel.Concn.(A.U.)','FontSize',14);legend('aRAF','aMEK','Location','southwest');set(gca,'FontSize',14);xlim([0 50]);ylim([0 1]);
export_fig(gcf,'model_prediction','-jpg','-r300','-q100','-transparent');



%% Response to ligand blocker
A.timespan=[0 5 10 20 30 40 60];
stimulants=[100 50];
Y1=zeros(size(Q1,1),length(A.timespan));
Y2=zeros(size(Q1,1),length(A.timespan));
parfor i=1:size(Q1,1)
    M=SimpleMAPKModel_blocked;
    M.timespan=A.timespan;
    M.EGF=stimulants(1);M.NGF=stimulants(2);    
    % Simulate for EGF
    M.params=[exp(Q1(i,:)) 1 0] ;
    y=simulate_model(M);
    Y1(i,:)=y(3,:)/max(y(3,:));
    M.params=[exp(Q1(i,:)) 0 1] ;
    y=simulate_model(M); 
    Y2(i,:)=y(3,:)/max(y(3,:));  
end
[Ym1,Ys1]=mean_standard_error(Y1,V1);
[Ym2,Ys2]=mean_standard_error(Y2,V1);

erk1=[0, -0.004255319148936065
4.80836236933798, 0.619403464551363
10.034843205574916, 0.9990510786566833
20.069686411149828, 0.810868114760175
29.89547038327527, 0.49788222502285817
39.93031358885018, 0.37778436750932864
60.000000000000014, 0.3702127659574467]

erk2=[0, 0.06950354609929077
5.017421602787458, 0.9995255393283417
10.034843205574916, 0.1423134900042009
20.069686411149828, -0.00047940297032145196
29.89547038327527, -0.0014085551189857437
39.93031358885018, 0.031685076729186745
60.000000000000014, 0.02411347517730489]

fs=12;
figure('Position',[0 0 600 175]);subplot(1,2,1);hold on;plot(erk2(:,1),erk2(:,2),'o--','Linewidth',1.5);plot(erk1(:,1),erk1(:,2),'s--','Linewidth',1.5);hold off;xlabel('time(minutes)','FontSize',fs);ylabel('pERK(A.U.)','FontSize',fs);legend({'EGF','NGF'},'FontSize',fs);xlim([0 60]);ylim([0 1]);
subplot(1,2,2);errorbar(A.timespan,Ym1,Ys1*sqrt(1000),'o-','MarkerSize',8,'Linewidth',1.5);hold on,errorbar(A.timespan,Ym2,Ys2*sqrt(1000),'s-','MarkerSize',8,'LineWidth',1.5);hold off;xlabel('time(minutes)','FontSize',fs);ylabel('aERK(A.U.)','FontSize',fs);legend({'EGF','NGF'},'FontSize',fs);xlim([0 60]);ylim([0 1]);
export_fig(gcf,'model_fit_blocked','-jpg','-r300','-q100','-transparent');
%% Dose response of pERK
load Q1;
load V1;


egf_doses=[0 0.0001 0.001 0.01 0.03  0.05 0.08 0.1 0.3 0.5 1.0 5 50];%*stimulant/100;
%doses=[0 0.001 0.01 0.03 0.05 0.08 0.1 0.3 0.5 1 5  50]*stimulant/100;
Ym_all=[];
Ys_all=[];

for j=1:length(egf_doses)
 Y=zeros(size(Q1,1),1);
 flag1=zeros(size(Q1,1),1);
 
parfor i=1:size(Q1,1)
    M=SimpleMAPKModel;
    M.EGF=egf_doses(j);
    M.timespan=[0 2 5];
    M.NGF=0;
    M.y0=[0 0 0];
    
    % Simulate for EGF
    M.params=[exp(Q1(i,:)) 1 0];
    y0=simulate_model(M);
  
   if size(y0,2)==length(M.timespan)
       Y(i)=y0(3,3); % ERK level at 5 min       
       flag1(i)=1;
   end
   
  
end
I1=flag1==1;
Y_1=Y(I1);
V_1=V1(I1);

[Ym,Ys]=mean_standard_error(Y_1,V_1);

%[rm,rs]=mean_standard_error(r1,V1);
Ym_all=[Ym_all;Ym];
Ys_all=[Ys_all;Ys];

end

% Simulate NGF levels
ngf_doses=[ 0 0.01 0.1 0.3 1 3 10 50];
%doses=[0 0.001 0.01 0.03 0.05 0.08 0.1 0.3 0.5 1 5  50]*stimulant/100;
Ym_all1=[];
Ys_all1=[];

for j=1:length(ngf_doses)
 Y=zeros(size(Q1,1),1);
 flag1=zeros(size(Q1,1),1);
 
parfor i=1:size(Q1,1)
    M=SimpleMAPKModel;
    M.EGF=0;
    M.timespan=[0 2 5 8 10 15 30 40 50];
    M.NGF=ngf_doses(j);
    M.y0=[0 0 0];
    
    % Simulate for EGF
    M.params=[exp(Q1(i,:)) 0 1];
    y0=simulate_model(M);
  
   if size(y0,2)==length(M.timespan)
       Y(i)=y0(3,3); % ERK level at 5 min       
       flag1(i)=1;
   end
   
  
end
I1=flag1==1;
Y_1=Y(I1);
V_1=V1(I1);

[Ym,Ys]=mean_standard_error(Y_1,V_1);


Ym_all1=[Ym_all1;Ym];
Ys_all1=[Ys_all1;Ys];

end



EGF_data=[0.0001, 0.0042
0.001, 0.001
0.01, 0.016
0.03, 0.065
0.05, 0.14
0.08, 0.21
0.1, 0.27
0.3, 0.56
0.5, 0.7
1, 0.85
5, 1.00
50, 1
];



NGF_data=[
0.01, 0.016
0.1, 0.07
0.3, 0.095
1, 0.37
3, 0.54
10, 1
50, 0.91
];

%D1=(D-min(D))/max(D-min(D));

hill_fit = @(b,x)  b(1).*(x)./(b(2)+x);
b0 = [1; 1];                                  % Initial Parameter Estimates

B1 = lsqcurvefit(hill_fit, b0, egf_doses, Ym_all'/max(Ym_all));x1=linspace(min(egf_doses),max(egf_doses),500000);y1=hill_fit(B1,x1);
B2 = lsqcurvefit(hill_fit, b0, EGF_data(:,1), EGF_data(:,2));x2=linspace(min(EGF_data(:,1)),max(EGF_data(:,1)),500000);y2=hill_fit(B2,x2);
B3 = lsqcurvefit(hill_fit, b0,ngf_doses, Ym_all1'/max(Ym_all1));x3=linspace(min(ngf_doses),max(ngf_doses),500000);y3=hill_fit(B3,x3);
B4 = lsqcurvefit(hill_fit, b0, NGF_data(:,1), NGF_data(:,2));x4=linspace(min(NGF_data(:,1)),max(NGF_data(:,1)),500000);y4=hill_fit(B4,x4);

figure('Position',[0 0 880 185]);subplot(1,4,1);semilogx(egf_doses,Ym_all/max(Ym_all),'or','MarkerSize',5,'LineWidth',2);hold on;
                                                semilogx(egf_doses,(Ym_all+Ys_all*sqrt(1000))/max(Ym_all),'r--','MarkerSize',5,'LineWidth',0.5);
                                                semilogx(egf_doses,(Ym_all-Ys_all*sqrt(1000))/max(Ym_all),'r--','MarkerSize',5,'LineWidth',0.5);
                                                semilogx(x1,y1,'r-','LineWidth',2);hold off;
                                                set(gca,'FontSize',12,'XTick',[0.01 1 100]);xlabel('EGF(ng/ml)','FontSize',12);ylabel('pERK(A.U.)','FontSize',12);xlim([0.0001 100]);ylim([0 1]);
                                 subplot(1,4,2);semilogx(EGF_data(:,1),EGF_data(:,2),'sk','MarkerSize',5,'LineWidth',2);hold on;
                                                semilogx(x2,y2,'k-','LineWidth',2);hold off;
                                                set(gca,'FontSize',12,'XTick',[0.01 1 100]);xlabel('EGF(ng/ml)','FontSize',12);ylabel('pERK(A.U.)','FontSize',12);xlim([0.0001 100]);ylim([0 1]);%semilogx(egf_doses,(Ym_all+Ys_all*sqrt(1000))/max(Ym_all),'r--','LineWidth',0.5);semilogx(egf_doses,(Ym_all-Ys_all*sqrt(1000))/max(Ym_all),'r--','LineWidth',0.5); hold off;ylim([0 1.5]);set(gca,'FontSize',14);xlabel('EGF(ng/ml)','FontSize',14);ylabel('pERK(A.U.)','FontSize',14);legend({'Simulation','Experimental data'},'Location','NorthWest','FontSize',14);%xlim([0.01 100]);%plot(log(doses),D1,'sk','MarkerSize',10,'LineWidth',2);
                                 subplot(1,4,3);semilogx(ngf_doses,Ym_all1/max(Ym_all1),'or','MarkerSize',5,'LineWidth',2);hold on;
                                                semilogx(ngf_doses,(Ym_all1+Ys_all1*sqrt(1000))/max(Ym_all1),'r--','MarkerSize',5,'LineWidth',0.5);
                                                semilogx(ngf_doses,(Ym_all1-Ys_all1*sqrt(1000))/max(Ym_all1),'r--','MarkerSize',5,'LineWidth',0.5);
                                                semilogx(x3,y3,'r-','LineWidth',2);hold off;
                                                set(gca,'FontSize',12,'XTick',[0.01 1 100]);xlabel('EGF(ng/ml)','FontSize',12);ylabel('pERK(A.U.)','FontSize',12);xlim([0.0001 100]);ylim([0 1]);
                                 subplot(1,4,4);semilogx(NGF_data(:,1),NGF_data(:,2),'sk','MarkerSize',5,'LineWidth',2);hold on;
                                                semilogx(x4,y4,'k-','LineWidth',2);hold off;
                                                set(gca,'FontSize',12,'XTick',[0.01 1 100]);ylim([0 1]);xlabel('EGF(ng/ml)','FontSize',12);ylabel('pERK(A.U.)','FontSize',12);xlim([0.0001 100]);ylim([0 1]);
      %semilogx(ngf_doses,(Ym_all1+Ys_all1*sqrt(1000))/max(Ym_all1),'r--','LineWidth',0.5);
      %semilogx(ngf_doses,(Ym_all1-Ys_all1*sqrt(1000))/max(Ym_all1),'r--','LineWidth',0.5);
      %hold off;
      %ylim([0 1.5]);set(gca,'FontSize',14);xlabel('NGF(ng/ml)','FontSize',14);ylabel('pERK(A.U.)','FontSize',14);legend({'Simulation','Experimental data'},'Location','NorthWest','FontSize',14);%xlim([0.01 100]);%plot(log(doses),D1,'sk','MarkerSize',10,'LineWidth',2);
      export_fig(gcf,'model_fit_dose_response','-jpg','-r300','-q100','-transparent');
%figure;hold on;semilogx(doses,Ym_all1/max(Ym_all1),'or','MarkerSize',10,'LineWidth',2);%plot(log(doses),D1,'sk','MarkerSize',10,'LineWidth',2);
%figure;plot(Ym_all')

%%
% Simulate Single cell ERK response to NGF
load Q1;
load V1;
%[Qm,Qse]=mean_standard_
ngf_doses=[0 0.01 0.1 1 3 5 10 30 50 100];
%ngf_doses=[0 0.00001 0.0001 0.001 0.01 0.1 1 10 100 1000 10000];%*stimulant/100;
%doses=[0 0.001 0.01 0.03 0.05 0.08 0.1 0.3 0.5 1 5  50]*stimulant/100;
Ym_all1=[];
Ys_all1=[];
Y_all={};
Y_all1={};
 %timespan=[0 2 5 8 10 15 30 40 50];
timespan=0:120;
 for j=1:length(ngf_doses)
 Y=zeros(size(Q1,1),length(timespan));
 Y1=zeros(size(Q1,1),length(timespan));
 flag1=zeros(size(Q1,1),1);

parfor i=1:size(Q1,1)
    M=SimpleMAPKModel;
    M.EGF=0;
    M.timespan=timespan;
    M.NGF=ngf_doses(j);%ngf_doses(j);
    M.y0=[0 0 0];    
    % Simulate for EGF
    M.params=[exp(Q1(i,:)) 0 1];
    %M.n=2;
    y0=simulate_model(M);  
   if size(y0,2)==length(M.timespan)
       Y(i,:)=y0(3,:); % ERK level 
       Y1(i,:)=y0(2,:);% RAF level
       flag1(i)=1;
   end
   
  
end
I1=(flag1==1)&(sum(Y<0,2)==0);
Y_1=Y(I1,:);
Y1_1=Y1(I1,:);
V_1=V1(I1);

% [fi,xi]=ksdensity(Y_1,'bandwidth',0.05);
% f=fit(xi',fi','gauss2');
Y_all{j}=Y_1;
Y_all1{j}=Y1_1;
[Ym,Ys]=mean_standard_error(Y_1,V_1);
Ym_all1=[Ym_all1;Ym];
Ys_all1=[Ys_all1;Ys];
%Ym_all1=[Ym_all1;[f.b1, f.b2 Ym]];
end

figure('Position',[0 0 700 500]);%hold on;
C=linspecer(4);
for j=2:length(ngf_doses)
    Yss=Y_all{j};Yss1=Y_all1{j};
    Yss_t=Yss(:,60);%trapz(Yss')';
    Yss1_t=Yss1(:,60);%trapz(Yss1')';
    [f,x]=ksdensity(Yss_t);
    [s g]=fit(x.',f.','gauss2');
    [s1 g1]=fit(x.',f.','gauss1');
    f1=s.a1*exp(-((x-s.b1).^2)/((s.c1^2)));
    f2=s.a2*exp(-((x-s.b2).^2)/((s.c2^2)));
    f4=s1.a1*exp(-((x-s1.b1).^2)/((s1.c1^2)));
    f3=f1+f2;
    e1=g.sse;
    e2=g1.sse;

%     lw1=1;
%     ls1='--';
%     lw2=1;
%     ls2='--';
    flag=0;
    if e1>e2
        flag=1;
%     else
%         lw2=2;
%         ls2='-';
    end
    subplot(3,3,j-1);hold on; l1=plot(x,f,'LineWidth',2,'Color',C(1,:));xlim([0 max(x)]);xlabel('pERK level'); ylabel('probability');
    if flag==0
    l2=plot(x,f1,'-','LineWidth',2,'Color',C(2,:));l3=plot(s.b1,s.a1,'o','LineWidth',1,'Color',C(2,:),'MarkerFaceColor',C(2,:),'MarkerSize',6);xlim([0 max(x)]); ylabel('p(pERK)'); 
    l4=plot(x,f2,'-','LineWidth',2,'Color',C(3,:));l5=plot(s.b2,s.a2,'o','LineWidth',1,'Color',C(3,:),'MarkerFaceColor',C(3,:),'MarkerSize',6);xlim([0 max(x)]); ylabel('p(pERK)'); 
    else
        l6=plot(x,f4,'-','LineWidth',2,'Color',C(4,:)); l7=plot(s1.b1,s1.a1,'o','LineWidth',1,'Color',C(4,:),'MarkerFaceColor',C(4,:),'MarkerSize',6);xlim([0 max(x)]); ylabel('probability'); 
    end
    hold off;%N=hist3([Yss1_t,Yss_t],[10,10]);surf(N);%hold on;imagesc(N);contour(N);hold off;
end

export_fig(gcf,'model_fit_bimodality','-jpg','-r300','-q100','-transparent');
