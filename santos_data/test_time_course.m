% parpool(4);
pctRunOnAll warning off;
Nq=27;
A=ABC_SMC_WEST;
A.Nq=Nq;
stimulants=[100 50];
A.stimulants=stimulants;
T=[0 5 10 15 30 60];
 E1=[0 0.76 0.33 0.10 0.03 0.02];E1=E1/max(E1);% ERK
 E2=[0 0.8 0.55 0.30 0.30 0.20];E2=E2/max(E2);
 R1=[0 0.79 0.47 0.34 0.15 0.09];R1=R1/max(R1);
M1=[0 0.76 0.4 0.3 0.14 0.09];M1=M1/max(M1);
R2=[0 0.78 0.31 0.39 0.26 0.35];R2=R2/max(R2);
M2=[0 0.75 0.34 0.42 0.35 0.30];M2=M2/max(M2);

% load Q1;load V1;[Qm,Qs]=mean_standard_error(Q1,V1);
% A.prior_params=struct('mu',Qm,'sig',Qs,'Indexes',1:Nq);%prior distribution of the parameter

 MU=[2*ones(1,Nq-3) 20 100 500];
 A.prior_params=struct('mu',log(MU),'sig',[4*ones(1,Nq-3) 1 1 1],'Indexes',1:Nq);%prior distribution of the parameter


A.N=40;
A.Nb=30;

Yobs=[R1 M1 E1 R2 M2 E2];


Yobs_params=eye(36);

A.D_obs=Yobs;
A.D_params=Yobs_params;
A.timespan=T;%[0 2 5 8 10 15 30 40 50];
[Q1,V1]=adaptive_weights_abc_smc(A);
%%
load Q1;
load V1;

    
stimulants=[100 50];
T=[0 5 10 15 30 60];
 E1=[0 0.76 0.33 0.10 0.03 0.02];E1=E1/max(E1);% ERK
 E2=[0 0.8 0.55 0.30 0.30 0.20];E2=E2/max(E2);

% E1=[0 0.7 1 0.405 0.11 0.02 0.02 0.02 0.02]; % time course EGF
% E2=[0 0.2 0.9 1 0.68 0.52 0.58 0.51 0.51];% time course NGF
r1=[[-1 -0.10 -0.53];[1.11 -1 -0.35];[0.09 1.00 -1]];r1=r1(:);r1=r1([2:4 6:8]);% LRC EGF
r2=[[-1 -0.17 0.40];[6.18 -1 -3.73];[0.96 0.63 -1]];r2=r2(:);r2=r2([2:4 6:8]); % LRC NGF 
% 
% %stimulants=[100 50];
% A.timespan=[0 2 5 8 10 15 30 40 50];%[0 2 5 8 10 15 30 40 50];
% E=[0 0.7 1 0.405 0.11 0.02 0.02 0.02 0.02];
% r=[[-1 -0.10 -0.53];[1.11 -1 -0.35];[0.09 1.00 -1]];r=r(:);r=r([2:4 6:8]);
% Yobs=[r1' E1 0 0 0 r2' E2 0 0 0 ];


Y=zeros(size(Q1,1),36);
%Y2_r=zeros(size(Q1,1),length(A.timespan));
%Y2_m=zeros(size(Q1,1),length(A.timespan));
flag1=zeros(size(Q1,1),1);
di=zeros(size(Q1,1),1);
parfor i=1:size(Q1,1)
    M=SimpleMAPKModel;
    M.timespan=T;%A.timespan;
    M.EGF=stimulants(1);M.NGF=stimulants(2);    
    Oi=M.simulate_observed_metric(M,exp(Q1(i,:)));
    if(length(Oi)>1)
        flag1(i)=1;
        Y(i,:)=Oi;
        diff=(Yobs-Oi);        
        d(i)=sqrt(diff*Yobs_params*diff');
    end
end
disp('...........')
I1=flag1==1;
Y1_1=Y(I1,:);
V1_1=V1(I1);

[Ym,Ys]=mean_standard_error(Y1_1,V1_1);

%%
figure;
subplot(2,3,1);plot(T,R1,'s--');hold on; plot(T,Ym(1:6),'*-');hold off;
subplot(2,3,2);plot(T,M1,'s--');hold on; plot(T,Ym(7:12),'*-');hold off;
subplot(2,3,3);plot(T,E1,'s--');hold on; plot(T,Ym(13:18),'*-');hold off;
subplot(2,3,4);plot(T,R2,'s--');hold on; plot(T,Ym(19:24),'*-');hold off;
subplot(2,3,5);plot(T,M2,'s--');hold on; plot(T,Ym(25:30),'*-');hold off;
subplot(2,3,6);plot(T,E2,'s--');hold on; plot(T,Ym(31:36),'*-');hold off;


%% RAF and MEK phosphorylation
load Q1;
load V1;

stimulants=[15.62 1.88];


T=[0 5 10 15 30 60];
%pr1=zeros(size(Q1,1),6);
%Y1=zeros(size(Q1,1),length(A.timespan));
Y1_r=zeros(size(Q1,1),length(T));
Y1_m=zeros(size(Q1,1),length(T));
%pr2=zeros(size(Q1,1),6);
%Y2=zeros(size(Q1,1),length(A.timespan));
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
    %[EGFt,NGFt]=getGrowthFactorLevels(M);
    
%     EGFt1=EGFt(3);
%     NGFt1=NGFt(3);
    %[EGFt1 NGFt1]
    y=simulate_model(M);
    if size(y,2)==length(T)
    %y1=y(:,3);
    %ri=M.r_from_Jacobian_custom(M,y1,EGFt1,NGFt1);ri=ri(:);ri=ri([2:4 6:8]);  
    
    %y3=y(3,:)/max(y(3,:));
    %Y1(i,:)=y3;
    Y1_r(i,:)=y(1,:)/max(y(1,:));
    Y1_m(i,:)=y(2,:)/max(y(2,:));
    %pr1(i,:)=ri;
    flag1(i)=1;    
    end
    
%     EGFt1=EGFt(6);
%     NGFt1=NGFt(6);
    % Simulate for NGF
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

figure('Position',[0 0 500 400]);
subplot(2,2,1);plot(T,R1,'ko-','MarkerSize',8,'LineWidth',1.2);hold on;errorbar(T,Ymr,Ysr,'r*-','MarkerSize',8,'LineWidth',1.2);xlabel('time(minutes)','Fontsize',14);ylabel('pRAF','FontSize',14);legend('Experimental data','Model fit');set(gca,'FontSize',14);xlim([0 50]),ylim([0 1]);
subplot(2,2,2);plot(T,M1,'ko-','MarkerSize',8,'LineWidth',1.2);hold on;errorbar(T,Ymm,Ysm,'r*-','MarkerSize',8,'LineWidth',1.2);xlabel('time(minutes)','Fontsize',14);ylabel('ppMEK','FontSize',14);legend('Experimental data','Model fit');set(gca,'FontSize',14);xlim([0 50])
subplot(2,2,3);plot(T,R2,'ko-','MarkerSize',8,'LineWidth',1.2);hold on;errorbar(T,Ymr1,Ysr1,'r*-','MarkerSize',8,'LineWidth',1.2);xlabel('time(minutes)','Fontsize',14);ylabel('pRAF','FontSize',14);legend('Experimental data','Model fit');set(gca,'FontSize',14);xlim([0 50]),ylim([0 1]);
subplot(2,2,4);plot(T,M2,'ko-','MarkerSize',8,'LineWidth',1.2);hold on;errorbar(T,Ymm1,Ysm1,'r*-','MarkerSize',8,'LineWidth',1.2);xlabel('time(minutes)','Fontsize',14);ylabel('ppMEK','FontSize',14);legend('Experimental data','Model fit');set(gca,'FontSize',14);xlim([0 50])
export_fig(gcf,'model_prediction','-jpg','-r300','-q100','-transparent');



%% blocked
A.timespan=[0 5 10 20 30 40 60];
stimulants=[15.62 1.88];
Y1=zeros(size(Q1,1),length(A.timespan));
Y2=zeros(size(Q1,1),length(A.timespan));
for i=1:size(Q1,1)
    M=SimpleMAPKModel_blocked;
    M.timespan=A.timespan;
    M.EGF=stimulants(1);M.NGF=stimulants(2);    
    % Simulate for EGF
    M.params=[exp(Q1(i,:)) 1 0] ;
    y=simulate_model(M);
    Y1(i,:)=y(3,:);
    M.params=[exp(Q1(i,:)) 0 1] ;
    y=simulate_model(M); 
    Y2(i,:)=y(3,:);  
end
[Ym1,Ys1]=mean_standard_error(Y1,V1);
[Ym2,Ys2]=mean_standard_error(Y2,V1);

figure;plot(A.timespan,Ym1/max(Ym1));hold on,plot(A.timespan,Ym2/max(Ym2));hold off;

%% Dose response of pERK
load Q1;
load V1;


egf_doses=[0 0.0001 0.001 0.01 0.03  0.05 0.08 0.1 0.3 0.5 1.0 5 50]*15.62/100;%*stimulant/100;
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
ngf_doses=[ 0 0.01 0.1 0.3 1 3 10 50]*1.88/100;
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



% EGF_data=[0, 0.0047
% 0.001, 0.0037
% 0.01, 0.0126
% 0.03, 0.062
% 0.05, 0.136
% 0.08, 0.211
% 0.1, 0.27
% 0.3, 0.5535
% 0.5, 0.712
% 1, 0.84647
% 5, 1
% 50, 0.99
% ];

% EGF_data=[0.001, 0
% 0.01, 0.0164
% 0.03, 0.05
% 0.1, 0.06
% 0.4, 0.33
% 1, 0.72
% 10, 1];

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

% NGF_data=[0, 0
% 0.01, 0.0106
% 0.1, 0.0635
% 0.3, 0.09524
% 1, 0.37037
% 3, 0.55
% 10, 1
% 50, 0.915];

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
%%
hill_fit = @(b,x)  b(1).*(x)./(b(2)+x);
b0 = [1; 1];                                  % Initial Parameter Estimates
egf_doses=egf_doses*100/15.62;
ngf_doses=ngf_doses*50/1.88;
B1 = lsqcurvefit(hill_fit, b0, egf_doses, Ym_all'/max(Ym_all));x1=linspace(min(egf_doses),max(egf_doses),500000);y1=hill_fit(B1,x1);
B2 = lsqcurvefit(hill_fit, b0, EGF_data(:,1), EGF_data(:,2));x2=linspace(min(EGF_data(:,1)),max(EGF_data(:,1)),500000);y2=hill_fit(B2,x2);
B3 = lsqcurvefit(hill_fit, b0,ngf_doses, Ym_all1'/max(Ym_all1));x3=linspace(min(ngf_doses),max(ngf_doses),500000);y3=hill_fit(B3,x3);
B4 = lsqcurvefit(hill_fit, b0, NGF_data(:,1), NGF_data(:,2));x4=linspace(min(NGF_data(:,1)),max(NGF_data(:,1)),500000);y4=hill_fit(B4,x4);

figure('Position',[0 0 500 400]);subplot(2,2,1);semilogx(egf_doses,Ym_all/max(Ym_all),'or','MarkerSize',5,'LineWidth',2);hold on;semilogx(x1,y1,'r-','LineWidth',2);hold off;set(gca,'FontSize',12);xlabel('EGF(ng/ml)','FontSize',12);ylabel('pERK(A.U.)','FontSize',12);xlim([0.0001 100]);ylim([0 1]);
                                 subplot(2,2,2);semilogx(EGF_data(:,1),EGF_data(:,2),'sk','MarkerSize',5,'LineWidth',2);hold on;semilogx(x2,y2,'k-','LineWidth',2);hold off;set(gca,'FontSize',12);xlabel('EGF(ng/ml)','FontSize',12);ylabel('pERK(A.U.)','FontSize',12);xlim([0.0001 100]);ylim([0 1]);%semilogx(egf_doses,(Ym_all+Ys_all*sqrt(1000))/max(Ym_all),'r--','LineWidth',0.5);semilogx(egf_doses,(Ym_all-Ys_all*sqrt(1000))/max(Ym_all),'r--','LineWidth',0.5); hold off;ylim([0 1.5]);set(gca,'FontSize',14);xlabel('EGF(ng/ml)','FontSize',14);ylabel('pERK(A.U.)','FontSize',14);legend({'Simulation','Experimental data'},'Location','NorthWest','FontSize',14);%xlim([0.01 100]);%plot(log(doses),D1,'sk','MarkerSize',10,'LineWidth',2);
                                 subplot(2,2,3);semilogx(ngf_doses,Ym_all1/max(Ym_all1),'or','MarkerSize',5,'LineWidth',2);hold on;semilogx(x3,y3,'r-','LineWidth',2);hold off;set(gca,'FontSize',12);xlabel('EGF(ng/ml)','FontSize',12);ylabel('pERK(A.U.)','FontSize',12);xlim([0.0001 100]);ylim([0 1]);
                                 subplot(2,2,4);semilogx(NGF_data(:,1),NGF_data(:,2),'sk','MarkerSize',5,'LineWidth',2);hold on;semilogx(x4,y4,'k-','LineWidth',2);hold off;set(gca,'FontSize',12);ylim([0 1]);xlabel('EGF(ng/ml)','FontSize',12);ylabel('pERK(A.U.)','FontSize',12);xlim([0.0001 100]);ylim([0 1]);
      %semilogx(ngf_doses,(Ym_all1+Ys_all1*sqrt(1000))/max(Ym_all1),'r--','LineWidth',0.5);
      %semilogx(ngf_doses,(Ym_all1-Ys_all1*sqrt(1000))/max(Ym_all1),'r--','LineWidth',0.5);
      %hold off;
      %ylim([0 1.5]);set(gca,'FontSize',14);xlabel('NGF(ng/ml)','FontSize',14);ylabel('pERK(A.U.)','FontSize',14);legend({'Simulation','Experimental data'},'Location','NorthWest','FontSize',14);%xlim([0.01 100]);%plot(log(doses),D1,'sk','MarkerSize',10,'LineWidth',2);
      export_fig(gcf,'model_fit_dose_response','-jpg','-r300','-q100','-transparent');
%figure;hold on;semilogx(doses,Ym_all1/max(Ym_all1),'or','MarkerSize',10,'LineWidth',2);%plot(log(doses),D1,'sk','MarkerSize',10,'LineWidth',2);
%figure;plot(Ym_all')

%%
% Simulate NGF levels
load Q1;
load V1;
%[Qm,Qse]=mean_standard_
ngf_doses=[0 0.01 0.1 1 3 5 10 30 50 100]*1.88/50;
%ngf_doses=[0 0.00001 0.0001 0.001 0.01 0.1 1 10 100 1000 10000];%*stimulant/100;
%doses=[0 0.001 0.01 0.03 0.05 0.08 0.1 0.3 0.5 1 5  50]*stimulant/100;
Ym_all1=[];
Ys_all1=[];
Y_all={};
Y_all1={};
 %timespan=[0 2 5 8 10 15 30 40 50];
timespan=0:60;
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
%%
figure('Position',[0 0 700 500]);%hold on;
C=linspecer(4);
for j=2:length(ngf_doses)
    Yss=Y_all{j};Yss1=Y_all1{j};
    Yss_t=Yss(:,60);%trapz(Yss')';
    Yss1_t=Yss1(:,60);%trapz(Yss1')';
    [f,x]=ksdensity(Yss_t,'support','positive');
    s=fit(x.',f.','gauss2');
    s1=fit(x.',f.','gauss1');
    f1=s.a1*exp(-((x-s.b1).^2)/(2*(s.c1^2)));
    f2=s.a2*exp(-((x-s.b2).^2)/(2*(s.c2^2)));
    f4=s1.a1*exp(-((x-s1.b1).^2)/(2*(s1.c1^2)));
    f3=f1+f2;
    e1=sum((f-f3).^2);
    e2=sum((f-f4).^2);
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
    l2=plot(x,f1,'-','LineWidth',2,'Color',C(2,:));l3=plot(s.b1,s.a1,'o','LineWidth',1,'Color',C(2,:),'MarkerFaceColor',C(2,:),'MarkerSize',6);xlim([0 max(x)]); ylabel('probability'); 
    l4=plot(x,f2,'-','LineWidth',2,'Color',C(3,:));l5=plot(s.b2,s.a2,'o','LineWidth',1,'Color',C(3,:),'MarkerFaceColor',C(3,:),'MarkerSize',6);xlim([0 max(x)]); ylabel('probability'); 
    else
        l6=plot(x,f4,'-','LineWidth',2,'Color',C(4,:)); l7=plot(s1.b1,s1.a1,'o','LineWidth',1,'Color',C(4,:),'MarkerFaceColor',C(4,:),'MarkerSize',6);xlim([0 max(x)]); ylabel('probability'); 
    end
    hold off;%N=hist3([Yss1_t,Yss_t],[10,10]);surf(N);%hold on;imagesc(N);contour(N);hold off;
end
% hl=legend([l1 l4 l2 l3 l4 l5],{'Emperical','Gauss1','Gauss1-peak','Gauss2-dist1','Gauss2-dist1-peak','Gauss2-dist2','Gauss2-dist2-peak'},'Orientation','Horizontal');
% newPosition = [0.4 0.4 0.2 1];
% newUnits = 'normalized';
% set(hl,'Position', newPosition,'Units', newUnits);

export_fig(gcf,'model_fit_bimodality','-jpg','-r300','-q100','-transparent');
%hold off;
%%
%Ym_all2=Ym_all1/max(Ym_all1(:));
%figure('Position',[0 0 450 400]);semilogx(ngf_doses,Ym_all2(:,1),'sr','MarkerSize',10,'LineWidth',2);hold on; semilogx(ngf_doses,Ym_all2(:,2),'ok','MarkerSize',10,'LineWidth',2);semilogx(ngf_doses,Ym_all2(:,3),'db--','MarkerSize',8,'LineWidth',1.5);hold off;legend({'Peak1','Peak2','Average'},'FontSize',14,'Location','NorthWest');
%xlabel('NGF(ng/ml)','FontSize',14);ylabel('pERK(A.U.)','FontSize',14);set(gca,'FontSize',14);
%export_fig(gcf,'model_fit_bimodality','-jpg','-r300','-q100','-transparent');
%%
%figure('Position',[0 0 400 300]);semilogx(ngf_doses,Ym_all1/max(Ym_all1),'or-','MarkerSize',10,'LineWidth',2);