%parpool(4);
% Generate simulated
% data....................................................
% sigmas=[0,2,5,10,15,20];
% params=[1 10 2 10 1 10 1 10 0.1 10 10 10 1 10  1 10];
% 
% M=SimpleMAPKModel;
% M.params=params;
% D=[];
% D_params=[];
% NEGF=length(M.EGFs);
% M.kf=0.3;
% for i=1:length(sigmas)
%     M.sigma=sigmas(i);
%     [r,y]=simulate_local_responses_from_noisy_data(M);
%     y1=reshape(y,3,NEGF);y0=repmat(y1(:,1),1,NEGF);y1=y1./y0;y1=y1(:,2:NEGF);
%     D=[D;[r;y1(:)]'];
%     Dw=[1/norm(r) 1/norm(y1(:))];
%     Dw=Dw/sum(Dw);
%     D_params=[D_params; [repmat(Dw(1),length(r),1); repmat(Dw(2),length(y1(:)),1)]'];
%     
% end
% 
% figure;imagesc(D(:,1:18))
% save('Data','D');
% save('Data_params','D_params');
%[r,y]=r_from_Jacobian_for_all_EGFs(M);


%% Fit parameters
load Data;
load Data_params;
A=ABC_SMC_WEST;
A.Nb=600;
A.N=1000;
prior_params=struct('mu',log(10)*ones(1,16),'sig',2*ones(1,16),'Indexes',1:16);%prior distribution of the parameter
A.prior_params=prior_params;


for i=3:length(sigmas)

A.Q1_name=['Q' num2str(i) '.mat'];
A.V1_name=['V' num2str(i) '.mat'];
A.D_obs=D(i,:);
A.D_params=diag(D_params(i,:)');
[Q1,V1]=adaptive_weights_abc_smc(A);
fprintf('************ at noise levels %f ****************',sigmas(i));
end
%MQ1=sum(repmat(V1',1,size(Q1,2)).*exp(Q1)); % mean paramters

%%
%Simulated local response coefficients and steady states


sigmas=[0,2,5,10,15,20];
params=[1 10 2 10 1 10 1 10 0.1 10 10 10 1 10  1 10];
M=SimpleMAPKModel;
M.params=params;


%[rj,yj]=r_from_Jacobian_for_all_EGFs(M);
[Yobs]=simulate_observed_metric(M);


Ym_all=[];
Ys_all=[];
for j=1:(length(sigmas))
Y=[];
load(['Q' num2str(j) '.mat']);
load(['V' num2str(j) '.mat']);
for i=1:size(Q1,1)
    M.params=exp(Q1(i,:));
    [Yi]=simulate_observed_metric(M);
    %Yi=[ri;yi];
    Y=[Y Yi'];
end

%V1=ones(1,size(Q1,1))/size(Q1,1);
Ym=sum(Y.*repmat(V1,size(Y,1),1),2);
Ys=sqrt(sum(((Y-repmat(Ym,1,size(Y,2))).^2).*repmat(V1,size(Y,1),1),2));
Ym_all=[Ym_all;Ym'];
Ys_all=[Ys_all;Ys'];
end

NLRC=(length(M.EGFs)*6);
NSS=length(Yi)-NLRC;
MKs={'o','s','d','x','v','>','<','p','h'};
C=lines(12);
%legends({'\sigma=0','\sigma=2','\sigma=5','\sigma=10','\sigma=15','\sigma=20',});
figure('position',[0 0 1500 1000]);subplot(2,1,1);hold on; bar(Yobs(1:NLRC),'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',1.5);
for i=2:(length(sigmas))%[2:4 6]%length(sigmas)
    Ym=Ym_all(i,:);
    Ys=Ys_all(i,:);
    errorbar(Ym(1:NLRC),Ys(1:NLRC),'.','Marker',MKs{i},'color',C(i,:),'LineWidth',1.5,'MarkerSize',12);set(gca,'FontSize',14);xlim([0.5 NLRC+0.5]);
end
    Ym=Ym_all(1,:);
    Ys=Ys_all(1,:);
    errorbar(Ym(1:NLRC),Ys(1:NLRC),'.','Marker','o','color','r','LineWidth',1.5,'MarkerSize',12);set(gca,'FontSize',14);xlim([0.5 NLRC+0.5]);
hold off;
legend({'Input data','\sigma=2','\sigma=5','\sigma=10','\sigma=15','\sigma=20','\sigma=0'},'Orientation','horizontal','Location','NorthOutside');
subplot(2,1,2),hold on; bar(Yobs((1+NLRC):end),'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',1.5);
for i=2:(length(sigmas))
    Ym=Ym_all(i,:);
    Ys=Ys_all(i,:);
    errorbar(Ym((1+NLRC):end),Ys((1+NLRC):end),'.','Marker',MKs{i},'color',C(i,:),'LineWidth',1.5,'MarkerSize',15);set(gca,'FontSize',14);xlim([0.5 NSS+0.5]);
end
Ym=Ym_all(1,:);
Ys=Ys_all(1,:);
errorbar(Ym((1+NLRC):end),Ys((1+NLRC):end),'.','Marker','o','color','r','LineWidth',1.5,'MarkerSize',15);set(gca,'FontSize',14);xlim([0.5 NSS+0.5]);
hold off;
legend({'Real','\sigma=2','\sigma=5','\sigma=10','\sigma=15','\sigma=20','\sigma=0'},'Orientation','horizontal','Location','NorthOutside');
export_fig(gcf,'model_fit','-jpg','-r300','-q100','-transparent');
%%
%Simulate dose response

sigmas=[0,2,5,10,15,20];
params=[1 10 2 10 1 10 1 10 0.1 10 10 10 1 10  1 10];
EGFs=[0.2 0.4 0.8 1.6 3.2 6.4 12.8 25.6 51.2 102.4];
M=SimpleMAPKModel;
M.params=params;
M.EGFs=EGFs;
N=1000;

%[rj,yj]=r_from_Jacobian_for_all_EGFs(M);
[Yobs]=simulate_model_steadystate_for_all_EGFs(M);


Ym_all=[];
Ys_all=[];
for j=1:(length(sigmas))
Y=[];
load(['Q' num2str(j) '.mat']);
load(['V' num2str(j) '.mat']);
for i=1:size(Q1,1)
    M.params=exp(Q1(i,:));
    [Yi]=simulate_model_steadystate_for_all_EGFs(M);%simulate_observed_metric(M);
    %Yi=[ri;yi];
    Y(:,:,i)=Yi;
end

%V1=ones(1,size(Q1,1))/size(Q1,1);
N=size(Q1,1);
W=reshape(V1,1,1,N);
W=repmat(W,size(Y,1),size(Y,2),1);
Ym=sum(Y.*W,3);
size(Ym)
Ys=sqrt(sum(((Y-repmat(Ym,1,1,size(Y,3))).^2).*W,3));
Ym_all(:,:,j)=Ym;
Ys_all(:,:,j)=Ys.*(1/sqrt(N));
end

MKs={'o','s','d','x','v','>','<','p','h'};
C=lines(12);
figure('Position',[0,0,750,500]);h=semilogx(EGFs,Yobs(3,:)/max(Yobs(3,:)),'-ko');set(h,'MarkerFaceColor',[0.7 0.7 0.7], 'MarkerSize',8, 'LineWidth',1.5);
hold on; 

for i=1:size(Ym_all,3)
h1= semilogx(EGFs,Ym_all(3,:,i)/max(Ym_all(3,:,i)));set(h1,'Color',C(i,:),'Marker',MKs{i},'MarkerEdgeColor',C(i,:)/1.5,'MarkerFaceColor',C(i,:), 'MarkerSize',8, 'LineWidth',1.5)
end

xlim([min(EGFs) max(EGFs)]);
xlabel('EGF(ng/ml)','FontSize',16);ylabel('Normalized aERK (A.U.)','FontSize',16);
legend({'TrueValues','\sigma=0','\sigma=2','\sigma=5','\sigma=10','\sigma=15','\sigma=20'},'Orientation','Horizontal','Location','NorthOutside')
set(gca,'FontSize',20);
export_fig(gcf,'dose_response','-jpg','-r300','-q100','-transparent');



%%

%Effect of prior params

   A=ABC_SMC_WEST; 
   params=[1 10 2 10 1 10 1 10 0.1 10 10 10 1 10  1 10];
   %A.sigma=5;
   priors=[1 5 10 20 30 50 100 1000];
   %priors=[20 30 50 100 1000];
   
   load Data;
   load Data_params;
   for j=1:size(D,1)
   Ym_all=zeros(1,length(priors));
   Ys_all=zeros(1,length(priors));   
   for i=[1:length(priors)]%3:length(sigmas)
       %A.sigma=sigmas(i);
       A.D_obs=D(j,:);
       A.D_params=diag(D_params(j,:)');
       prior_params=struct('mu',log(priors(i))*ones(1,16),'sig',3*ones(1,16),'Indexes',1:16);%prior distribution of the parameter
       A.Q1_name=['Q_p' num2str(j) '_' num2str(i) '.mat'];
       A.V1_name=['V_p' num2str(j) '_' num2str(i) '.mat'];
       [Q1,V1]=adaptive_weights_abc_smc(A);
       load(A.Q1_name);
       load(A.V1_name);
       M=SimpleMAPKModel;
       M.params=params;
       [Ym,Ys]=M.posterior_error(M,exp(Q1),V1,D(3,:));
       Ym_all(i)=Ym;
       Ys_all(i)=Ys;
       %err(i)=sum((Y_obs-Ym).^2,2);
       fprintf('************ at prior level %f ****************',priors(i));
   end
  
   figure('Position',[0 0 600 300]);errorbar(Ym_all,Ys_all,'b-','LineWidth',2);xlim([0.5 length(priors)]);set(gca,'XTick',1:length(priors),'XTickLabel',priors,'FontSize',18);ylabel('Sum of squared error','FontSize',14), xlabel('Prior Mean','FontSize',14)
    export_fig(gcf,['sensitivity_to_prior_mean' num2str(j)],'-jpg','-r300','-q100','-transparent');
   end