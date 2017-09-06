% % Simulate NGF levels
% load Q1;
% load V1;
% %doses=[0 0.001 0.01 0.03 0.05 0.08 0.1 0.3 0.5 1 5  50]*stimulant/100;
% timespan=[0:60];
% 
% 
% Y1=zeros(size(Q1,1),length(timespan));
% Y2=zeros(size(Q1,1),length(timespan));
% flag1=zeros(size(Q1,1),1);
% parfor i=1:size(Q1,1)
%     M=SimpleMAPKModel;
%     M.EGF=100;
%     %M.n=1;
%     M.timespan=timespan;
%     M.NGF=0;%ngf_doses(j);
%     M.y0=[0 0 0];    
%     % Simulate for EGF
%     M.params=[exp(Q1(i,:)) 1 0];
%     y0=simulate_model(M);  
%    if size(y0,2)==length(M.timespan)
%        Y1(i,:)=y0(1,:)/max(y0(1,:)); % ERK level at 5 min 
%        Y2(i,:)=y0(2,:)/max(y0(2,:));
%        flag1(i)=1;
%    end 
%   
% end
% 
% 
% Y31=zeros(size(Q1,1),length(timespan));
% Y21=zeros(size(Q1,1),length(timespan));
% flag2=zeros(size(Q1,1),1);
% flag3=zeros(size(Q1,1),1);
% parfor i=1:size(Q1,1)
%     M=SimpleMAPKModel;
%     %M.n=1;
%     M.EGF=0;
%     M.timespan=timespan;
%     M.NGF=200;%ngf_doses(j);
%     M.y0=[0 0 0];    
%     % Simulate for EGF
%     M.params=[exp(Q1(i,:)) 0 1];
%     y0=simulate_model(M);  
%    if size(y0,2)==length(M.timespan)           
%        Y21(i,:)=y0(2,:)/max(y0(2,:));
%        flag2(i)=1;
%    end 
%    
%    M.NGF=50;
%     y0=simulate_model(M);  
%    if size(y0,2)==length(M.timespan)           
%        Y31(i,:)=y0(2,:)/max(y0(2,:));
%        flag3(i)=1;
%    end 
%   
% end
% 
% MEK_EGF=[0
% 1
% 0.378490899
% 0.321084852
% 0.186799571
% 0.134931336
% ];
% T1=[0 2 5 30 60 120];
% 
% MEK_NGF=[-9.7239E-08
% 0.999999999
% 1.029222406
% 0.864747207
% 0.557374398
% 0.691585405
% 0.88203768
% 0.482617818
% 0.431627178
% 0.121503999
% ];
% T2=[0
%     2
%     10
%     20
%     40
%     60
%     90
%     180
%     300
%     420];
% 
% 
% % MEK_NGF1=[0,0
% %     6, 1.0165289256198347
% % 15, 0.47107438016528924
% % 30, 0.28925619834710736
% % 60, 0.3223140495867769];
% 
% 
% 
% I1=flag1==1;
% Y2_1=Y2(I1,:);
% V1_1=V1(I1);
% [Ym1,Ys1]=mean_standard_error(Y2_1,V1_1);
% 
% I2=flag2==1;
% Y2_2=Y21(I2,:);
% V1_2=V1(I2);
% [Ym2,Ys2]=mean_standard_error(Y2_2,V1_2);
% 
% I3=flag3==1;
% Y3_2=Y31(I3,:);
% V1_3=V1(I3);
% [Ym3,Ys3]=mean_standard_error(Y3_2,V1_3);
% 
% 
% 
% figure;hold on;
% C=parula(10);
% plot(timespan,Ym1,'b-','Linewidth',2);
% %plot(timespan,Ym1+Ys1*sqrt(sum(I1)),'b--','Linewidth',1);
% %plot(timespan,Ym1-Ys1*sqrt(sum(I1)),'b--','Linewidth',2);
% plot(T1(1:5),MEK_EGF(1:5),'sb--','Linewidth',2,'MarkerSize',10);%hold off;
% xlim([0 60]);ylim([0 1.05]);
% 
% %figure;hold on;
% plot(timespan,Ym2,'r','Linewidth',2);
% %plot(timespan,Ym2+Ys2*sqrt(sum(I2)),'g--','Linewidth',1);
% %plot(timespan,Ym2-Ys2*sqrt(sum(I2)),'g--','Linewidth',1);
% % plot(timespan,Ym3,'Linewidth',2);
% % plot(MEK_NGF1(:,1),MEK_NGF1(:,2),'og--','Linewidth',2,'MarkerSize',8);
% plot(T2(1:6),MEK_NGF(1:6),'o--','color','r','Linewidth',2,'MarkerSize',10);
% hold off;%plot(timespan,mean(Y21(flag1==1,:)),'Linewidth',2);
% xlim([0 60]);ylim([0 1.1]);
% 
% legend({'Simulation:EGF 100ng/ml','Data:EGF 100ng/ml','Simulation:NGF 200ng/ml','Data:NGF 200ng/ml'},'FontSize',12);
% xlabel('time(minutes)','FontSize',14);
% ylabel('pMEK1/2 (A.U.)','FontSize',14);
% set(gca,'FontSize',14);
% export_fig(gcf,'results/model_fit_pMEK','-jpg','-r300','-q100','-transparent');

%% Pulsed treatment
load Q1;
load V1;



%doses=[0 0.001 0.01 0.03 0.05 0.08 0.1 0.3 0.5 1 5  50]*stimulant/100;
timespan=0:120;
EGF=25;
NGF=50;
% N=size(Q1,1);
N=40;
Y1=zeros(N,length(timespan));
Y2=zeros(N,length(timespan));
flag1=zeros(N,1);
flag2=zeros(N,1);
for i=1:N
    M=SimpleMAPKModel_pulsed;
    M.params=[exp(Q1(i,:)) 1 0];
    M.EGF=EGF;
    M.NGF=NGF;
    M.timespan=timespan;
    y1=simulate_model(M);
    
    if(length(y1)==length(timespan))
        Y1(i,:)=y1(3,:)/max(y1(3,:));
        flag1(i)=1;
    end
    M.params=[exp(Q1(i,:)) 0 1];
    y2=simulate_model(M);
    if(length(y2)==length(timespan))
        Y2(i,:)=y2(3,:)/max(y2(3,:));
        flag2(i)=1;
    end
end

figure;plot(timespan,mean(Y1(flag1==1,:)));hold on; plot(timespan,mean(Y2(flag2==1,:)));legend({'EGF','NGF'});
%%
% [Qm,Qse]=mean_standard_error(Q1,V1);
% Qs=Qse*sqrt(size(Q1,1));%standard deviation;
% 
% %sample from posterior
 Np=100;
% Q11=repmat(Qm,Np,1)+repmat(Qs,Np,1).*randn(Np,size(Qm,2));
% V11=ones(1,Np)/Np;
load Q1;
load V1;

timespan=[0:300];
width=3;
gap=3;
egf.high=25;
egf.low=1;
ngf.high=50;
ngf.low=2;
Il=1:Np;
r_3_3=pulsed_stimulation(Q1(Il,:),V1(Il),timespan,width,gap,egf,ngf);

save('r_3_3.mat','r_3_3');

Xi=0.2*rect_pulse(width,gap,timespan);
I=1:(length(timespan)-1);
Ym1=r_3_3.means(1,:);
Ym2=r_3_3.means(2,:);
Ym11=r_3_3.means(3,:);
Ym21=r_3_3.means(4,:);
figure('Position',[0 0 900 500]);subplot(2,2,1);hold on;plot(timespan(I),Ym1(I),'Linewidth',2);plot(timespan(I),Ym2(I),'Linewidth',2);legend({'pERK(EGF 25ng/ml)','pERK(EGF 1 ng/ml)'},'Location','North','FontSize',12);xlabel('time(minutes)','FontSize',14);ylabel('pERK','FontSize',14);set(gca,'FontSize',14);
subplot(2,2,2);hold on;plot(timespan(I),Ym11(I),'Linewidth',2);plot(timespan(I),Ym21(I),'Linewidth',2); hold off;legend({'pERK(NGF 50ng/ml)','pERK(NGF 2 ng/ml)'},'Location','North','FontSize',12);xlabel('time(minutes)','FontSize',14);ylabel('pERK','FontSize',14);set(gca,'FontSize',14);
subplot(2,2,3);plot(timespan(I),Xi(I),'k');ylim([0,1]);
subplot(2,2,4);plot(timespan(I),Xi(I),'k');ylim([0,1]);
%export_fig(gcf,'model_fit_pulstrain1','-jpg','-r300','-q100','-transparent');

disp('done 3 3...')


%width=3;
gap=10;

r_3_10=pulsed_stimulation(Q1(Il,:),V1(Il),timespan,width,gap,egf,ngf);

save('r_3_10.mat','r_3_10');

Xi=0.2*rect_pulse(width,gap,timespan);
I=1:(length(timespan)-1);
Ym1=r_3_10.means(1,:);
Ym2=r_3_10.means(2,:);
Ym11=r_3_10.means(3,:);
Ym21=r_3_10.means(4,:);
figure('Position',[0 0 900 500]);subplot(2,2,1);hold on;plot(timespan(I),Ym1(I),'Linewidth',2);plot(timespan(I),Ym2(I),'Linewidth',2);legend({'pERK(EGF 25ng/ml)','pERK(EGF 1 ng/ml)'},'Location','North','FontSize',12);xlabel('time(minutes)','FontSize',14);ylabel('pERK','FontSize',14);set(gca,'FontSize',14);
subplot(2,2,2);hold on;plot(timespan(I),Ym11(I),'Linewidth',2);plot(timespan(I),Ym21(I),'Linewidth',2); hold off;legend({'pERK(NGF 50ng/ml)','pERK(NGF 2 ng/ml)'},'Location','North','FontSize',12);xlabel('time(minutes)','FontSize',14);ylabel('pERK','FontSize',14);set(gca,'FontSize',14);
subplot(2,2,3);plot(timespan(I),Xi(I),'k');ylim([0,1]);
subplot(2,2,4);plot(timespan(I),Xi(I),'k');ylim([0,1]);

disp('done 3 10...')


%width=3;
gap=20;

r_3_20=pulsed_stimulation(Q1(Il,:),V1(Il),timespan,width,gap,egf,ngf);

save('r_3_20.mat','r_3_20');

Xi=0.2*rect_pulse(width,gap,timespan);
I=1:(length(timespan)-1);
Ym1=r_3_20.means(1,:);
Ym2=r_3_20.means(2,:);
Ym11=r_3_20.means(3,:);
Ym21=r_3_20.means(4,:);
figure('Position',[0 0 900 500]);subplot(2,2,1);hold on;plot(timespan(I),Ym1(I),'Linewidth',2);plot(timespan(I),Ym2(I),'Linewidth',2);legend({'pERK(EGF 25ng/ml)','pERK(EGF 1 ng/ml)'},'Location','North','FontSize',12);xlabel('time(minutes)','FontSize',14);ylabel('pERK','FontSize',14);set(gca,'FontSize',14);
subplot(2,2,2);hold on;plot(timespan(I),Ym11(I),'Linewidth',2);plot(timespan(I),Ym21(I),'Linewidth',2); hold off;legend({'pERK(NGF 50ng/ml)','pERK(NGF 2 ng/ml)'},'Location','North','FontSize',12);xlabel('time(minutes)','FontSize',14);ylabel('pERK','FontSize',14);set(gca,'FontSize',14);
subplot(2,2,3);plot(timespan(I),Xi(I),'k');ylim([0,1]);
subplot(2,2,4);plot(timespan(I),Xi(I),'k');ylim([0,1]);

disp('done 3 20...')


%width=3;
gap=60;

r_3_60=pulsed_stimulation(Q1(Il,:),V1(Il),timespan,width,gap,egf,ngf);

save('r_3_60.mat','r_3_60');

Xi=0.2*rect_pulse(width,gap,timespan);
I=1:(length(timespan)-1);
Ym1=r_3_60.means(1,:);
Ym2=r_3_60.means(2,:);
Ym11=r_3_60.means(3,:);
Ym21=r_3_60.means(4,:);
figure('Position',[0 0 900 500]);subplot(2,2,1);hold on;plot(timespan(I),Ym1(I),'Linewidth',2);plot(timespan(I),Ym2(I),'Linewidth',2);legend({'pERK(EGF 25ng/ml)','pERK(EGF 1 ng/ml)'},'Location','North','FontSize',12);xlabel('time(minutes)','FontSize',14);ylabel('pERK','FontSize',14);set(gca,'FontSize',14);
subplot(2,2,2);hold on;plot(timespan(I),Ym11(I),'Linewidth',2);plot(timespan(I),Ym21(I),'Linewidth',2); hold off;legend({'pERK(NGF 50ng/ml)','pERK(NGF 2 ng/ml)'},'Location','North','FontSize',12);xlabel('time(minutes)','FontSize',14);ylabel('pERK','FontSize',14);set(gca,'FontSize',14);
subplot(2,2,3);plot(timespan(I),Xi(I),'k');ylim([0,1]);
subplot(2,2,4);plot(timespan(I),Xi(I),'k');ylim([0,1]);

disp('done 3 60...')