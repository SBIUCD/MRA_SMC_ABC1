A=ABC_SMC_WEST;
N=10000;
%Q=zeros(N,A.Nq);
Y=zeros(N,A.Nd);
%MU=ones(1,16)*15;
MU=[5 20 10 20 5 30 5 20 5 15 10 20 1 20 1 20];
Q=A.propose(log(MU),1*ones(1,16),N);
d=zeros(N,1);
S=zeros(N,1);
%count =0;
parfor i=1:N
    M=SimpleMAPKModel;
    M.EGFs=[5 10 20 30 50];
    M.params=exp(Q(i,:));
    [r,y]=r_from_Jacobian_for_all_EGFs(M);
    [r1,y1,s1]=r_from_simulation_for_all_EGFs(M);
    
    %size([r' y'])
    d(i)=sum((r-r1).^2);
    %s1
    S(i)=s1;
    Y(i,:)=[r' y'];
    if s1>100 && d(i)<5
        [s1 d(i)]
    end
    %count =count+1;
    
%       fprintf(repmat('\b',1,numel(S1)));
%       S1=sprintf('%f percent done',count);
%       fprintf(S1);
end

[ds, Ids]=sort(d,'ascend');
Q=Q(Ids,:);

ds(1:5)
%%
save('Q_test.mat','Q');
save('d.mat','ds');

%%
% Find params with the minumum differences between jacobian and simulation
% based r
% Id=(abs(d)<5)&(abs(d)>2);
% Y1=Y(Id,:);
% Q1=Q(Id,:);
% S1=S(Id);
% d1=d(Id);
% %% sort in-terms of difference
% [S2,Is]=sort(S1,'descend');
% % S=S(Is(1:500));
% % Y=Y(Is(1:500),:);
% % Q=Q(Is(1:500),:);
% % d=d(Is);
% S1=S1(Is);
% Y1=Y1(Is,:);
% Q1=Q1(Is,:);
% d1=d1(Is);
% %%
% % I=sum(Y>100,2)==0;
% % Y=Y(I,:);
% % Q=Q(I,:);
% % Yb=Y;
% % Y=(Y-repmat(mean(Y),size(Y,1),1));
% % Y=Y./repmat(max(abs(Y)),size(Y,1),1);
% % Y(isnan(Y))=0;
% 
% 
% % cobj=clustergram(Y,'cluster','column');
% % %%
% % 
% % RL=get(cobj,'RowLabels');
% % save('Yb.mat','Yb');
% % save('Q_test.mat','Q');
% % save('RL.mat','RL');
% % save('d.mat','d');
% % %%
% % M=SimpleMAPKModel;
% % M.params=exp(Q(8902,:));
% % count=1;
% % figure('Position',[0 0 1200 1500]);
% %     for j=1:length(M.EGFs)
% %         M.timespan=timespan;
% %         M.EGF=M.EGFs(j);
% %         
% %         y=simulate_model(M);
% %         %[Ym,Ys]=M.posterior_simulation(M,exp(Q1),V1,timespan);
% %         %size(Ym)
% %         %size(Ys)
% %         subplot(5,3,count);hold on;plot(timespan,y(1,:),'r-','LineWidth',1.5);hold off;%errorbar(timespan,Ym(1,:),Ys(1,:),'b-','LineWidth',1.5);hold off;xlim([0 timespan(end)]); ylim([0 max(max([Ym(1,:)+Ys(1,:);y(1,:)]))+5]);xlabel('time(mins)');ylabel('aRAF');
% %         count=count+1;
% %         subplot(5,3,count);hold on;plot(timespan,y(2,:),'r-','LineWidth',1.5);hold off;%errorbar(timespan,Ym(2,:),Ys(2,:),'b-','LineWidth',1.5);hold off;xlim([0 timespan(end)]); ylim([0 max(max([Ym(2,:)+Ys(2,:);y(2,:)]))+5]);xlabel('time(mins)');ylabel('aMEK');
% %         count=count+1;
% %         subplot(5,3,count);hold on;plot(timespan,y(3,:),'r-','LineWidth',1.5);hold off;%errorbar(timespan,Ym(3,:),Ys(3,:),'b-','LineWidth',1.5);hold off;xlim([0 timespan(end)]); ylim([0 max(max([Ym(3,:)+Ys(3,:);y(3,:)]))+5]);xlabel('time(mins)');ylabel('aERK');
% %         count=count+1;
% % 
% %     end
% % 
