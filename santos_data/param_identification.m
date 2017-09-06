A=ABC_SMC_WEST;
N=1000000;
%Q=zeros(N,A.Nq);
Y=zeros(N,A.Nd);
MU=ones(1,16)*10;
%MU=[2 20 10 20 2 20 10 20 2 20 10 20 2 20 2 20];
%MU=[2.0000 6.7397 2.0000 44.7311 1.9815 20.7545 0.4179 3.9714 2.7497 1.2514 20.0000 10.0000 0.7853 3.8448 0.7264 90.0000]; 
Q=A.propose(log(MU),2*ones(1,16),N);
TOT=zeros(N,3);
Y=zeros(N,10);
d=zeros(N,1);
S=zeros(N,1);
flag=zeros(N,1);
count =0;
parfor i=1:N
    M=SimpleMAPKModel;
    M.EGF=15.62; % EGF molecular weigth =6400 dalton. concentration was 100 ng/ml = 100000 ng/L=100000/6400 nM/L=15.62 nM/L
    M.tot=exp(log([50 50 50])+1*randn(1,3));
    M.timespan=1:10;
    M.params=exp(Q(i,:));
    y=simulate_model(M);
    %size(y)
    if size(y,2)==10
        y5=y(3,5);
        ym=max(y(3,:),[],2);
        %disp(ym);
        if y5==ym && y5>1.1*y(3,end)
            disp('****************************found********************* \n');
            flag(i)=1;
            Y(i,:)=y(3,:);
            TOT(i,:)=M.tot;
            %count =count+1;
            %disp(count);
        end
    end
end
%%
Y1=Y(flag==1,:);
Q1=Q(flag==1,:);
T1=TOT(flag==1,:);
I=sum(Y1<0,2)==0;

Y1=Y1(I,:);
T1=T1(I,:);
Q1=Q1(I,:);

I1=sum((Y1>=repmat(T1(:,3),1,10)),2)==0;
Y1=Y1(I1,:);
T1=T1(I1,:);
Q1=Q1(I1,:);

save('Y1.mat','Y1');
save('Q1.mat','Q1');
save('T1.mat','T1');
%% Further narrow down the parameter space
flag1=zeros(size(Q1,1),1);


parfor i=1:size(Q1,1)
    M=SimpleMAPKModel;
    M.EGF=15.62; % EGF molecular weigth =6400 dalton. concentration was 100 ng/ml = 100000 ng/L=100000/6400 nM/L=15.62 nM/L
    M.tot=T1(i,:);
    M.timespan=0:0.25:10;
    M.params=exp(Q1(i,:));
    y=simulate_model(M);
    d1=abs((y(1,22)-y(1,20))/0.5);
    d2=abs((y(2,22)-y(2,20))/0.5);
    if d1<0.001 && d2<0.001
        flag1(i)=1;
    end
end

sum(flag1)
%%


