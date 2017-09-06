function [r] = pulsed_stimulation( Q1,V1,timespan,width,gap,egf,ngf)
Y1=zeros(size(Q1,1),length(timespan));
Y2=zeros(size(Q1,1),length(timespan));
flag1=zeros(size(Q1,1),1);
flag2=zeros(size(Q1,1),1);
parfor i=1:size(Q1,1)
    M=SimpleMAPKModel_pulsed;
    M.EGF=egf.high;
    %M.n=1;
    M.pulse.width=width;
    M.pulse.gap=gap;
    M.timespan=timespan;
    M.NGF=0;%ngf_doses(j);
    M.y0=[0 0 0];    
    % Simulate for EGF
    M.params=[exp(Q1(i,:)) 1 0];
    y0=simulate_model(M);  
   if size(y0,2)==length(M.timespan)
       Y1(i,:)=y0(3,:);%/max(y0(3,:)); % ERK level at 5 min        
       flag1(i)=1;
   end 
   
   M.EGF=egf.low;
    y0=simulate_model(M);  
   if size(y0,2)==length(M.timespan)
       Y2(i,:)=y0(3,:);%/max(y0(3,:)); % ERK level at 5 min        
       flag2(i)=1;
   end
   
  
end


Y11=zeros(size(Q1,1),length(timespan));
Y21=zeros(size(Q1,1),length(timespan));
flag11=zeros(size(Q1,1),1);
flag21=zeros(size(Q1,1),1);
parfor i=1:size(Q1,1)
    M=SimpleMAPKModel_pulsed;
    %M.n=1;
    M.EGF=0;
    M.timespan=timespan;
    M.NGF=ngf.high;%ngf_doses(j);
    M.y0=[0 0 0];    
    % Simulate for EGF
    M.pulse.width=width;
    M.pulse.gap=gap;
    M.params=[exp(Q1(i,:)) 0 1];
    y0=simulate_model(M);  
   if size(y0,2)==length(M.timespan)
       Y11(i,:)=y0(3,:);%/max(y0(3,:)); % ERK level at 5 min       
       flag11(i)=1;
   end 
   
   M.NGF=ngf.low;
   y0=simulate_model(M);
   if size(y0,2)==length(M.timespan)
       Y21(i,:)=y0(3,:);%/max(y0(3,:)); % ERK level at 5 min       
       flag21(i)=1;
   end
  
end
%%
I1=flag1==1;
Y1_1=Y1(I1,:);
V1_1=V1(I1);
%Yf1=pERK_2_FRET.pERK2FRET_timecourse_bulk(Y1_1);
%[Ym1,Ys1]=mean_standard_error(Yf1,V1_1);
[Ym1,Ys1]=mean_standard_error(Y1_1,V1_1);

I2=flag2==1;
Y2_1=Y2(I2,:);
V2_1=V1(I2);
%Yf2=pERK_2_FRET.pERK2FRET_timecourse_bulk(Y2_1);
%[Ym2,Ys2]=mean_standard_error(Yf2,V2_1);
[Ym2,Ys2]=mean_standard_error(Y2_1,V2_1);

I11=flag11==1;
Y1_2=Y11(I11,:);
V1_2=V1(I11);
%Yf11=pERK_2_FRET.pERK2FRET_timecourse_bulk(Y1_2);
%[Ym11,Ys11]=mean_standard_error(Yf11,V1_2);
[Ym11,Ys11]=mean_standard_error(Y1_2,V1_2);

I21=flag21==1;
Y2_2=Y21(I21,:);
V2_2=V1(I21);
%Yf21=pERK_2_FRET.pERK2FRET_timecourse_bulk(Y2_2);
%[Ym21,Ys21]=mean_standard_error(Yf21,V2_2);
[Ym21,Ys21]=mean_standard_error(Y2_2,V2_2);

% r.EGF_high=Yf1;
% r.EGF_low=Yf2;
% r.NGF_high=Yf11;
% r.NGF_low=Yf21;

if size(Ym1,2)==1
    Ym1=zeros(size(timespan));
end
if size(Ym2,2)==1
    Ym2=zeros(size(timespan));
end
if size(Ym11,2)==1
    Ym11=zeros(size(timespan));
end
if size(Ym21,2)==1
    Ym21=zeros(size(timespan));
end

if size(Ys1,2)==1
    Ys1=zeros(size(timespan));
end
if size(Ys2,2)==1
    Ys2=zeros(size(timespan));
end
if size(Ys11,2)==1
    Ys11=zeros(size(timespan));
end
if size(Ys21,2)==1
    Ys21=zeros(size(timespan));
end

r.means=[Ym1;Ym2;Ym11;Ym21];
r.stds=[Ys1;Ys2;Ys11;Ys21];
r.Y1=Y1_1;
r.Y2=Y2_1;
r.Y11=Y1_2;
r.Y22=Y2_2;
%time=[0 2 10 20 40 60 120 180 240 300]+1;


end

