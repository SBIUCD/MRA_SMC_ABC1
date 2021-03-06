classdef SimpleMAPKModel
    %This class defines all propertis and functions needed for describing a
    %simple three stage MAPK model
    properties
        %tot=[20 100 300];
        tot=[100 500 1000];
        y0=zeros(3,1);
        perts=[0.2 0.5 1.5 2];
        EGFs=[0.5 1 2 5];%[0.1 0.5 1 2 5 7.5];%;
        EGF=1;
        kf=0.7; % knockdown factor for perturbation experiments
        timespan=1:30;% time span for simulation
        %params=[5 20 10 20 5 30 5 20 5 15 10 20 5 20 5 20];%model parameters
        params=[]% model parameters
        Y=[];%% Time course simulation results;
        Y_SS=[]; % Seteady state simulation results
        R=[]; % global response coefficients
        r=[]; % local response coefficients
        sigma=0;% percent noise
        RS=2500; % Reporter sensitivity
    end
    
    methods
        % rate equations
        
        %         function [r,y0] = simulate_noisy_local_responses(O)
        %             % simulate local response coefficients from noisy data
        %             NR=6;
        %             R=[];
        %             y0=[];
        %             for i=1:NR
        %             [Ri,y0i]=simulate_global_response_matrix_noisy(O);
        %             R=[R Ri];
        %             y0=[y0 yi];
        %             %O.R=R;
        %             %r=MRA(O);
        %             end
        %
        %         end
        
        %         function [r,y0] = noisy_r_from_simulation_for_all_EGFs(O)
        %             % simulate local response coefficients from noisy data for all
        %             % levels of EGFs
        %             r=[];
        %             y0=[];
        %             for i=1:length(O.EGFs)
        %                 O.EGF=O.EGFs(i);
        %                 [Ri,yi]=simulate_global_response_matrix_noisy(O);
        %                 O.R=Ri;
        %                 ri=MRA(O);
        %                 ri=ri(:);
        %                 r=[r; ri([2:4 6:8])];
        %                 y0=[y0;yi];
        %             end
        %         end
        
        function [r,y0] = simulate_local_responses(O)
            % simulate local response coefficients from noise-free data
            [R,y0]=simulate_global_response_matrix(O);
            O.R=R;
            r=MRA(O);
        end
        
        function [r,y0,S] = r_from_simulation_for_all_EGFs(O)
            % simulate local response coefficients from noise-free data
            r=[];
            y=[];
            for i=1:length(O.EGFs)
                O.EGF=O.EGFs(i);
                [R,y0,S]=simulate_global_response_matrix(O);
                O.R=R;
                ri=MRA(O);
                ri=ri(:);
                r=[r;ri([2:4 6:8])];
                y=[y;y0];
            end
        end
        
        
        function [y] = simulate_model_steadystate_for_all_EGFs(O)
            y=[];
            for i=1:length(O.EGFs)
                O.EGF=O.EGFs(i);
                y=[y simulate_model_steadystate(O)];
            end
        end
        
        function [y] = simulate_model_steadystate(O)
            % Model simulation steady state
            f=@(y) O.dydt(1,y,O);
            y=fsolve(f,O.y0,optimoptions('fsolve','Display','off'));
        end
        
%         function [y] = simulate_model_steadystate(O)
%             % Model simulation steady state
%             f=@(y) O.dydt(1,y,O);
%             y=fsolve(f,O.y0,optimoptions('fsolve','Display','off'));
%         end
        
        function [y] = simulate_model_steadystate_noisy(O)
            % Model simulation steady state
            f=@(y) O.dydt(1,y,O);
            y=fsolve(f,O.y0,optimoptions('fsolve','Display','off'));
            %y1=y;
            y=y*O.RS;
            y=y+randn(size(y)).*O.sigma;
            %[y1 y]
        end
        
        function [Y] = simulate_model_for_all_EGFs(O)
            %Simulate model
            
           % option=odeset('AbsTol',0.01,'RelTol',0.01);
            Y=[];
            for i=1:length(O.EGFs)
                O.EGF=O.EGFs(i);
                y=simulate_model(O);
                Y=[Y;y];
            end
            
        end
        
        function [y] = simulate_model(O)
            %Simulate model
            option=odeset('AbsTol',0.01,'RelTol',0.01);
            [t,y]=ode15s(@O.dydt,O.timespan,O.y0,option,O);
            y=y';
            
        end
        function [y] = simulate_model_steadystate_ode15s(O)
            %Simulate model
            option=odeset('AbsTol',0.01,'RelTol',0.01);
            [t,y]=ode15s(@O.dydt,O.timespan,O.y0,option,O);
            y=y';
            y=y(:,end);
            
        end
        function [R,y0,S] = simulate_global_response_matrix(O)
            % simulate global response matrix for the model
            [ss0, ss1]=SimulatePerturbationExperiments(O);
            %R=2*(ss1-ss0)./(ss1+ss0);
            D=ss1-ss0;%Effects of perturbations
            D=D(:);
            S=sum(D.^2);
            R1=log(ss1)-log(ss0);
            
            d=diag(R1)';
            d=repmat(d,length(d),1);
            R=R1./d;
            y0=ss0(:,1);
        end
        
        function [D] = simulate_observed_metric(O)
            NEGF=length(O.EGFs);
            [r,y] =  r_from_Jacobian_for_all_EGFs(O);
            y1=reshape(y,3,NEGF);y0=repmat(y1(:,1),1,NEGF);y1=y1./y0;y1=y1(:,2:NEGF);
            D=[r;y1(:)]';
        end
        
        function [r,y] = simulate_local_responses_from_noisy_data(O)
            % simulate global response matrix for the model
            Nr=6; % Number of replicates
            Ns=length(O.tot); % Number of species
            Ne=length(O.EGFs);
            %R=zeros(Ns,Ns*Nr,Ne);
            r=[];%zeros(Ns,Ns,Ne);
            y=[];%zeros(Ns,Nr,Ne);
            %SS0=zeros(Ns,Ns*Nr,Ne);
            %SS1=zeros(Ns,Ns*Nr,Ne);
            for j=1:Ne
                Rj=[];
                y0j=[];
                %S0=[];
                %S1=[];
                O.EGF=O.EGFs(j);
                for i=1:Nr
                    [ss0, ss1]=SimulatePerturbationExperiments_noisy(O);
                    R1=log(ss1)-log(ss0);
                    d=diag(R1)';
                    d=repmat(d,length(d),1);
                    Rj=[Rj R1./d];
                    y0j=[y0j ss0(:,1)];
                end
                Rjm=mean(reshape(Rj,Ns,Ns,Nr),3);
                O.R=Rjm;
                r1=MRA(O);r1=r1(:);r1=r1([2:4 6:8]);
                r=[r;r1];
                %SS0(:,:,j)=S0;
                %SS1(:,:,j)=S1;
                %R(:,:,j)=Rj;
                y=[y;mean(y0j,2)];
            end
        end
        
        function [R,y0,SS0,SS1] = simulate_noisy_global_responses_for_all_EGFs(O)
            % simulate global response matrix for the model
            Nr=20; % Number of replicates
            Ns=length(O.tot); % Number of species
            Ne=length(O.EGFs);
            R=zeros(Ns,Ns*Nr,Ne);
            y0=zeros(Ns,Nr,Ne);
            SS0=zeros(Ns,Ns*Nr,Ne);
            SS1=zeros(Ns,Ns*Nr,Ne);
            for j=1:Ne
                Rj=[];
                y0j=[];
                S0=[];
                S1=[];
                O.EGF=O.EGFs(j);
                for i=1:Nr
                    [ss0, ss1]=SimulatePerturbationExperiments_noisy(O);
                    %                       Rij=2*(ss1-ss0)./(ss1+ss0);
                    %                       Rj=[Rj Rij];
                    S0=[S0 ss0];
                    S1=[S1 ss1];
                    R1=log(ss1)-log(ss0);
                    
                    d=diag(R1)';
                    d=repmat(d,length(d),1);
                    Rj=[Rj R1./d];
                    y0j=[y0j ss0(:,1)];
                end
                SS0(:,:,j)=S0;
                SS1(:,:,j)=S1;
                R(:,:,j)=Rj;
                y0(:,:,j)=y0j;
            end
        end
        
        
        function [R,y0] = simulate_global_response_matrix_noisy(O)
            % simulate global response matrix for the model
            
            [ss0, ss1]=SimulatePerturbationExperiments_noisy(O);
            %R=2*(ss1-ss0)./(ss1+ss0);
            R1=log(ss1)-log(ss0);
            
            d=diag(R1)';
            d=repmat(d,length(d),1);
            R=R1./d;
            y0=ss0(:,1);
        end
        
        function [ss0,ss1] = SimulatePerturbationExperiments_noisy(O)
            
            %timespan=1:30;
            y=simulate_model_steadystate_noisy(O);%simulate_model(params,y0,tot,EGF,timespan);
            %ss0=repmat(y(:,end),1,length(tot));
            ss0=repmat(y,1,length(O.tot));
            ss1=[];
            %kf=0.2;
            for i=1:length(O.tot)
                tot1=O.tot;
                tot1(i)=O.tot(i)*O.kf;
                O1=O;
                O1.tot=tot1;
                y1=simulate_model_steadystate_noisy(O1);%simulate_model(params,y0,tot1,EGF,timespan);
                %ss1=[ss1,y1(:,end)];
                ss1=[ss1,y1];
            end
            % Z=2*(ss1-ss0)./(ss1+ss0);
            % Z1=Z\eye(size(Z));
            % r=-pinv(diag(diag(Z1)))*Z1;
        end
        
        function [ss0,ss1] = SimulatePerturbationExperiments(O)
            
            %timespan=1:30;
            y=simulate_model_steadystate(O);%simulate_model(params,y0,tot,EGF,timespan);
            %ss0=repmat(y(:,end),1,length(tot));
            ss0=repmat(y,1,length(O.tot));
            ss1=[];
            %kf=0.2;
            for i=1:length(O.tot)
                tot1=O.tot;
                tot1(i)=O.tot(i)*O.kf;
                O1=O;
                O1.tot=tot1;
                y1=simulate_model_steadystate(O1);%simulate_model(params,y0,tot1,EGF,timespan);
                %ss1=[ss1,y1(:,end)];
                ss1=[ss1,y1];
            end
            % Z=2*(ss1-ss0)./(ss1+ss0);
            % Z1=Z\eye(size(Z));
            % r=-pinv(diag(diag(Z1)))*Z1;
        end
        
        function [ Cr] = MRA(O)
            %Calculates the regulatory coefficients using MRA
            R_1=O.R\eye(size(O.R));
            Cr=-((diag(diag(R_1)))\eye(size(R_1)))*R_1;
        end
        
        function [ J,y] = Jacobian(O)
            %Calculates the Jacobian of the model
            
            y=simulate_model_steadystate(O);
            %y=y(:,end);
            
            y=y(:,end);
            
            kf4=O.params(1);
            Kmf4=O.params(2);
            Vm4=O.params(3);
            Km4=O.params(4);
            
            kf5=O.params(5);
            Kmf5=O.params(6);
            Vm5=O.params(7);
            Km5=O.params(8);
            
            kf6=O.params(9);
            Kmf6=O.params(10);
            Vm6=O.params(11);
            Km6=O.params(12);
            
            kr4=O.params(13);
            Kmr4=O.params(14);
            kr5=O.params(15);
            Kmr5=O.params(16);
            
            RF0=O.tot(1);
            MK0=O.tot(2);
            EK0=O.tot(3);
            
            aRF=y(1);
            aMK=y(2);
            aEK=y(3);
            
            iRF=RF0-y(1);
            iMK=MK0-y(2);
            iEK=EK0-y(3);
            
            J=zeros(3);
            J(1,1)=-O.dMMa_dsub(kr4,Kmr4,aRF,aEK)-O.dMM0_dsub(Vm4,Km4,aRF);
            J(1,3)=-O.dMMa_dmod(kr4,Kmr4,aRF,aEK);
            J(2,1)=O.dMMa_dmod(kf5,Kmf5,iMK,aRF);
            J(2,2)=-O.dMMa_dsub(kr5,Kmr5,aMK,aEK)-O.dMM0_dsub(Vm5,Km5,aMK);
            J(2,3)=-O.dMMa_dmod(kr5,Kmr5,aMK,aEK);
            J(3,2)=O.dMMa_dmod(kf6,Kmf6,iEK,aMK);
            J(3,3)=-O.dMM0_dsub(Vm6,Km6,aEK);
            
        end
        
        function [r,y] = r_from_Jacobian(O)
            % Caculate r from Jacobian
            [J,y]=Jacobian(O);
            % r=J;
            
            s=repmat(y',length(y),1).*repmat(1./y,1,length(y));
            r1=repmat(diag(J),1,size(J,2));
            r=-J./r1;
            r=r.*s;
        end
        
        function [r,y] = r_from_Jacobian_for_all_EGFs(O)
            % Caculate r from Jacobian
            r=[];
            y=[];
            for i=1:length(O.EGFs)
                O.EGF=O.EGFs(i);
                [ri,yi]=r_from_Jacobian(O);
                ri=ri(:);
                r=[r;ri([2:4 6:8])];
                y=[y;yi];
            end
            
        end
        function [r,y] = r_for_all_EGFs(O)
            % Caculate r from Jacobian
            r=[];
            y=[];
            for i=1:length(O.EGFs)
                O.EGF=O.EGFs(i);
                [ri,yi]=simulate_local_responses(O);
                ri=ri(:);
                r=[r;ri([2:4 6:8])];
                y=[y;yi];
            end
        end
        
        
        function [r,y0] = Estimate_r_FromMultipleNoisyPerturbationExperiments_all_EGFs(O)
            r=[];
            y0=[];
            for i=1:length(O.EGFs)
                O.EGF=O.EGFs(i);
                [ri,S0]=EstimateLocalResponsesFromMultipleNoisyPerturbationExperiments(O);
                ri=ri(:);
                r=[r;ri([2:4 6:8])];
                %size(S0)
                y0=[y0;100*S0./(O.tot'*O.RS)];
            end
            
        end
        
        
        
        
        
        
        function [r,y0] = EstimateLocalResponsesFromMultipleNoisyPerturbationExperiments(O)
            Nr=6;
            S0=[];
            S1=[];
            for i=1:Nr
                [S0i,S1i]=SimulateMultiplePerturbationExperiments_noisy(O);
                S0(:,:,i)=S0i;
                S1(:,:,i)=S1i;
            end
            
            S00=mean(S0,3);
            S11=mean(S1,3);
            
            S1=log(S11);
            S0=log(S00);
            Rx=zeros(3);
            for i=1:3
                Ii=i:3:12;
                X1=S1(:,Ii);
                X2=S0(:,Ii);
                % average of mutiple three point numerical derivation using
                % lagrange polynomial
                
                X1_1=[X1(:,1) X2(:,1) X1(:,4)];
                X1_2=[X1(:,2) X2(:,1) X1(:,3)];
                X1_3=[X1(:,1) X2(:,1) X1(:,3)];
                X1_4=[X1(:,2) X2(:,1) X1(:,4)];
                Ii=setdiff(1:3,i);
                d21_1=three_point_numerical_derivative(X1_1(Ii(1),:),X1_1(i,:));
                d21_2=three_point_numerical_derivative(X1_2(Ii(1),:),X1_2(i,:));
                d21_3=three_point_numerical_derivative(X1_3(Ii(1),:),X1_3(i,:));
                d21_4=three_point_numerical_derivative(X1_4(Ii(1),:),X1_4(i,:));
                d21=(d21_1+d21_2+d21_3+d21_4)/4;
                
                d31_1=three_point_numerical_derivative(X1_1(Ii(2),:),X1_1(i,:));
                d31_2=three_point_numerical_derivative(X1_2(Ii(2),:),X1_2(i,:));
                d31_3=three_point_numerical_derivative(X1_3(Ii(2),:),X1_3(i,:));
                d31_4=three_point_numerical_derivative(X1_4(Ii(2),:),X1_4(i,:));
                d31=(d31_1+d31_2+d31_3+d31_4)/4;
                
                Rx(i,i)=1;
                Rx(Ii,i)=[d21;d31];
                
            end
            O.R=Rx;
            r=MRA(O);
            y0=S00(:,1);
        end
        
        
        function [r] = EstimateLocalResponsesFromMultiplePerturbationExperiments(O)
            [S0,S1]=SimulateMultiplePerturbationExperiments(O);
            S1=log(S1);
            S0=log(S0);
            Rx=zeros(3);
            for i=1:3
                Ii=i:3:12;
                X1=S1(:,Ii);
                X2=S0(:,Ii);
                X1_1=[X1(:,1) X2(:,1) X1(:,4)];
                X1_2=[X1(:,2) X2(:,1) X1(:,3)];
                X1_3=[X1(:,1) X2(:,1) X1(:,3)];
                X1_4=[X1(:,2) X2(:,1) X1(:,4)];
                Ii=setdiff(1:3,i);
                d21_1=three_point_numerical_derivative(X1_1(Ii(1),:),X1_1(i,:));
                d21_2=three_point_numerical_derivative(X1_2(Ii(1),:),X1_2(i,:));
                d21_3=three_point_numerical_derivative(X1_3(Ii(1),:),X1_3(i,:));
                d21_4=three_point_numerical_derivative(X1_4(Ii(1),:),X1_4(i,:));
                d21=(d21_1+d21_2+d21_3+d21_4)/4;
                
                d31_1=three_point_numerical_derivative(X1_1(Ii(2),:),X1_1(i,:));
                d31_2=three_point_numerical_derivative(X1_2(Ii(2),:),X1_2(i,:));
                d31_3=three_point_numerical_derivative(X1_3(Ii(2),:),X1_3(i,:));
                d31_4=three_point_numerical_derivative(X1_4(Ii(2),:),X1_4(i,:));
                d31=(d31_1+d31_2+d31_3+d31_4)/4;
                
                Rx(i,i)=1;
                Rx(Ii,i)=[d21;d31];
                
            end
            O.R=Rx;
            r=MRA(O);
            
        end
        function [S0,S1] = SimulateMultiplePerturbationExperiments(O)
            S0=[];
            S1=[];
            for i=1:length(O.perts)
                O.kf=O.perts(i);
                [ss0, ss1]=SimulatePerturbationExperiments(O);
                S0=[S0 ss0];
                S1=[S1 ss1];
            end
        end
        
        function [S0,S1] = SimulateMultiplePerturbationExperiments_noisy(O)
            S0=[];
            S1=[];
            for i=1:length(O.perts)
                O.kf=O.perts(i);
                [ss0, ss1]=SimulatePerturbationExperiments_noisy(O);
                S0=[S0 ss0];
                S1=[S1 ss1];
            end
        end
        
        
        
    end
    
    methods(Static)
        function [ y ] = HM( kcat,Km,sub,mod,kcat1,Km1,mod1)
            %Hyperbolic modifier function
            y=((kcat*sub*mod)*(Km1+kcat1*mod1))/((Km+sub)*(Km1+mod1));
        end
        
        function [ y] = MM0(Vm,Km,sub)
            % 0 order Michaelis Menten
            y=Vm*sub/(Km+sub);
            
        end
        
        function [ y ] = MMa( kcat,Km,sub,mod )
            %First order Michaelis Menten
            y=kcat*sub*mod/(Km+sub);
            
        end
        
        function [r] = dMM0_dsub(Vm,Km,sub)
            % Differential of MM0 w.r.t. sub
            r=Vm*Km/((Km+sub)^2);
        end
        
        function [ r] = dMMa_dmod(kcat,Km,sub,mod)
            % Differential of MMa w.r.t. mod
            r=(kcat*sub)/(Km+sub);
        end
        function [ r] = dMMa_dsub(kcat,Km,sub,mod)
            % Differential of MMa w.r.t. sub
            r=kcat*mod*Km/((Km+sub)^2);
        end
        
        function [ r] = dHM_dmod(kcat,Km,sub,mod,kcat1,Km1,mod1)
            % Differential of HM w.r.t. mod
            %r=((kcat*mod*Km)*(Km1+kcat1*mod1))/(((Km+sub)^2)*(Km1+mod1));
            r=(kcat*sub*(Km1+kcat1*mod1))/((Km+sub)*(Km1+mod1));
        end
        
        function [ r] = dHM_dmod1(kcat,Km,sub,mod,kcat1,Km1,mod1)
            % Differential of HM w.r.t. mod1
            %r=((kcat*mod*Km)*(Km1+kcat1*mod1))/(((Km+sub)^2)*(Km1+mod1));
            r1=Km1*(kcat1-1)/((Km1+mod1)^2);
            r=((kcat*sub*mod)/(Km+sub))*r1;
        end
        
        function [ r] = dHM_dsub(kcat,Km,sub,mod,kcat1,Km1,mod1)
            % Differential of HM w.r.t. sub
            r=((kcat*mod*Km)*(Km1+kcat1*mod1))/(((Km+sub)^2)*(Km1+mod1));
        end
        
        function [dydt1] = dydt(t,y,O)
            %Rate function for MAPK model
            
            kf4=O.params(1);
            Kmf4=O.params(2);
            Vm4=O.params(3);
            Km4=O.params(4);
            
            kf5=O.params(5);
            Kmf5=O.params(6);
            Vm5=O.params(7);
            Km5=O.params(8);
            
            kf6=O.params(9);
            Kmf6=O.params(10);
            Vm6=O.params(11);
            Km6=O.params(12);
            
            kr4=O.params(13);
            Kmr4=O.params(14);
            kr5=O.params(15);
            Kmr5=O.params(16);
            
            RF0=O.tot(1);
            MK0=O.tot(2);
            EK0=O.tot(3);
            
            aRF=y(1);
            aMK=y(2);
            aEK=y(3);
            
            
            iRF=RF0-y(1);
            iMK=MK0-y(2);
            iEK=EK0-y(3);
            dydt1=zeros(3,1);
            dydt1(1)=O.MMa(kf4,Kmf4,iRF,O.EGF)-O.MMa(kr4,Kmr4,aRF,aEK)-O.MM0(Vm4,Km4,aRF);
            dydt1(2)=O.MMa(kf5,Kmf5,iMK,aRF)-O.MMa(kr5,Kmr5,aMK,aEK)-O.MM0(Vm5,Km5,aMK);
            dydt1(3)=O.MMa(kf6,Kmf6,iEK,aMK)-O.MM0(Vm6,Km6,aEK);
            dydt1=dydt1(:);
        end
        
        function [Ym,Ys]=posterior_simulation(O,post_params,l,timespan)
            % simulate for a set of parameters
            Y1=[];
            Y2=[];
            Y3=[];
            l1=[];
            EGF=O.EGF;
            
            parfor i=1:size(post_params,1)
                M=SimpleMAPKModel;
                M.EGF=EGF;
                M.params=post_params(i,:);
                M.timespan=timespan;
                y=simulate_model(M);
                y=y/max(y(:));
                %size(y)
                if size(y,2)==length(timespan)
                    Y1=[Y1;y(1,:)];
                    Y2=[Y2;y(2,:)];
                    Y3=[Y3;y(3,:)];
                    l1=[l1 l(i)];
                end
            end
            l1=l1/sum(l1);
            Y1=Y1.*repmat(l1',1,size(Y1,2));
            Ym1=sum(Y1);
            Ys1=sqrt(sum(((Y1-repmat(Ym1,size(Y1,1),1)).^2).*repmat(l1',1,size(Y1,2))))/sqrt(size(Y1,1));
            
            Y2=Y2.*repmat(l1',1,size(Y2,2));
            Ym2=sum(Y2);
            Ys2=sqrt(sum(((Y2-repmat(Ym2,size(Y2,1),1)).^2).*repmat(l1',1,size(Y2,2))))/sqrt(size(Y2,1));
            
            Y3=Y3.*repmat(l1',1,size(Y3,2));
            Ym3=sum(Y3);
            Ys3=sqrt(sum(((Y3-repmat(Ym3,size(Y3,1),1)).^2).*repmat(l1',1,size(Y3,2))))/sqrt(size(Y3,1));
            
            Ym=[Ym1;Ym2;Ym3];
            Ys=[Ys1;Ys2;Ys3];
        end
        
        function [Ym,Ys]=posterior_r_Y(O,post_params,l)
            % simulate for a set of parameters
            Y=zeros(size(post_params,1),54);
            %l1=[];
            parfor i=1:size(post_params,1)
                M=O;
                M.params=post_params(i,:);
                %O.timespan=timespan;
                [ri,yi]=r_from_Jacobian_for_all_EGFs(M);
                Y(i,:)=[ri;yi]';
            end
            l=l/sum(l);
            Y=Y.*repmat(l',1,size(Y,2));
            Ym=sum(Y);
            Ys=sqrt(sum(((Y-repmat(Ym,size(Y,1),1)).^2).*repmat(l',1,size(Y,2))))/sqrt(size(Y,1));
        end
        
        function [Ym,Ys]=posterior_error(O,post_params,l,Yobs)
            % simulate for a set of parameters
            Y=zeros(size(post_params,1),1);
            %l1=[];
            parfor i=1:size(post_params,1)
                M=O;
                M.params=post_params(i,:);
                %O.timespan=timespan;
                %[ri,yi]=r_from_Jacobian_for_all_EGFs(M);
                Yi=simulate_observed_metric(M);%[ri;yi];  
                %size(Yi)
                %size(Yobs)
                Y(i)=sum((Yobs-Yi).^2);
                
            end
            l=l/sum(l);
            Y=Y.*repmat(l',1,size(Y,2));
            Ym=sum(Y);
            Ys=sqrt(sum(((Y-repmat(Ym,size(Y,1),1)).^2).*repmat(l',1,size(Y,2))))/sqrt(size(Y,1));
        end
              
        
    end
end

