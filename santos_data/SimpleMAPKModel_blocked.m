classdef SimpleMAPKModel_blocked
    %This class defines all propertis and functions needed for describing a
    %simple three stage MAPK model
    properties
       
        y0=zeros(3,1);
       
       
        EGF=0;
        NGF=0;
       
        timespan=1:30;% time span for simulation
        %params=[5 20 10 20 5 30 5 20 5 15 10 20 5 20 5 20];%model parameters
        params=[];
        pulse=struct('width',0,'gap',0);
    end
    
    methods
       
        
     function [y,y1] = simulate_model(O)
            option=odeset('AbsTol',0.000001,'RelTol',0.000001);
            [t,y]=ode15s(@O.dydt,O.timespan,O.y0,option,O);
            y=y';     
           
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
             %Rate function for MAPK model
            kf11=O.params(1);
            Kmf11=O.params(2);
            
            kf12=O.params(3);
            Kmf12=O.params(4);
            
            kf13=O.params(5);
            Kmf13=O.params(6);
            
            kf14=O.params(7);
            Kmf14=O.params(8);
            
            Vm1=O.params(9);
            Km1=O.params(10);
            
            kf21=O.params(11);
            Kmf21=O.params(12);
            
            kf22=O.params(13);
            Kmf22=O.params(14);
            
            Vm2=O.params(15);
            Km2=O.params(16);
            
            kf31=O.params(17);
            Kmf31=O.params(18);
            
            kf32=O.params(19);
            Kmf32=O.params(20);
            
            Vm3=O.params(21);
            Km3=O.params(22);            
            
            kd_egf=O.params(23);
            kd_ngf=O.params(24);
            
            RF0=O.params(25);
            MK0=O.params(26);
            EK0=O.params(27);
            
            
            
            be=O.params(28);
            bn=O.params(29);
            
            
            aRF=y(1);
            aMK=y(2);
            aEK=y(3);
            
            
            iRF=RF0-aRF;
            iMK=MK0-aMK;
            iEK=EK0-aEK;
            
            %[iRF iMK iEK]
            dydt1=zeros(3,1);
            
            %EGF=(O.EGF/(1+O.EGF))*t*exp(-kd_egf*t);
            %NGF=(O.NGF/(1+O.NGF))*t*exp(-kd_ngf*t);
            [EGF,NGF]=O.get_GF_pulse(t,O);
            dydt1(1)=be*(O.MMa(kf11,Kmf11,iRF,EGF)-O.MMa(kf12,Kmf12,aRF,aEK)) + bn*(O.MMa(kf13,Kmf13,iRF,NGF) + O.MMa(kf14,Kmf14,iRF,aEK))- O.MM0(Vm1,Km1,aRF);
            dydt1(2)=O.MMa(kf21,Kmf21,iMK,aRF)-O.MMa(kf22,Kmf22,aMK,aEK)-O.MM0(Vm2,Km2,aMK);
            dydt1(3)=O.MMa(kf31,Kmf31,iEK,aMK) + bn*O.MMa(kf32,Kmf32,iEK,aRF) - O.MM0(Vm3,Km3,aEK);

%             dydt1(1)=be*(O.MMa(kf11,Kmf11,iRF,EGF)*(Kmf12/(Kmf12+aEK))) + bn*(O.MMa(kf13,Kmf13,iRF,NGF) + O.MMa(kf14,Kmf14,iRF,aEK))- O.MM0(Vm1,Km1,aRF);
%             dydt1(2)=O.MMa(kf21,Kmf21,iMK,aRF)*(Kmf22/(Kmf22+aEK))-O.MM0(Vm2,Km2,aMK);
%             dydt1(3)=O.MMa(kf31,Kmf31,iEK,aMK) + bn*O.MMa(kf32,Kmf32,iEK,aRF) - O.MM0(Vm3,Km3,aEK);
%             dydt1(1)=be*O.MMa(kf11,Kmf11,iRF,O.EGF)+bn*O.MMa(kf13,Kmf13,iRF,O.NGF);
%             dydt1(2)=0;%O.MMa(kf21,Kmf21,iMK,aRF);
%             dydt1(3)=0;%O.MMa(kf31,Kmf31,iEK,aMK);
            dydt1=dydt1(:);
            %dydt1
        end
        
        function [Ym,Ys]=posterior_simulation(O,post_params,l,timespan)
            % simulate for a set of parameters
            Y1=[];
            Y2=[];
            Y3=[];
            l1=[];
            for i=1:size(post_params,1)
                O.params=post_params(i,:);
                O.timespan=timespan;
                y=simulate_model(O);
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
        
        
        function [r,y] = r_from_Jacobian_custom(O,y)
            % Caculate r from Jacobian
            [J,y]=O.Jacobian_custom(O,y);
            % r=J;
            
            s=repmat(y',length(y),1).*repmat(1./y,1,length(y));
            r1=repmat(diag(J),1,size(J,2));
            r=-J./r1;
            r=r.*s;
        end
        function [ J,y] = Jacobian_custom(O,y)
            %Calculates the Jacobian of the model
                        %Rate function for MAPK model
            kf11=O.params(1);
            Kmf11=O.params(2);
            
            kf12=O.params(3);
            Kmf12=O.params(4);
            
            kf13=O.params(5);
            Kmf13=O.params(6);
            
            kf14=O.params(7);
            Kmf14=O.params(8);
            
            Vm1=O.params(9);
            Km1=O.params(10);
            
            kf21=O.params(11);
            Kmf21=O.params(12);
            
            kf22=O.params(13);
            Kmf22=O.params(14);
            
            Vm2=O.params(15);
            Km2=O.params(16);
            
            kf31=O.params(17);
            Kmf31=O.params(18);
            
            kf32=O.params(19);
            Kmf32=O.params(21);
            
            Vm3=O.params(22);
            Km3=O.params(23);            
            
            kd_egf=O.params(24);
            kd_ngf=O.params(25);
            
            RF0=O.params(26);
            MK0=O.params(27);
            EK0=O.params(28);
            
            be=O.params(29);
            bn=O.params(30);
            
            
            aRF=y(1);
            aMK=y(2);
            aEK=y(3);
            
            
            iRF=RF0-y(1);
            iMK=MK0-y(2);
            iEK=EK0-y(3);
            
%             dydt1(1)=be*(O.MMa(kf11,Kmf11,iRF,O.EGF)-O.MMa(kf12,Kmf12,aRF,aEK)) + bn*(O.MMa(kf13,Kmf13,iRF,O.NGF) + O.MMa(kf14,Kmf14,iRF,aEK))- O.MM0(Vm1,Km1,aRF);
%             dydt1(2)=O.MMa(kf21,Kmf21,iMK,aRF)-O.MMa(kf22,Kmf22,aMK,aEK)-O.MM0(Vm2,Km2,aMK);
%             dydt1(3)=O.MMa(kf31,Kmf31,iEK,aMK) + bn*O.MMa(kf32,Kmf32,iEK,aRF) - O.MM0(Vm3,Km3,aEK);
%             dydt1=dydt1(:);
            
            J=zeros(3);
            J(1,1)=-be*O.dMMa_dsub(kf12,Kmf12,aRF,aEK)-O.dMM0_dsub(Vm1,Km1,aRF);
            J(1,3)=-be*O.dMMa_dmod(kf12,Kmf12,aRF,aEK)+bn*O.dMMa_dmod(kf14,Kmf14,iRF,aEK);
            J(2,1)=O.dMMa_dmod(kf21,Kmf21,iMK,aRF);
            J(2,2)=-O.dMMa_dsub(kf22,Kmf22,aMK,aEK)-O.dMM0_dsub(Vm2,Km2,aMK);
            J(2,3)=-O.dMMa_dmod(kf22,Kmf22,aMK,aEK);
            J(3,1)=bn*O.dMMa_dmod(kf32,Kmf32,iEK,aRF);
            J(3,2)=O.dMMa_dmod(kf31,Kmf31,iEK,aMK);
            J(3,3)=-O.dMM0_dsub(Vm3,Km3,aEK);
            
        end
        
        function [J]=Jacobian_custom_numerical(O,y,t)
            J=[];
            d=0.0001;
            for i=1:length(y)
                yp=y;yn=y;
                yp(i)=yp(i)+d;
                yn(i)=yn(i)-d;
                dyp=O.dydt(t,yp,O);
                dyn=O.dydt(t,yn,O);
                Ji=(dyp-dyn)/(2*d);
                J=[J Ji];
            end
        end
        
        function [EGFt,NGFt]= get_GF_pulse(t,O)
            EGFt=0;
            NGFt=0;
            if(t<=10)
                kd_egf=O.params(23);
                kd_ngf=O.params(24);
                EGFt=(O.EGF/(1+O.EGF))*(t^5)*exp(-kd_egf*t);
                NGFt=(O.NGF/(1+O.NGF))*(t^5)*exp(-kd_ngf*t);
            end
        end
    end
        
        
    end


