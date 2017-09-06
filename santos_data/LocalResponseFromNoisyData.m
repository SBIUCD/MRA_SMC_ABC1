classdef LocalResponseFromNoisyData
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods(Static)
        function [r]=multiple_MRA(R,O)
            Nr=size(R,2)/size(R,1); % Number of replicate experiments
            Ns=size(R,1); % number of species
            r=zeros(Ns);
            for i=1:Nr
                is=(i-1)*Ns+1;
                ie=i*Ns;
                Ri=R(:,is:ie);
                r=r+O.MRA(Ri);
            end
            r=r/Nr;
        end
        
        function [r]=multiple_MRA_for_all_EGFs(R,O)
            r=[];
            for i=1:size(R,3)
                Ri=R(:,:,i);
                ri=O.multiple_MRA(Ri,O);
                ri=ri(:);
                r=[r;ri([2:4 6:8])];                
            end
        end
        
        function [r]=MRA_on_average_global_response(R,O)
            Nr=size(R,2)/size(R,1); % Number of replicate experiments
            Ns=size(R,1); % number of species
            R1=zeros(Ns,Ns,Nr);
            for i=1:Nr
                is=(i-1)*Ns+1;
                ie=i*Ns;
                Ri=R(:,is:ie);
                R1(:,:,i)=Ri;
            end
            r=O.MRA(mean(R1,3));
        end
        
         function [r]=MRA_on_averaged_data(SS0,SS1,O)
            Nr=size(SS0,2)/size(SS0,1); % Number of replicate experiments
            Ns=size(SS0,1); % number of species
            S0=zeros(Ns,Ns,Nr);
            S1=zeros(Ns,Ns,Nr);
            for i=1:Nr
                is=(i-1)*Ns+1;
                ie=i*Ns;
                
                S0i=SS0(:,is:ie);
                S0(:,:,i)=S0i;
                S1i=SS1(:,is:ie);
                S1(:,:,i)=S1i;
                %R1(:,:,i)=Ri;
            end
            S0m=mean(S0,3);
            S1m=mean(S1,3);
            
            R1=log(S1m)-log(S0m);
            
            d=diag(R1)';
            d=repmat(d,length(d),1);
            R1=R1./d;
            r=O.MRA(R1);
         end
        
        
         function [r]=MRA_on_averaged_data_for_all_EGFs(SS0,SS1,O)
            r=[];
            for i=1:size(SS0,3)
                ri=O.MRA_on_averaged_data(SS0(:,:,i),SS1(:,:,i),O);
                ri=ri(:);
                r=[r;ri([2:4 6:8])];
            end
        end
         
         function [r]=BootstrapMRA(R,O)
            Nr=size(R,2)/size(R,1); % Number of replicate experiments
            Ns=size(R,1); % number of species
            R1=zeros(Ns,Ns,Nr);
            for i=1:Nr
                is=(i-1)*Ns+1;
                ie=i*Ns;
                Ri=R(:,is:ie);
                R1(:,:,i)=Ri;
            end
            Rm=mean(R1,3);
            Rs=std(R1,[],3);
            Nsam=100000;
            r=zeros(Ns);
            for i=1:Nsam
                Ri= Rm+randn(size(Rm)).*Rs;
                r=r+O.MRA(Ri);
            end
            r=r/Nsam;
         end
        
        
         function [r]=Bootstrap_MRA_for_all_EGFs(R,O)
            r=[];
            for i=1:size(R,3)
                ri=O.BootstrapMRA(R(:,:,i),O);
                ri=ri(:);
                r=[r;ri([2:4 6:8])];
            end
        end
         
        function [r]=MRA_on_average_global_response_for_all_EGFs(R,O)
            r=[];
            for i=1:size(R,3)
                ri=O.MRA_on_average_global_response(R(:,:,i),O);
                ri=ri(:);
                r=[r;ri([2:4 6:8])];
            end
        end
        
        
        
        
        function [ Cr] = MRA(R)
            %Calculates the regulatory coefficients using MRA
            R_1=R\eye(size(R));
            Cr=-((diag(diag(R_1)))\eye(size(R_1)))*R_1;
        end
        
    end
    
end

