classdef ABC_SMC_WEST
    %This class implements a generic version of the ABC_SMC algorith with
    %adapative weights as suggested by Mike West.
    
    properties
        %e=[75,60,50,40,30,25,20,15,10,8,6,5,4,3.5,3,2.5,2.25,2,1.75,1.5,1.35,1.25,1.125,1,0.95:-0.025:0]; %Threshold schedule
%         e=[15,10,7,6.5,6,5.5,5,4.5,4,3.5,3,2.75,2.5:-0.1:1,0.95:-0.025:0]; %Threshold schedule
       e=[5,4.5,4,3.5,3,2.75,2.5:-0.1:1,0.95:-0.025:0]; %Threshold schedule
        N=1000;% lower bound on the number of samples in each population, 
        Nq=17;% Dimention of the paramteres;
        Nb=512;%Blocks of parameters that will be evaluated in parallel
        Nd=45;% Dimension of the data
        prior_params=struct('mu',log(5)*ones(1,12),'sig',2*ones(1,12),'Indexes',1:11);%prior distribution of the parameter
        sigma=1;% Noise level;
        target=[];% the true paramters
        Q1_name='Q1.mat';
        V1_name='V1.mat';
        D_obs=[];
        D_params=[];
        timespan=[];
        stimulants=[];
    end
%     properties(Dependent)
%         Q_t1; % Sampled parameters in interation t-1;
%         Q_t;% Sampled paramteres in iteration t;
%         M; % The MAPK model
%         D_obs; % Observed data
%         W; % Initial weights
%     end
    
    methods
       
 % The abc_smc_west algorithm       
        function [Q1,V1]=adaptive_weights_abc_smc(O)
            % population t=1
%             M=SimpleMAPKModel;
%             M.sigma=O.sigma;
%             M.params=O.target;
            %[r,y]=r_from_Jacobian_for_all_EGFs(M);
            %D_params=eye(length(D_obs));
            W1=ones(1,O.N)/O.N;
            [Q1,d1]=O.find_N_params_within_e_ball(O.D_obs,O.D_params,O.e(1),O);% find the parameters within the first e-ball
            %d1
            %size(W1)
            %size(d1)
            V1=W1.*O.convert_distance_into_kernel(d1');
            V1=V1/sum(V1);
            fprintf('\r at error level %f \r',O.e(1));
            for i=2:length(O.e)
                [Q2,d2,interrupted]=O.find_N_params_within_e_ball_from_prev_samples(Q1,V1,O.D_obs,O.D_params,O.e(i),O);
                
%                 size(Q1)
%                 size(Q2)
%                 size(V1)
                if interrupted
                    fprintf('\n Search interrupted at error level %f due to slow convergence.. \n',O.e(i));
                    break;
                else
                W1=O.calculate_weights(Q1,Q2,V1,O);
                %[d2';O.convert_distance_into_kernel(d2')]
                V1=W1.*O.convert_distance_into_kernel(d2');
                V1=V1/sum(V1);
                Q1=Q2;
                save(O.Q1_name,'Q1');%Save the last set of params
                save(O.V1_name,'V1');%Save the last set of params
                fprintf('\r at error level %f \r',O.e(i));
                end
            end
            %size(V1)
            
            
        end       
        
    end
    
    
    methods(Static)
        
        
            function [Q1,V1]=adaptive_weights_abc_smc_restart_from_saved_parameters(O,Q1,V1,e)
            % This function continues running the sampling algorithms from
            % last saved parameters, Q1 = last saved parameters, V1=weights
            % of the parameters, e= new error schedule
            for i=1:length(e)
                [Q2,d2,interrupted]=O.find_N_params_within_e_ball_from_prev_samples(Q1,V1,O.D_obs,O.D_params,e(i),O);
                
%                 size(Q1)
%                 size(Q2)
%                 size(V1)
                if interrupted
                    fprintf('\n Search interrupted at error level %f due to slow convergence.. \n',e(i));
                    break;
                else
                W1=O.calculate_weights(Q1,Q2,V1,O);
                %[d2';O.convert_distance_into_kernel(d2')]
                V1=W1.*O.convert_distance_into_kernel(d2');
                V1=V1/sum(V1);
                Q1=Q2;
                save(O.Q1_name,'Q1');%Save the last set of params
                save(O.V1_name,'V1');%Save the last set of params
                fprintf('\r at error level %f \r',e(i));
                end
            end
            %size(V1)
            
            
        end

        
        
        
        function [k]=GaussianKernel(X,X_obs,Sig_t,diag)
            %Returns the Gaussian functions in the form of matrix, if you
            %need only the diagonal element set diag =true;
            x1=X-X_obs;
            if diag
            k=diag(exp(-x1'*Sig_t*x1));
            else
            k=exp(-x1'*Sig_t*x1);
            end
        end
        
        function [V_t1]=data_based_weight(W_t1,Kxt)
            %Calculate the data based weight vector
            V_t1=W_t1.*Kxt;
        end
        
        function [W_t]=new_weight(P_Q,V_t1,K_Qt)
            %Calculate the second weight vector
            W_t=P_Q./(V_t1*K_Qt);
            W_t=W_t/sum(W_t);
        end
        function [Q_s]=propose(Q,Sig,N)
            % Proposes a new Qs based on Gaussian kernel;Assumes Sig is a vector
            % Q is a row vector;Note that Q(13) and Q(15) are uniformely
            % distributed between (0,1) and therefore needs to be sampled
            % accordingly
           
            if size(Q,1)==1
                Q_s=[repmat(Q,N,1)+randn(N,length(Q)).*repmat(Sig,N,1)];
            else
                Q_s=[Q+randn(N,size(Q,2)).*repmat(Sig,N,1)];
            end
        end
        
        function [p]=prior(vars,prior_params)
            % log normal prior, assumes that vars is in log-space
        %p=exp(-sum((vars-prior_params.mu).^2)./(2*(prior_params.sig^2)));
        Nv=size(vars,1);
        M=repmat(prior_params.mu,Nv,1);
        V=(2*repmat(prior_params.sig,Nv,1).^2);
        p=exp(-sum(((vars-M).^2)./V ,2));
        
        end
        
        function [d]=distance(O,O_obs,D_params)
            % log normal prior, assumes that vars is in log-space
        I=1:length(O_obs);%[7:12 22:27];%1:length(O_obs);%[7:15 25:33];%[1:15 19:33];%1:length(O_obs);%[7:15 25:33];%[1:6 16:18 19:24 34:36];%1:length(O_obs);%[1:6 16:18 19:24 34:36];%1:length(O_obs);%[1:15 19:33];
        %I=19:24;
        O1=O(I)-O_obs(I);
        %O(I)
        d=sqrt(O1*D_params(I,I)*O1');
        end
        
        function [d]=sample_N_within_error_threshold(O)
            % log normal prior, assumes that vars is in log-space
        O1=O-O_obs;
        d=sqrt(O1*D_params*O1');
        end
        
        function [d,count]=simulate_model_outputs(Q,D_obs,D_params,e,O)
            % simulates the r and steady states response for a block of Nb
            % paramteres, to save time Nb should be multiple of the numbers
            % of cores active in your parpool
            Nb=size(Q,1); % Blocksize
            %Ot=zeros(Nb,O.Nd); % Initialize output
            ME=max(O.e)+1;
            d=ME+zeros(1,Nb);
            count=0;
            parfor i=1:Nb
                M=SimpleMAPKModel;
                M.EGF=O.stimulants(1);%;100;%15.62;
                M.NGF=O.stimulants(2);
                M.timespan=O.timespan;
                p=exp(Q(i,:));
                Oi=M.simulate_observed_metric(M,p);
                if length(Oi)>1
                    d(i)= O.distance(Oi,D_obs,D_params);
                    if d(i)<=e
                        count=count+1;
                    end
                end
                %                 if count==Nr
                %                     break;
                %                 end
            end
        end
        
        function [Qs]=sample(Q,W,Ns)
            % Draws N samples for parameter set Q according to
            % probabilities in W; We assume that W is normalized and a row
            % vector
            Wc=repmat(cumsum(W),Ns,1); % replicate the cumulative sum Ns times
            r=repmat(rand(Ns,1),1,length(W)); % generate random number
            %size(r)
            %size(Wc)
            Is=sum(Wc<r,2)+1; % indexes of samples    
            Qs=Q(Is,:); % samples
        end
        
        function [K]=convert_distance_into_kernel(d)
            
            c1=prctile1(d,10); % rule of thumb
            if c1==0
                I=d~=0;
                d1=d(I);
                c1=prctile1(d1,10);
            end
            K=exp(-(d.^2)./(2*(c1^2)));% we use gaussian kernel, we chan change it to different kernels
        end
        
        function [W_t]= calculate_weights(Q_t1,Q_t,v_t1,O)
            %Calculate transition kernel
            Q_1=Q_t1;% Take only those paramteres which have Gaussian Kernel, the last two paramters have uniform kernel
            Q0=Q_t;% Take only those paramteres which have Gaussian Kernel, the last two paramters have uniform kernel
            C=diag(var(Q_1));
            D=(pdist2(Q_1,Q0,'mahalanobis',C).^2)/2;
            D=exp(-D);
            P_t=O.prior(Q0,O.prior_params);
            %size(P_t)
            %size(v_t1*D)
            W_t=O.new_weight(P_t',v_t1,D);
            
        end
        
        
         function [Q_s,d]=find_N_params_within_e_ball(D_obs,D_params,e,O)
             %O is the object
             count =0; % numbers of hits (cases within the e-ball)
             Nb=O.Nb ;% calculate stuff in blocks of 64 for the purpose of parallelizing
             Q_s=zeros(O.N,O.Nq); % Initialize the output matrix
             d=zeros(O.N,1);%Initialize the distance metrics
             S='';
             while count<O.N % Until all N, has been sampled
                 Qtb=O.propose(O.prior_params.mu,O.prior_params.sig,Nb); % propose a block of parameters from the prior distribution
                 [d1,count1]=O.simulate_model_outputs(Qtb,D_obs,D_params,e,O); % simulate model output and calculate distance from the observed values
                 if count+count1<=O.N % the last samples do not fullfill the total rewuirement
                     Ie=find(d1<=e); % Indexes of samples withing e ball
                     Q_s((count+1):(count+length(Ie)),:)=Qtb(Ie,:); % Store the samples within e ball
                     d((count+1):(count+length(Ie)))=d1(Ie); % store the corresponding distance metrics
                     count=count+count1; % increment the count
                 else % if more than required samples are within e ball
                     Nr=O.N-count; % count how many are required 
                     Ie=find(d1<=e);% Indexes of samples withing e ball
                     Ie=Ie(1:Nr); % Take only as many as required
                     Q_s((count+1):(count+Nr),:)=Qtb(Ie,:); % Store the samples
                     d((count+1):(count+Nr))=d1(Ie);% store the distances
                     count=count+count1; % increase count
                 end % end-if
                 fprintf(repmat('\b',1,numel(S)));
                 S=sprintf('%d of %d samples found',count,O.N);
                 fprintf(S);
             end % end while
         end % end-func
         
         
         function [Q_s,d,interrupted]=find_N_params_within_e_ball_from_prev_samples(Qs,W,D_obs,D_params,e,O)
             %O is the object
             % Here Qs is previous sample and W is the weight vector
             
             count =0; % numbers of hits (cases within the e-ball)
             total_count=0; % number of params tested
             interrupted=false; % was the loop interrupted because of poor acceptance rates
             Nb=O.Nb ;% calculate stuff in blocks of 64 for the purpose of parallelizing
             Q_s=zeros(O.N,O.Nq); % Initialize the output matrix
             sig=std(Qs)/(O.N^(1/6));% standard deviation of the input samples Qs
             d=zeros(O.N,1);%Initialize the distance metrics
             S='';
             %disp(sig);
             while count<O.N % Until all N, has been sampled
                 Qs1=O.sample(Qs,W,Nb);
                 Qtb=O.propose(Qs1,sig,Nb); % propose a block of parameters from the prior distribution
                 [d1,count1]=O.simulate_model_outputs(Qtb,D_obs,D_params,e,O); % simulate model output and calculate distance from the observed values
                 
                 
                 if count+count1<=O.N % the last samples do not fullfill the total requirement
                     Ie=find(d1<=e); % Indexes of samples withing e ball
                     Q_s((count+1):(count+length(Ie)),:)=Qtb(Ie,:); % Store the samples within e ball
                     d((count+1):(count+length(Ie)))=d1(Ie); % store the corresponding distance metrics
                     count=count+count1; % increment the count
                     total_count=total_count+Nb;% increase the number of total parameter sets tested
                     acc_rate=count*100/total_count;%calculate acceptance rate
                     if total_count>50000 && acc_rate<=0.01
                       interrupted=true;
                       break;
                     end
                 else % if more than required samples are within e ball
                     Nr=O.N-count; % count how many are required 
                     Ie=find(d1<=e);% Indexes of samples withing e ball
                     Ie=Ie(1:Nr); % Take only as many as required
                     Q_s((count+1):(count+Nr),:)=Qtb(Ie,:); % Store the samples
                     d((count+1):(count+Nr))=d1(Ie);% store the distances
                     count=count+count1; % increase count                     
                 end % end-if
                 fprintf(repmat('\b',1,numel(S)));
                 S=sprintf('%d of %d samples found',count,O.N);
                 fprintf(S);
             end % end while
         end % end-func
        
    end
    
        
        
    end
    
    
    


