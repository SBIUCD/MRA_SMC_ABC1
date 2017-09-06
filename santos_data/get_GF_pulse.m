function [EGFt,NGFt]= get_GF_pulse(t,O)
            kd_egf = O.k1;
            kd_ngf = O.k2;
            width=O.width;
            gap=O.gap;
            T=width+gap;
            t1=mod(t,T);
            t1(t1==0)=T;
            I1=t1<=width;
            I2=t1>width;
            
            EGFt=t1;            
            EGFt(I1)=O.EGF*exp(-kd_egf*(t1(I1)-1));
            %EGFt(I2)=O.EGF*exp(-kd_egf*(t1(I2)-width));
            EGFt(I2)=0;
            
            NGFt=t1;            
            NGFt(I1)=O.NGF*exp(-kd_ngf*(t1(I1)-1));
            NGFt(I2)=0;
            %NGFt(I2)=O.NGF*exp(-kd_ngf*(t1(I2)-width));     
            
        end