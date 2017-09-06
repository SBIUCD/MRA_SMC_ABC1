function [X]= rect_pulse(width,gap,t)
            
           
            T=width+gap;%period
            t1=mod(t,T);
            t1(t1==0)=T;
            I1=t1<=width;
            I2=t1>width;
            X=zeros(1,length(t));
            X(I1)=1;
            X(I2)=0;
        end

