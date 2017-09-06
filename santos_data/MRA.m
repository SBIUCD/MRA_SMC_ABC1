function [ Cr] = MRA(R)
            %Calculates the regulatory coefficients using MRA
            R_1=R\eye(size(R));
            Cr=-((diag(diag(R_1)))\eye(size(R_1)))*R_1;
        end

