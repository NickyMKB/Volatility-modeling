function [ sigmasquared, residuals ] = ASGarchFilter(mu,omega,alpha_pos,alpha_neg,beta,returns)
% Extract the sample size (make sure returns are a column vector)
T     = size(returns,1);
% Run the GARCH filter
for t=1:(T+1)
    if t==1;
        % Initialise at the unconditional mean of sigmasquared
        sigmasquared(t,1)     = omega/(1-alpha_pos/2-alpha_neg/2-beta) ;      
    else
        if (returns(t-1,1)-mu)>=0;
            %inidicator functions returns 1
            indicator =1;
        else
            indicator = 0;
        end
        residuals(t,1) = (returns(t-1,1)-mu)/sqrt(sigmasquared(t-1,1));
        sigmasquared(t,1) = omega + (alpha_pos*indicator + alpha_neg*(1-indicator) )* (returns(t-1,1)-mu)^2 + ...
            beta* sigmasquared(t-1,1);
    end
    end
% Close the function
end