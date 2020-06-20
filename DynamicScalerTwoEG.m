function [ sigma,short_run,long_run ] = DynamicScalerTwoEG(mu, theta_4, theta_5, theta_6, theta_7, theta_8, theta_9,theta_10,returns)
%Instead of capturing long memory by a factionally integrated process, two
%components may be used

% Extract the sample size (make sure returns are a column vector)
T     = size(returns,1);
% Run the dynamic scale
for t=1:(T+1)
    if t==1;
        % Initialise sigma at unconditional mean of variance
        sigma(t,1)     = sqrt(var(returns));
        %Initalise short_run at 0 , whereas the long-term starts at
        %unconditional mean of variance
        short_run(t,1) = 0;
        long_run(t,1) = log(sigma(t,1));
    else
        residual = (returns(t-1,1) - mu)/sigma(t-1,1); 
        
        short_run(t,1) = theta_4*short_run(t-1,1) + theta_5*residual + theta_6*(abs(residual) - sqrt(2/pi));
        
        long_run(t,1) = theta_7 + theta_8*long_run(t-1,1) + theta_9*residual + theta_10*(abs(residual) - sqrt(2/pi));
        
        sigma(t,1) = exp(short_run(t,1) + long_run(t,1));
        
    end
% Close the function
end