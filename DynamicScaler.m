function [ lambda , residuals ] = DynamicScaler(mu,lambdaS,phi,kappa,kappa_tilde,nu,returns)
% Extract the sample size (make sure returns are a column vector)
T     = size(returns,1);
% Run the dynamic scale
for t=1:(T+1)
    if t==1;
        % Initialise (lambda starting point will be ln(var(return)))
        lambda(t,1)     = lambdaS;      
    else
        residual = (returns(t-1,1) - mu)*exp(-lambda(t-1,1) ); 
        residuals(t,1) = residual;
        lambda(t,1) = lambdaS*(1-phi) + phi*lambda(t-1,1) + ... 
            kappa*(sqrt(nu+3)/sqrt(nu*2))*((nu+1)/(nu-2+ residual^2))*(residual^2 -1 ) + ...
            kappa_tilde*( sqrt( (nu-2)*(nu+3) )/ sqrt(nu*(nu+1)) )*((nu+1)/(nu-2+ residual^2))*residual ;
    end
    end
% Close the function
end