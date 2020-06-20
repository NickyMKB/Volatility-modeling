function [ lambda,lambda1,lambda2 , residuals] = DynamicScalerTwo(omega,phi1,kappa1,kappastar1,phi2,kappa2,kappastar2,nu,mu,returns)
%Instead of capturing long memory by a factionally integrated process, two
%components may be used

% Extract the sample size (make sure returns are a column vector)
T     = size(returns,1);
% Run the dynamic scale
for t=1:(T+1)
    if t==1;
        % Initialise (lambda starting point will be ln(var(return)))
        lambda(t,1)     = omega; 
        %Initalise lambda1,lamba2 TODO
        lambda1(t,1) = 0;
        lambda2(t,1) = 0;
    else
        residual = (returns(t-1,1) - mu)*exp(-lambda(t-1,1) ); 
        residuals(t,1) = residual;
        
        lambda1(t,1) = phi1 * lambda1(t-1,1) + ...
        kappa1* (sqrt(nu+3)/sqrt(nu*2))*((nu+1)/(nu-2+ residual^2))*(residual^2 -1 ) + ...
        +kappastar1* ( sqrt( (nu-2)*(nu+3) )/ sqrt(nu*(nu+1)) )*((nu+1)/(nu-2+ residual^2))*residual ;
    
        lambda2(t,1) = phi2 * lambda2(t-1,1) + ...
        kappa2* (sqrt(nu+3)/sqrt(nu*2))*((nu+1)/(nu-2+ residual^2))*(residual^2 -1 ) + ...
        +kappastar2* ( sqrt( (nu-2)*(nu+3) )/ sqrt(nu*(nu+1)) )*((nu+1)/(nu-2+ residual^2))*residual ;
        
        lambda(t,1) = omega + lambda1(t,1) + lambda2(t,1);
    
        
    end
% Close the function
end