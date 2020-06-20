function [NIC] = NICBTEG(kappa_ML,kappastar_ML,nu_ML,residuals)
%Instead of capturing long memory by a factionally integrated process, two
%components may be used
T= size(residuals,1);
% Extract the sample size (make sure returns are a column vector)
% Run the dynamic scale
for t=1:(T)
    if t==1;
        NIC(t,1)=0;
    else
    NIC(t,1) = 2*kappa_ML* (sqrt(nu_ML+3)/sqrt(nu_ML*2))*((nu_ML+1)/(nu_ML-2+ residuals(t,1)^2))*(residuals(t,1)^2 -1 ) + ...
+2*kappastar_ML* ( sqrt( (nu_ML-2)*(nu_ML+3) )/ sqrt(nu_ML*(nu_ML+1)) )*((nu_ML+1)/(nu_ML-2+ residuals(t,1)^2))*residuals(t,1);
end

end
% Close the function
end