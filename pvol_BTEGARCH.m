function[prediction_h] = pvol_BTEGARCH(mu,lambdaS,phi,kappa,kappa_tilde,h,lambda,returns)

T = size(returns,1);
for t = 1:T
    sum = 0;
        for d =1:h 
            temp = 2*lambdaS + 2*phi^(d-1)*(lambda(t+1,1) - lambdaS ) + ...
    2*(kappa^2 + kappa_tilde^2)* ( (1-phi^(2*d-2) ) )/ (1- phi^2);
            sum = sum + exp(temp);
        end
prediction_h(t,1) = sqrt( h*mu^2 + sum);
sum
end
