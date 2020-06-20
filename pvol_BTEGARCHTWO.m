function[prediction_h] = pvol_BTEGARCH(omega,phi1,kappa1,kappastar1,phi2,kappa2,kappastar2,nu,mu,returns)

T = size(returns,1);
for t = 1:T
    sum = 0;
        for d =1:h
            temp = 2*omega + 2*phi1^(d-1)*(lambda1(t+1,1) - omega ) + ...
    2*(kappa^2 + kappa_tilde^2)* ( (1-phi^(2*d-2) ) )/ (1- phi^2);
            sum = sum + exp(temp);
        end
prediction_h(t,1) = sqrt( h*mu^2 + sum);
sum
end
