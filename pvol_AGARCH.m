function[prediction_h] = pvol_AGARCH(mu,omega,alpha_pos,alpha_neg,beta,h,sigmasquared,returns)

T = size(returns,1);

for t = 1:(T)
            
prediction_h(t,1) = sqrt ( h*(mu^2 + sigmasquared(1)) + ...
    ( (1- (alpha_pos/2 + alpha_neg/2 + beta)^h )/ ( 1 - alpha_pos/2 - alpha_neg/2 - beta ) ) * ...
     ( sigmasquared(t+1) - sigmasquared(1) ) );
end