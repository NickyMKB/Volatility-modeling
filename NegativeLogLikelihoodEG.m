function [negativeLL]=NegativeLogLikelihoodEG(parameter_vector,returns)

%% Extract the stuff we need from the input arguments
mu    = parameter_vector(1,1);
theta_4 = parameter_vector(2,1);
theta_5 = parameter_vector(3,1);
theta_6 = parameter_vector(4,1);
theta_7 = parameter_vector(5,1);
theta_8 = parameter_vector(6,1);
theta_9 = parameter_vector(7,1);
theta_10 = parameter_vector(8,1);
nu = parameter_vector(9,1);
T     = size(returns,1);

%% Run the GARCH filter
[ sigma,~,~] = DynamicScalerTwoEG(mu,theta_4,theta_5,theta_6,theta_7,theta_8,theta_9,theta_10,returns);

% Collect a row vector of log likelihood per observation (this is the log
% of the pdf of a normal distribution)
LL = - log(sigma(1:T)) + log( gamma( (nu+1) /2 ) / ( gamma( nu/2 ) * sqrt(pi*(nu-2)) )  * ...
    ( 1 + (((returns - mu)./sigma(1:T)).^2 /(nu-2) )).^-((nu+1)/2) );

% Put a negative sign in front and sum over all obserations
negativeLL = - sum( LL(1:end) )    ;            

% Close the function
end

