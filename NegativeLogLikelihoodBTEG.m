function [negativeLL]=NegativeLogLikelihoodBTEG(parameter_vector,returns)

%% Extract the stuff we need from the input arguments
mu    = parameter_vector(1,1);
lambdaS = parameter_vector(2,1);
phi = parameter_vector(3,1);
kappa = parameter_vector(4,1);
kappa_tilde  = parameter_vector(5,1);
nu = parameter_vector(6,1);
T     = size(returns,1);

%% Run the GARCH filter
[ sigma ] = exp(DynamicScaler(mu,lambdaS,phi,kappa,kappa_tilde,nu,returns));

% Collect a row vector of log likelihood per observation (this is the log
% of the pdf of a t distribution)
LL = - log(sigma(1:T)) + log( gamma( (nu+1) /2 ) / ( gamma( nu/2 ) * sqrt(pi*(nu-2)) )  * ...
    ( 1 + (((returns - mu)./sigma(1:T)).^2 /(nu-2) )).^-((nu+1)/2) );

% Put a negative sign in front and sum over all obserations
negativeLL = - sum( LL(1:end) )    ;            

% Close the function
end

