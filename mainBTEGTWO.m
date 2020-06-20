%% Clear all data 
clear
close all

%% Load the directory (change this to where you have saved this file)
%pad = 'C:/Users/61849rla/Dropbox/Basis Weco/';
pad = '/Users/Nicky/Documents/Case3';
cd(pad);
load('data')

%% Get dates in correct Matlab format
dates  = num2str(dates);
dates  = datenum(dates,'yyyymmdd');
%% Subset data
SSreturns = returns(1:5114); 

%% Specify starting values for the paramsters mu, omega, alpha, beta

% startingvalues = [ lambda: sqrt of long variance returns ,, nu]
%omega,phi1,kappa1,kappastar1,phi2,kappa2,kappastar2,nu,returns
startingvalues=[log(var(SSreturns))/20; 0.1 ; 0.9 ; -0.05 ;0.5; 0.1 ; -0.5 ; 3 ; mean(SSreturns)];

%% Check the negative log likelihood at the starting values
NegativeLogLikelihoodBTEGTWO(startingvalues,SSreturns)

%% Do maximum likelihood optimsiation 
% Matlab likes to minimize, so we minimise the *negative* log likelihood

% Clear any pre-existing options
clearvars options

% Load some options
options  =  optimset('fmincon');
options  =  optimset(options , 'TolFun'      , 1e-6);
options  =  optimset(options , 'TolX'        , 1e-6);
options  =  optimset(options , 'Display'     , 'on');
options  =  optimset(options , 'Diagnostics' , 'on');
options  =  optimset(options , 'LargeScale'  , 'off');
options  =  optimset(options , 'MaxFunEvals' , 10^6) ;
options  =  optimset(options , 'MaxIter'     , 10^6) ;

%omega,phi1,kappa1,kappastar1,phi2,kappa2,kappastar2,nu,returns free to roam ; nu >2
lowerbound = [-inf,-inf,-inf,-inf,-inf,-inf,-inf,2,-inf];
upperbound = [inf,inf,inf,inf,inf,inf,inf,inf,inf];

% Perform ML maximisation (we actually minimize the negative likelihood)
[ML_parameters,ML_NegativeLogL,~,~,~,~,hessian]=fmincon('NegativeLogLikelihoodBTEGTWO', startingvalues ,[],[],[],[],lowerbound,upperbound,[],options,SSreturns )

%% Check if the optimal negative log likelihood is better (lower) than the starting point 
% This number should be positive
NegativeLogLikelihoodBTEGTWO(startingvalues,SSreturns) - ML_NegativeLogL

%% Save the parameters and compute GARCH filter at these parameters

omega_ML = ML_parameters(1);
phi1_ML=ML_parameters(2);
kappa1_ML=ML_parameters(3);
kappastar1_ML=ML_parameters(4);
phi2_ML= ML_parameters(5);
kappa2_ML=ML_parameters(6);
kappastar2_ML=ML_parameters(7);
nu_ML=ML_parameters(8);
mu_ML= ML_parameters(9);
se_ML =  sqrt(diag(inv(hessian)));

[ lambda_ML , lambda1_ML, lambda2_ML, residuals ] = DynamicScalerTwo(omega_ML,phi1_ML,kappa1_ML,kappastar1_ML,phi2_ML,kappa2_ML,kappastar2_ML,nu_ML,mu_ML,returns);

[sigma_ML] = exp(lambda_ML);
T = size(returns,1); 
%% Plot the absolute returns along with one times the GARCH-implied stdev

figure 
plot(dates,abs(returns-mu_ML),'ok')
hold on 
plot(dates,2*sigma_ML(1:T),'b-','Linewidth',2)
legend({'$|r_t-\mu|$','$2\times \widehat{\sigma}_{t}$'},'Interpreter','latex')
set(gca, 'XTick', datenum(['20000101';'20010101';'20020101';'20030101';'20040101';'20050101';'20060101';'20070101';'20080101';'20090101';'20100101';'20110101';'20120101';'20130101';'20140101';'20150101';'20160101';'20170101';'20180101';'20190101';'20200101';'20210101'],'yyyymmdd') )
dateFormat = 'yy';
datetick('x', dateFormat,'keepticks')
xtickangle(0)
%axis([datenum('20000101','yyyymmdd') datenum('20210101','yyyymmdd') -inf inf])
set(gca,'FontName','Times','fontsize',15)
set(gca,'TickDir','out')

%% Plot NIC
NIC1 = NICBTEG(kappa1_ML,kappastar1_ML,nu_ML,residuals); %long term
NIC2 = NICBTEG(kappa2_ML,kappastar2_ML,nu_ML,residuals); %short term
figure 
plot(residuals(2:5115),NIC1(2:5115),'ok')
figure 
plot(residuals(2:5115),NIC2(2:5115),'ok')


%% Check that the implied shocks have approximately the desired characteristics (mean zero, variance one)

implied_epsilon=(returns-mu_ML)./(sigma_ML(1:T));

[mean(implied_epsilon);var(implied_epsilon)]

