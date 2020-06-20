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
SSreturns = returns(1:4112); 
%% Specify starting values for the paramsters mu, omega, alpha, beta

% startingvalues = [ mu, theta_4, theta_5, theta_6, theta_7, theta_8, theta_9,theta_10,nu,returns]
startingvalues=[mean(SSreturns);0.333;-0.07;-0.02;-0.02;0.989;-0.032;0.061;6];

%% Check the negative log likelihood at the starting values
NegativeLogLikelihoodEG(startingvalues,SSreturns)

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

% Parameter lower bound and upper bound ( mu, theta_4, theta_5, theta_6, theta_7, theta_8, theta_9,theta_10
lowerbound = [-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,2];
upperbound = [inf,inf,inf,inf,inf,inf,inf,inf,inf,];

% Perform ML maximisation (we actually minimize the negative likelihood)
[ML_parameters,ML_NegativeLogL,~,~,~,~,hessian]=fmincon('NegativeLogLikelihoodEG', startingvalues ,[],[],[],[],lowerbound,upperbound,[],options,SSreturns )

%% Check if the optimal negative log likelihood is better (lower) than the starting point 
% This number should be positive
NegativeLogLikelihoodEG(startingvalues,SSreturns) - ML_NegativeLogL

%% Save the parameters and compute GARCH filter at these parameters

mu_ML    = ML_parameters(1,1);
theta_4_ML = ML_parameters(2,1);
theta_5_ML = ML_parameters(3,1);
theta_6_ML = ML_parameters(4,1);
theta_7_ML = ML_parameters(5,1);
theta_8_ML = ML_parameters(6,1);
theta_9_ML = ML_parameters(7,1);
theta_10_ML = ML_parameters(8,1);
nu_ML = ML_parameters(9,1);
se_ML =  sqrt(diag(inv(hessian)));

[ sigma_ML ] = DynamicScalerTwoEG(mu_ML, theta_4_ML, theta_5_ML, theta_6_ML, theta_7_ML, theta_8_ML, theta_9_ML,theta_10_ML,returns);

T= size(returns,1);


T = size(returns,1)
for t=1:T
    if(t<= T)
        
pred_5(t,1) = sqrt(returns(t).^2);
pred_5(t,2) = sqrt(4/3 *rv5(t));
end
end

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

%% Check that the implied shocks have approximately the desired characteristics (mean zero, variance one)

implied_epsilon=(returns-mu_ML)./sigma_ML(1:T);

[mean(implied_epsilon);var(implied_epsilon)]
