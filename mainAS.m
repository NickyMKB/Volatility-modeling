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

% startingvalues = [ omega , sigma^2 , alpha_pos,alpha_neg , beta , nu]
startingvalues=[mean(SSreturns);var(SSreturns)/20;0.1;0.1;0.80;13];

%% Check the negative log likelihood at the starting values
NegativeLogLikelihoodAS(startingvalues,SSreturns)

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

% Parameter lower bound and upper bound (note that sigma^2,alpha_pos,alpha_neg,beta must be
% positive and nu >2)
lowerbound = [-inf,0,0,0,0,2];
upperbound = [inf,inf,inf,inf,inf,inf];

% Perform ML maximisation (we actually minimize the negative likelihood)
[ML_parameters,ML_NegativeLogL,~,~,~,~,hessian]=fmincon('NegativeLogLikelihoodAS', startingvalues ,[],[],[],[],lowerbound,upperbound,[],options,SSreturns )

%% Check if the optimal negative log likelihood is better (lower) than the starting point 
% This number should be positive
NegativeLogLikelihoodAS(startingvalues,SSreturns) - ML_NegativeLogL

%% Save the parameters and compute GARCH filter at these parameters

mu_ML    = ML_parameters(1);
omega_ML = ML_parameters(2);
alpha_pos_ML = ML_parameters(3);
alpha_neg_ML = ML_parameters(4);
beta_ML  = ML_parameters(5);
nu_ML  = ML_parameters(6);

T=size(returns,1);
se_ML =  sqrt(diag(inv(hessian)));
[ sigmasquared_ML , residuals] = ASGarchFilter(mu_ML,omega_ML,alpha_pos_ML,alpha_neg_ML,beta_ML,returns);

NIC = NICAS(alpha_pos_ML,alpha_neg_ML,residuals);

[prediction_5_AG ] = pvol_AGARCH(mu_ML,omega_ML,alpha_pos_ML,alpha_neg_ML,beta_ML,5,sigmasquared_ML,returns);
[prediction_21_AG ] = pvol_AGARCH(mu_ML,omega_ML,alpha_pos_ML,alpha_neg_ML,beta_ML,21,sigmasquared_ML,returns);

pred_5(:,1) = prediction_5_AG(:,1);
pred_5(:,2) = sqrt(5/250)*vix(:,1);

pred_21(:,1) = prediction_21_AG(:,1);
pred_21(:,2) = sqrt(21/250)*vix(:,1);

T = size(returns,1)
for t=1:T
    if(t< T- 6)
pred_5(t,3) = sqrt(sum(returns(t+1:t+6).^2));
pred_5(t,4) = sqrt(4/3 * sum(rv5(t+1:t+6)));
if(t < T -21)
pred_21(t,3) = sqrt(sum(returns(t+1:t+22).^2));
pred_21(t,4) = sqrt(4/3 * sum(rv5(t+1:t+22)));
end
    end
end

%% Plot the absolute returns along with one times the GARCH-implied stdev

figure
plot(residuals(2:5115),NIC(2:5115),'ok')

figure 
plot(dates,abs(returns-mu_ML),'ok')
hold on 
plot(dates,2*sqrt(sigmasquared_ML(1:T)),'b-','Linewidth',2)
legend({'$|r_t-\mu|$','$2\times \widehat{\sigma}_{t}$'},'Interpreter','latex')
set(gca, 'XTick', datenum(['20000101';'20010101';'20020101';'20030101';'20040101';'20050101';'20060101';'20070101';'20080101';'20090101';'20100101';'20110101';'20120101';'20130101';'20140101';'20150101';'20160101';'20170101';'20180101';'20190101';'20200101';'20210101'],'yyyymmdd') )
dateFormat = 'yy';
datetick('x', dateFormat,'keepticks')
xtickangle(0)
%axis([datenum('20000101','yyyymmdd') datenum('20210101','yyyymmdd') -inf inf])
set(gca,'FontName','Times','fontsize',15)
set(gca,'TickDir','out')

%% Check that the implied shocks have approximately the desired characteristics (mean zero, variance one)

implied_epsilon=(returns-mu_ML)./sqrt(sigmasquared_ML(1:T));

[mean(implied_epsilon);var(implied_epsilon)]
