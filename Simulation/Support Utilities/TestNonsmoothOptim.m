clear all
clc

global rgbtriplets
rgbtriplets = [ [0.0000 0.4470 0.7410];...
    [0.8500 0.3250 0.0980];...
    [0.4660 0.6740 0.1880];...
    [0.9290 0.6940 0.1250];...
    [0.3010 0.7450 0.9330];...
    [0.7000 0.1840 0.7760];...
    [0.9900 0.0000 0.0000]];

N       = 2000;
p       = 5;
X       = randn(N,p);
y       = zeros(N,1);

muX     = mean(X,1);
stdX    = std(X,[],1);
X       = (X-muX).*(3*rand(1,p));%./stdX;

global beta1 beta2 betaall
beta1   = zeros(p,1);
beta1([1,2,3]) = randn(3,1);
% beta1 = beta1/norm(beta1);

beta2   = zeros(p,1);
beta2([1,3,5]) = randn(3,1);
% beta2 = null(beta1'); beta2 = beta2(:,3);
% beta2 = beta2/norm(beta2);

betaall = [beta1,beta2];

nspc = null([beta1,beta2]');
beta3 = nspc(:,1);

NumCorupt = round(0.02*N);
y(1:N/2)    = (X(1:N/2,:) )*beta1;
y(1+N/2:end)= (X(1+N/2:end,:))*beta2;
% y(1:NumCorupt) = 2*randn(NumCorupt,1);

%%

% RegressWeight.Name      = '';
% RegressWeight.Name      = 'tukey';
%     RegressWeight.csq       = 4.6852^2;
% RegressWeight.Name      = 'huber';
%     RegressWeight.csq       = 1.3449^2;
% RegressWeight.Name      = 'welsch';
%     RegressWeight.csq       = 2.9843^2;
% RegressWeight.Name      = 'cauchy';
%     RegressWeight.csq       = 2.3848^2;
% RegressWeight.Name      = 't-robust';


RegressMethod.Name      = 'ls';
% RegressMethod.Name      = 'lars';
%     RegressMethod.
% RegressMethod.Name      = 'pls';
%     RegressMethod.XThreshold    = 0.99;
%     RegressMethod.YThreshold    = 0.99;
%     RegressMethod.NumComp       = p;

% FactorizationMethod.Name = 'tpm';
FactorizationMethod.Name = 'rpsf';

Option.NumMix           = 2;
% Option.RegressWeight    = RegressWeight;
Option.RegressMethod    = RegressMethod;
Option.FactorizationMethod = FactorizationMethod;

% Solution = Tensor_Chaganty2013(y,X,Option);
Solution = NonsmoothOptim_Bagirov2013(y,X,Option);

%%

delta = sum((beta1-Solution.beta).^2,1);
b1    = delta(1) < delta(2);

figure(1)
stem(beta1,'b.'); hold on
if b1
    plot(1:p,Solution.beta(:,1),'bo'); hold off
else
    plot(1:p,Solution.beta(:,2),'bo'); hold off
end
box off

figure(2)
stem(beta2,'r.'); hold on
if b1
    plot(1:p,Solution.beta(:,2),'ro'); hold off
else
    plot(1:p,Solution.beta(:,1),'ro'); hold off
end
box off

figure(3)
if b1
    yhat(1:round(N/2),1) = X(1:round(N/2),:)*Solution.beta(:,1);
    yhat(1+round(N/2):N,1) = X(1+round(N/2):N,:)*Solution.beta(:,2);
    plot(y(1:round(N/2)),yhat(1:round(N/2)),'b.',y(1+round(N/2):N),yhat(1+round(N/2):N),'r.')
else
    yhat(1:round(N/2),1) = X(1:round(N/2),:)*Solution.beta(:,2);
    yhat(1+round(N/2):N,1) = X(1+round(N/2):N,:)*Solution.beta(:,1);
    plot(y(1:round(N/2)),yhat(1:round(N/2)),'r.',y(1+round(N/2):N),yhat(1+round(N/2):N),'b.')
end
box off
grid on

