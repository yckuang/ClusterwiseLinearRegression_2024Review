clear
clc

TargetEff = 0.95;

% x   = linspace(0,6,1e8);
x   = linspace(0,6,1e7);

%% Huber

% c0  = linspace(1.34,1.35,100)';
% 
% eff = zeros(length(c0),1);
% parfor index = 1 : length(c0)
%         
%     c	= c0(index);
%     c2  = c^2;
%         
% %     bIn = x < c;
% %     a   = 4/(2*pi)*( trapz(x(bIn),exp(-0.5*x(bIn).^2)) )^2;
% %     b   = 2/sqrt(2*pi)* trapz(x(bIn), (x(bIn).^2).*exp(-0.5*x(bIn).^2));
% %     bOut= x>=c;
% %     b   = b + 2/sqrt(2*pi)* trapz(x(bOut),(c2).*exp(-0.5*x(bOut).^2));
% 
%     a = 4*(normcdf(c)^2-normcdf(c)+0.25);
%     b = (2-2*c2)*normcdf(c) + (2*c2-1) - 2*c*exp(-0.5*c2)/sqrt(2*pi);
%     eff(index) = a/b;
% 
% end

%% Tukey

% c0  = linspace(4.67,4.69,100)';
% 
% a = zeros(length(c0),1);
% b = zeros(length(c0),1);
% parfor index = 1 : length(c0)
%         
%     c	= c0(index);
%     
%     z = x/c;
%     bIn = x < c;
%     a(index)   = ( trapz(x(bIn),(1-6*(z(bIn).^2)+5*(z(bIn).^4)).*exp(-0.5*x(bIn).^2)) )^2;
%     b(index)   = trapz(x(bIn), (x(bIn).^2).*((1-z(bIn).^2).^4).*exp(-0.5*x(bIn).^2));
% 
% end
% eff = sqrt(2/pi)* a./b;

%% Welsch

% c0  = linspace(2.95,3.05,100)';
% 
% a = zeros(length(c0),1);
% b = zeros(length(c0),1);
% parfor index = 1 : length(c0)
% 
%     c	= c0(index);
% 
%     z = x/c;
%     a(index)   = ( trapz(x,(1-2*z.^2).*exp(-z.^2).*exp(-0.5*x.^2)) )^2;
%     b(index)   = trapz(x, (x.^2).*(exp(-2*z.^2)).*exp(-0.5*x.^2));
% 
% end
% eff = sqrt(2/pi)* a./b;


%% Cauchy

% c0  = linspace(2.3,2.4,100)';
% 
% a = zeros(length(c0),1);
% b = zeros(length(c0),1);
% parfor index = 1 : length(c0)
%         
%     c	= c0(index);
%     
%     z = x/c;
%     a(index)   = ( trapz(x,((1-z.^2)./((1+z.^2).^2)).*exp(-0.5*x.^2)) )^2;
%     b(index)   = trapz(x, ((x./(1+z.^2)).^2).*exp(-0.5*x.^2));
% 
% end
% eff = sqrt(2/pi)* a./b;


%% Output

[~,pt] = min(abs(eff-TargetEff));

figure(1)
plot(c0,eff)
grid on
title({['estimation efficiency'];...
       ['Target Efficiency ', num2str(TargetEff*100), '%, c = ', num2str(c0(pt)) ]})