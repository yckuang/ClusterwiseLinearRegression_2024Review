clear
clc

% d = 0.5;
% rho = @(z) min([(1-(1-z.^2).^3),ones(size(z))],[],2);
% c = linspace(1.5,1.6,100)';

d = 0.5;
rho = @(z) min([0.5*z.^2,abs(z)-0.5],[],2);
c = linspace(0.7,0.9,100)';

% d = 1/7; %0.243600636; %0.75;
% rho = @(z) log(1+z.^2);
% c = linspace(1,3,100)'; %linspace(1.5,2,100)';

% d = 0.25;
% rho = @(z) 1-exp(-z.^2);
% c = linspace(1.5,2,100)';

I   = zeros(length(c),1);
for index = 1 : length(c)
    x = linspace(-6,6,10000)';
    I(index) = 1/sqrt(2*pi)* trapz(x,rho(x/c(index)).*exp(-0.5*x.^2));
end
err = d-I;

figure(1)
subplot(2,1,1)
plot(c,err)
grid on
title('err')
subplot(2,1,2)
plot(c,I)
grid on
title('I')

% x = linspace(-4,4,10000)';
% 1/sqrt(2*pi)*trapz(x,rho(x).*exp(-0.5*x.^2))

%% Tukey

% c = 1.5476;
% c2 = c^2;
% c4 = c2^2;
% 
% % x	= 2*(rand(10000,1)-0.5);
% x	= randn(10000,1);
% x(1:500) = 4*randn(500,1);
% 
% w       = repmat(1/length(x),size(x));
% 
% s = zeros(20,1);
% s(1) = std(x);
% for index = 2 : length(s)
%     z       = x/s(index-1);
%     z2      = z.^2;
%     u       = 3-3*(z2/c2)+(z2.^2/c4);
%     bOut    = z2/c2>1;
%     u(bOut) = 1./(z2(bOut)/c2);
%     
%     s(index) = sqrt(sum(w.*u.*x.^2)/sum(w)/d/c2);
% end
% 
% figure(2)
% plot(s)
% grid on
% title(strcat('s = ',num2str(s(end))))

%% Huber

c = 0.7979;
c2= c^2;

% x	= 2*(rand(10000,1)-0.5);
x	= randn(10000,1);
x(1:500) = 4*randn(500,1);

w       = repmat(1/length(x),size(x));

s = zeros(20,1);
s(1) = std(x);
for index = 2 : length(s)
    z       = x/s(index-1);
    u       = ones(size(z));
    bOut    = z>c;
    u(bOut) = c./z(bOut);
    
    s(index) = sqrt(sum(w.*u.*x.^2)/sum(w)/d/c2);
end

figure(2)
plot(s)
grid on
title(strcat('s = ',num2str(s(end))))


%% Cauchy

% c = 2.4027;%2.335; %sqrt(3); %0.7568; %
% c2 = c^2;
% c4 = c2^2;
% 
% % x	= 2*(rand(10000,1)-0.5);
% x	= 2*randn(10000,1);
% x(1:500) = 4*randn(500,1);
% 
% u       = zeros(size(x));
% w       = repmat(1/length(x),size(x));
% 
% s = zeros(100,1);
% s(1) = std(x);
% for index = 2 : length(s)
%     z       = x/s(index-1);
%     v2      = z.^2/c2;
%     
%     bZero   = z < 0.1*c;
%     t       = v2(bZero);
%     u(bZero)= 1 - t/2 + t.^2/3 - t.^3/4 + t.^4/5;
%     t       = v2(~bZero);
%     u(~bZero)=log(1+t)./t;
%     
%     s(index) = sqrt(sum(w.*u.*x.^2)/sum(w)/d/c2);
% end
% 
% figure(2)
% plot(s)
% grid on
% title(strcat('s = ',num2str(s(end))))

%% Welsch

% c  = 1.6035;
% c2 = c^2;
% c4 = c2^2;
% 
% % x	= 2*(rand(10000,1)-0.5);
% x	= 2*randn(10000,1);
% x(1:500) = 4*randn(500,1);
% 
% u       = zeros(size(x));
% w       = repmat(1/length(x),size(x));
% 
% s = zeros(100,1);
% s(1) = std(x);
% for index = 2 : length(s)
%     z       = x/s(index-1);
%     v2      = z.^2/c2;
%     
%     bZero   = z < 0.1*c;
%     t       = v2(bZero);
%     u(bZero)= 1 - t/2 + t.^2/6 - t.^3/24 + t.^4/120;
%     t       = v2(~bZero);
%     u(~bZero)=(1-exp(-t))./t;
%     
%     s(index) = sqrt(sum(w.*u.*x.^2)/sum(w)/d/c2);
% end
% 
% figure(2)
% plot(s)
% grid on
% title(strcat('s = ',num2str(s(end))))