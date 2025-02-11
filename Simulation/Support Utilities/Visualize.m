function Visualize(Solution,X,y,N,beta,p)
% Visualisation of Solution
% 25/11/2022
%%

delta = sum((beta(:,1)-Solution.beta).^2,1);
b1    = delta(1) < delta(2);

figure(1)
stem(beta(:,1),'b.'); hold on
if b1
    plot(1:p,Solution.beta(:,1),'bo'); hold off
else
    plot(1:p,Solution.beta(:,2),'bo'); hold off
end
box off

figure(2)
stem(beta(:,2),'r.'); hold on
if b1
    plot(1:p,Solution.beta(:,2),'ro'); hold off
else
    plot(1:p,Solution.beta(:,1),'ro'); hold off
end
box off

figure(3)
if b1
    yhat(1:round(N/2),1) = X(1:round(N/2),:)*Solution.beta(:,1) + Solution.Offset(1);
    yhat(1+round(N/2):N,1) = X(1+round(N/2):N,:)*Solution.beta(:,2) + Solution.Offset(2);
    plot(y(1:round(N/2)),yhat(1:round(N/2)),'b.',y(1+round(N/2):N),yhat(1+round(N/2):N),'r.')
else
    yhat(1:round(N/2),1) = X(1:round(N/2),:)*Solution.beta(:,2) + Solution.Offset(2);
    yhat(1+round(N/2):N,1) = X(1+round(N/2):N,:)*Solution.beta(:,1) + Solution.Offset(1);
    plot(y(1:round(N/2)),yhat(1:round(N/2)),'r.',y(1+round(N/2):N),yhat(1+round(N/2):N),'b.')
end
box off
grid on
end