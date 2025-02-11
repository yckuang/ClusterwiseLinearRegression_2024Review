function Solution = Tensor_Sedghi2016(y0,X0,Option)

global N p K p2 p3 pK
[N,p]   = size(X0);
p2      = p*p;
p3      = p*p2;

K       = Option.NumMix;
pK      = p*K;

%%
if isfield(Option,'RegressWeight')              % Outlier weight-out method to achieve "robust" regression
    RegressWeight	= Option.RegressWeight;
else
    RegressWeight.Name= '';
end
if isfield(Option,'FactorizationMethod')            % Tensor factorization method
    FactorizationMethod = Option.FactorizationMethod;
else
    FactorizationMethod.Name= 'rpsf';
end

RegressWeight.Name      = lower(RegressWeight.Name);
FactorizationMethod.Name= lower(FactorizationMethod.Name);

%% Normalize X0 and y0 (X0 and y0 are assumed to have been centered (zero mean) )

sdy     = std(y0);
y       = y0/sdy;
y2      = y.^2;
y3      = y.*y2;

sdX     = std(X0,[],1);
X       = zeros(N,p);
bValid  = sdX > eps;
X(:,bValid) = X0(:,bValid)./sdX(bValid);

% Align to the priciple axes (removes predictor correlation)
[~,~,R] = svd(X);   % R = eye(p);
x       = X*R;


y3xavg	= mean(y3.*x,1)';

%% Form third order discriminative tensor

M2 = zeros(p,p);
for index = 1 : N
   v = x(index,:);
   M2 = M2 + y2(index)*(v'*v-eye(p,p));
end
M2 = M2/N;
[U,S,~] = svd(M2);
S   = diag(S);
W   = U(:,1:K)*diag( 1./sqrt(S(1:K)) );

M3 = zeros(p2,p);
for index = 1 : N
   v = x(index,:);
   M3 = M3 + reshape((y3(index)*v)'*v,[p2,1])*v;
end
M3 = reshape(M3/N,[p,p,p]);
for index = 1 : p
    M3(:,index,index) = M3(:,index,index)-y3xavg;
    M3(index,:,index) = M3(index,:,index)-y3xavg';
    M3(index,index,:) = M3(index,index,:)-reshape(y3xavg,[1,1,p]);
end

M3WWW = zeros(K,K,K);
for k = 1 : K
    TensorMap_k	= sum(M3 .* reshape(W(:,k),[1,1,p]),3);
    for j = 1 : K
        TensorMap_j	= sum(TensorMap_k .* W(:,j)', 2);
        M3WWW(:,j,k)= W' * TensorMap_j;
    end
end


%% Tensor factorization

switch FactorizationMethod.Name
    case 'tpm'
        [eigenv,eigenval]   = TPM(M3WWW);
    case 'rpsf'
        [eigenv,eigenval]   = RPSF(M3WWW);
end
v           = (R*pinv(W')*(eigenv));
v(bValid)   = (v(bValid)./sdX(bValid));%*sdy;
betahat     = v./sqrt(sum(v.^2,1));

% Estimate the magnitude of beta through constrained optimization
A           = (sdX.^2)*(betahat.^2);
B           = sdy^2;
prini       = repmat(1/K,[K,1]);
magini      = sqrt(A\(B/prini(1)));
param       = [prini;magini];
Aeq         = [ones(1,K),zeros(1,K)];
Beq         = 1;
LB          = [zeros(K,1);-inf(K,1)];
UB          = [ones(K,1);inf(K,1)];

fun         = @(param) EstimateNorm(param,y,X*betahat);
nonlincon   = @(param) NonlinCon(param,A,B);
param       = fmincon(fun,param,[],[],Aeq,Beq,LB,UB,nonlincon);

Solution.pr     = param(1:K)';
Solution.beta   = param(K+1:end)'.*betahat;


end

%% Private function
function errsq = EstimateNorm(param,y,Xb)
    global N K

    pr          = param(1:K)';
    mag         = param(K+1:end)';
    errmat      = y - mag.*Xb;
    [err,pt]	= min( abs(errmat), [],2);
    
    premp       = sum(pt==[1:K],1)/N;
        
    errsq       = sum(err.^2) + mean((pr-premp).^2)*N;

end

function [c,ceq] = NonlinCon(param,A,B)
    global K
    pr      = param(1:K);
    mag2    = param(K+1:end).^2;
    
    ceq     = B - A*(pr.*mag2);
    c       = [];
end
