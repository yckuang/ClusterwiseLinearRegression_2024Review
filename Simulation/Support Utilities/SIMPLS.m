% This implementation assumes that X and Y have zero mean

function [beta,varargout] = SIMPLS(X,Y,w,options)

% Maximum number of latent components
if isfield(options,'NumComp') 
    NumComp = options.NumComp;
else
    NumComp = min(size(X)-1);
end

% Threshold of variance explained in X. Terminate expansion when the
% threshold is reached
if isfield(options,'XThreshold')    
    thX = options.XThreshold;
else
    thX = 1;
end

% Threshold of variance explained in Y. Terminate expansion when the
% threshold is reached
if isfield(options,'YThreshold')    
    thY = options.YThreshold;
else
    thY = 1;
end

[n,p] = size(X);
[~,m] = size(Y);

%% Main iteration loop
if isempty(w)
    w = repmat(1/n,n,1);
else
    if length(w) == n
        w = abs(w);
        w = w/sum(w);
    else
        error('The weight vector must be of the same length as the observations. Supply empty matrix [] for unweighted data')
    end
end
sqrtw = sqrt(w);
X = sqrtw.*X;
Y = sqrtw.*Y;
S = (X'*Y);

ssXRef = ss(X);
ssYRef = ss(Y);

R = zeros(p,NumComp);
T = zeros(n,NumComp);
P = zeros(p,NumComp);
Q = zeros(m,NumComp);
U = zeros(n,NumComp);
V = zeros(p,NumComp);

bContinue = true;       c = 0;
while bContinue
    c = c + 1;
    
    [eigV,~] = eig(S'*S);
    Q(:,c) = eigV(:,end);
    
    r = S*Q(:,c);
    R(:,c) = r/norm(r);
    
    t = X*R(:,c);
    T(:,c) = t/norm(t);
    
    P(:,c) = X'*T(:,c);
    Q(:,c) = Y'*T(:,c);
    
    u = Y*Q(:,c);
    u = u - T*(T'*u);
    U(:,c) = u/norm(u);
    
    v = P(:,c);
    v = v - V*(V'*v);
    V(:,c) = v/norm(v);

    S = S - V(:,c)*(V(:,c)'*S);
    
    varX = ss(P) / ssXRef;      % percentage X-variance explained
    varY = ss(Q) / ssYRef;      % percentage Y-variance explained
    
    bContinue = (c < NumComp) & (varX < thX) & (varY < thY);
end
if c < NumComp
    R(:,c+1:end) = [];
    T(:,c+1:end) = [];
    P(:,c+1:end) = [];
    Q(:,c+1:end) = [];
    U(:,c+1:end) = [];
    V(:,c+1:end) = [];
end

beta = R*Q';
% beta = beta/norm(beta);

if nargout > 1
    OutVar.R	= R;
    OutVar.T	= T;
    OutVar.P	= P;
    OutVar.Q	= Q;
    OutVar.U	= U;
    OutVar.V	= V;
    
    varargout{1} = OutVar;
end

end


function T = ss(X)  % Compute sumsquare

T = sum(sum(X.*X));

end
