function [X, y, beta, varargout]   = GenData(Input)
% Generates data to be processed and plotted.
% 21/11/2022
%%

p               = Input.p;
K               = Input.K;
N               = Input.N;
nonzerodimratio = Input.nonzerodimratio;
dotprod         = Input.dotprod;
outlier         = Input.outlier;
noise           = Input.noise;
covariance      = Input.covariance;
offset          = Input.offset;
offsety         = Input.offsety;

NumCorupt = round(outlier*N);
blocksize = round((N-NumCorupt)/K);

X = zeros(N,p);
y = zeros(N,1);

beta0   = zeros(1,K);

% Assume that each component of X has unity variance
Sigma = repmat(covariance,[p,p]) + (1-covariance)*eye(p);
[eigvec,eigval] = eig(Sigma);       eigval = diag(eigval);
A = (eigvec*diag(sqrt(eigval)));    % Coordinate transformation to achieve the required covariance
if ~all(eigval > length(eigval)*eps(max(eigval)))
    error('The covariance value is not valid. Try a value between 0 and 1')
end
beta = GenBeta(p, K, dotprod, nonzerodimratio);      % p-dimensional mapping vector (without DC component)
% disp('! Not using GenBeta !')
% beta = randn(p,K);

% Generate offset mean vectors
u       = randn(p,1);
u       = u/norm(u)*offset;
uset    = [zeros(p,1),u,zeros(p,K-2)];
% Maximal tension model for unit variance and zero correlation
if K > p+1  % Strict solution doesn't exist, use optimisation to find approximate configuration
    optimopt = optimoptions(@fmincon,'Algorithm','sqp','OptimalityTolerance',1e-6,'Display','off');
    for k = 3 : K
        u0	= randn(p,1);
        u0	= u0/norm(u0)*5*offset;
        objfun      = @(u) sum(sum((u-uset(:,1:k-1)).^2));
        nonlcon     = @(u) offsetnonlcon(u,uset(:,1:k-1),offset);
        uset(:,k)   = fmincon(objfun,u0,[],[],[],[],[],[],nonlcon,optimopt);
    end
else        
    d2	 = sum((uset(:,2)).^2);     % should be sum((uset(:,2)-uset(:,1)).^2) but uset(:,1)==0
    for k = 3 : K
        centroid    = mean(uset(:,1:k-1),2);
        nsp         = null(uset(:,1:k-1)');
        b           = randn(size(nsp,2),1);  b = b/norm(b);     % random linear combination of the null space bases
        n           = nsp*b;
        n_length	= sqrt(d2 - sum(centroid.^2));
        uset(:,k)   = centroid + n_length*n;
    end
end
uset = uset - mean(uset,2);

shift = NumCorupt;
for k = 1:K
    GroupSz = min([blocksize,N-shift]);
    x = uset(:,k)' + randn(GroupSz,p);
    X(shift+1:shift+GroupSz,:) = x*A';      % Transformation to correlated predictors

    beta0(k)    = offsety * randn();        % Group dependent DC mapping noise
    y(shift+1:shift+GroupSz) = X(shift+1:shift+GroupSz,:)*beta(:,k) + beta0(k);
    %beta0(k)    = beta0(k) + (A*uset(:,k))'*beta(:,k);            % Recalculate overall DC offset
    % Inject gaussian mapping noise
    y(shift+1:shift+GroupSz) = y(shift+1:shift+GroupSz) + noise*std(y(shift+1:shift+GroupSz))*randn(GroupSz,1);
    
    shift = shift + blocksize;
end
y(1:NumCorupt) = std(y(NumCorupt+1:end))*randn(NumCorupt,1);

beta = [beta0;beta];

if nargout > 3
    % Estimate scale parameters
    shift = NumCorupt;
    for k = 1 : K
        GroupSz = min([blocksize,N-shift]);
        X1      = [ones(GroupSz,1),X(shift+1:shift+GroupSz,:)];
        yhat    = X1*beta(:,k);
        err     = yhat - y(shift+1:shift+GroupSz);
        s(k)	= std(err);
        shift = shift + blocksize;
    end
    
    yhat	= [ones(N-NumCorupt,1),X(NumCorupt+1:end,:)]*beta;
    yhat1   = yhat./s;
    yhat2   = yhat1./s;
    a2      = 1/sum(1./(s.^2));
    a       = sqrt(a2);
    b       = a2*sum(yhat2,2);
    b2      = b.^2;
    c2      = a2*sum(yhat1.^2,2);
    G       = sqrt(K)*a/prod(s.^(1/K))*exp( -0.5/a2*(c2-b2) );

    Q       = mean(G);
    R       = 1-Q;
    
    varargout{1} = R;
    
end
