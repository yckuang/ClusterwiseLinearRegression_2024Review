function beta   = GenBeta(p, K, rho, nonzerodimratio)

if p < K
    error('This code will only run if the input dimension is equal or greater than than the number of clusters')
end

numzero = round(p*(1-nonzerodimratio));
numfill = p-numzero;

%% Generate beta with the given rho

b   = randn(numfill,1);
b   = b./norm(b);

nullspc	= null(b');
ns      = zeros(size(nullspc));
for index = 1 : numfill-1
    ns(:,index)	= nullspc*(2*rand(size(nullspc,2),1)-1);
    ns(:,index)	= ns(:,index)/norm(ns(:,index));
    nullspc = null([b,ns]'); 
end
nullspc = ns;

for k = 2 : K
    v1      = nullspc(:,k-1);
    
    G       = b'*b;
    invG    = inv(G);
    v2      = b*sum(invG,2);
    v2      = v2/norm(v2);
    
    alpha     = rho/mean(b'*v2);
    if alpha >= 1
        alpha = 0.99;
        warning('can''t generate the required correlation vector')
    end
    b(:,k) = alpha*v2 + sqrt(1-alpha^2)*v1;
    b(:,k) = b(:,k)/norm(b(:,k));
end

%%
beta = zeros(p,K);
for index = 1 : size(b,2)
    ptfill = randperm(p,numfill);
    beta(ptfill,index) = b(:,index);
end
