% Random Projection Simultaneous Factorization
% Kuleshov 2015 Tensor Factorization via Matrix Factorization
%  - orthogonal decomposable symmetrical third order tensor 

function [eigenv,eigenval] = RPSF(Tensor,varargin)

[K,~,~] = size(Tensor);

NumTrial = 10*K;

if nargin > 1
    for index = 1 : 2 : length(varargin)
        switch lower(varargin{index})
            case 'NumTrial'
                NumTrial   = varargin{index+1};
        end
    end
end

%%

% STAGE 1: Random projections

u   = randn(K,NumTrial);
u   = u./sqrt(sum(u.^2,1));

M0 = zeros(K,K*NumTrial);
for index = 1 : NumTrial
    M0(:,(index-1)*K+1:index*K) = sum(Tensor.*reshape(u(:,index),[1,1,K]),3);
end
[Q, ~, ~] = jacobi(M0, 1e-4, K);
u = Q(:,1:K);

% STAGE 2: Plugin projections
M0 = zeros(K,K*K);
for index = 1 : K
    M0(:,(index-1)*K+1:index*K) = sum(Tensor.*reshape(u(:,index),[1,1,K]),3);
end
[Q, D, ~] = jacobi(M0, 1e-4, K);
eigenv = Q(:,1:K);

eigenval = zeros(1,K);
for k = 1 : K
    d       = diag(D(:,(k-1)*K+1:k*K));
    [maxd, pt] = max(abs(d));
    if d(pt) < 0
        eigenv(:,pt) = -eigenv(:,pt);
    end    
    eigenval(k) = maxd;
end


end

%% Private functions

