% Robust Tensor Power Method 
% Anandkumar 2014 Tensor Decomposition for Learning Latent Variable Models
%  - orthogonal decomposable symmetrical third order tensor 

function [theta,eigenval] = TPM(Tensor,varargin)

[K,~,~] = size(Tensor);

NumTrial    = 10*K;
MaxIter     = 100;
RelTol      = 1e-4;
MaxEmptyLoop= 10;

if nargin > 1
    for index = 1 : 2 : length(varargin)
        switch lower(varargin{index})
            case 'numtrial'
                NumTrial    = varargin{index+1};
            case 'maxiter'
                MaxIter     = varargin{index+1};
            case 'reltol'
                RelTol      = varargin{index+1};
            case 'maxemptyloop'
                MaxEmptyLoop= varargin{index+1};
        end
    end
end

%%
theta   = zeros(K,K);
eigenval= zeros(1,K);
for k = 1 : K
    
    zeroupdatecount = 0;
    for trial = 1 : NumTrial
        
        beta = randn(K,1);
        beta = beta/norm(beta);
        
        bContinue = true;	iter = 0;
        while bContinue
            iter = iter + 1;
            
            TensorApply = sum( Tensor .* reshape(beta,[1,1,K]), 3);
            TensorApply = sum( TensorApply .* beta', 2);
            betanew     = TensorApply / norm(TensorApply);
            
            bConverge = norm(beta-betanew) < RelTol;    % norm(beta)==1
            bContinue = ~bConverge && iter < MaxIter;
            
            beta = betanew;
        end
        TensorApply = sum( Tensor .* reshape(beta,[1,1,K]), 3);
        TensorApply = sum( TensorApply .* beta', 2);
        newscore = TensorApply'*beta;
        
        if newscore > eigenval(k)
            theta(:,k)  = beta;
            eigenval(k)	= newscore;
            zeroupdatecount = 0; % Reset zero update
        else
            zeroupdatecount = zeroupdatecount + 1;
        end        
        if zeroupdatecount > MaxEmptyLoop
            break;  % Early termination when there is no more new solutions
        end
        
    end
    
    % Deflate tensor for the next iteration
    EigenTensor = (theta(:,k)*theta(:,k)').*reshape(theta(:,k),[1,1,K]);
    Tensor = Tensor - eigenval(k)*EigenTensor;
    
end