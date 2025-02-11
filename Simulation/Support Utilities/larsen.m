% Title: larsen
% Author: Sanush K Abeysekera
% Date: 19/08/2022
% Description: Least Angle Regression Elastic Net (LARSEN) algorithm.
% Ref: B. Efron, T. Hastie, I. Johnstone, and R. Tibshirani. Least angle 
%      regression. Annals of Statistics. 32(2):407-499, 2004.
% Ref: H. Zou and T. Hastie. Regularization and variable selection via the
%      elastic net. J. Royal Stat. Soc. B. 67(2):301-320, 2005.

function [beta,stand] = larsen(X, y, w, options)
% X     -> Predictor covariates [n,p]
% y     -> Response [n,1]
% w     -> Weight vector (0~1) [n,1]
% beta  -> Regression coefficients (0~1) [n,p]
% stand -> Standardization structure [xmean,xstd,ymean,ystd]

% Options structure
%   singularthresh: Terminate lar when Cholesky Decomposition is impossible
%             because the matrix becomes singular (inverse do not exist)
%   stdtype:  Standardization type
%             0: No standardization, unweighted
%             1: Mean 0, standard deviation 1, weighted
%   delta:    Weight on the L2 penalty on the regression coefficients.
%             delta > 0: For p>n and prompt grouping of correlated features
%             delta = 0: LASSO solution
%   stop:     Stopping criteria
%             stop > 0: The L1 upper bound on BETA coefficient
%             stop < 0: Negative amount of active variables in BETA
%             stop = 0: Obtain full OLS solution without sparsity or until
%             lars terminates because insufficient singular value defined
%             by singularthresh

%% Define Inputs
% Default options parameters
singularthresh  = 0;
stdtype         = 1;
delta           = 0;
stop            = 0;
[n,p]           = size(X);

% User defined options
if (nargin < 4 || isempty(options))
    % Default options
else
    % Check existing fields
    exist_thr = isfield(options, 'singularthresh');
    exist_std = isfield(options, 'stdtype');
    exist_del = isfield(options, 'delta');
    exist_stp = isfield(options, 'stop');
    
    if (exist_thr), singularthresh  = options.singularthresh;   end
    if (exist_std), stdtype         = options.stdtype;          end
    if (exist_del), delta           = options.delta;            end
    if (exist_stp), stop            = options.stop;             end
end

% Default weight vector
if (nargin < 3 || isempty(w))
    w = ones(n,1) .* (1/n);
end

% Standardize data
[X,y,~,stand] = lars_standardize(X,y,w,stdtype);

% Correction of stopping criterion to fit Naïve Elastic Net
if (delta > 0 && stop > 0)
    stop      = stop/(1 + delta);
end

%% LARS
% I: Inactive set
% A: Active set
% R: Cholesky factorization R'R = X'X where R is upper triangular
% mu:        Total squared error with current BETA
% Step:      Current lars step
% maxStep:   lars does not work well for n>p, so stop at p or n/2
I           = 1:p;
A           = [];
R           = [];
mu          = 0;
Step        = 1;
maxStep     = min(p,round(n/2));

% bContinue: Violates the early stopping criteria
% bSkip:     Violates the sign restriction
bContinue   = true;
bSkip       = false;

% BETA storage variable
beta        = zeros(p,maxStep);

while (~isempty(I)) && (Step <= maxStep) && (bContinue)
    % Residue vector
    r    = y - mu;
    
    % Find maximum correlation
    c    = X'*r;
    [cmax,cidxI] = max(abs(c(I)));
    cidx = I(cidxI);    % Index of next active variable
    
    if (~bSkip)
        % Add variable to BETA
        % Cholesky Factorization: add variable
        R     = lars_cholinsert(R, X(:,cidx), X(:,A), delta);

        % Update active coefficients
        A     = [A, cidx];          % Add to active set
        I(cidxI)   = [];            % Drop from inactive set
    else
        % Sign restriction violates
        % If a variable has been dropped, do one step with this
        % configuration (don't add new one right away)
        % Continue with appended R without changing the active and inactive
        % sets
        bSkip = false;
    end
    
    if rcond(R) > singularthresh
        % OLS BETA solution
        b_OLS = R\(R'\(X(:,A)'*y)); % same as X(:,A)\y, but faster
        d     = X(:,A)*b_OLS - mu;
        
        % Compute length of walk along equiangular direction
        if isempty(I)
            % If all variables are active go to the OLS solution
            gamma     = 1;
        else
            % Calculate length of walk
            cd        = X'*d;
            gamma_hat = [(c(I) - cmax)./(cd(I) - cmax) ;... 
                (c(I) + cmax)./(cd(I) + cmax)];
            
            gamma_hat = sort(gamma_hat(gamma_hat>0));
            % Faster than min(gamma_hat(gamma_hat>0)) 
            
            if isempty(gamma_hat)
                bContinue = false;
                gamma     = 0;
            else
                gamma     = gamma_hat(1);
            end
        end
        
        % Compute smallest length of walk along equiangular direction until
        % sign change
        if (delta == 0) % For lasso
            gamma_tilde = beta(A(1:end-1),Step)./...
                (beta(A(1:end-1),Step) - b_OLS(1:end-1,1));
            gamma_tilde(gamma_tilde <= 0) = inf;
            [gamma_tilde, dropIdx]        = min(gamma_tilde);
            
            % Check if variable should be dropped
            if gamma_tilde < gamma
                bSkip     = true;
                gamma     = gamma_tilde;
            end
        end
        
        % Check if BETA storage variable must grow
        if (size(beta,2) < Step+1)
            beta  = [beta, zeros(p, size(beta,2))];
        end
        
        % Update BETA
        beta(A,Step+1) = beta(A,Step) + gamma*(b_OLS - beta(A,Step));
        
        % Update position
        mu    = mu + gamma*d;
        
    else
        % No propergatable path available
        bContinue  = false;
    end
    
    % If LASSO condition not satisfied, drop variable from active set
    if (bSkip)
        I     = [I A(dropIdx)]; % Add drop variable to inactive set
        A(dropIdx) = [];        % Remove drop variable from active set
        R     = lars_choldelete(R, dropIdx);
    end
    
    % Early stopping criteria: number of variables
    if stop < 0
        bContinue  = (length(A) < -stop);
    end
    
    % Early stopping criteria: bound on L1 norm
    if stop > 0
        t2 = sum(abs(beta(:,Step+1)));
        if (t2 >= stop)
            t1 = sum(abs(beta(:,Step)));
            s  = (stop - t1)/(t2 - t1); % Interpolation factor 0<s<1

            beta(:,Step+1) = beta(:,Step) +... 
                s*(beta(:,Step+1) - beta(:,Step));
            bContinue = false;
        end
    end
    
    % Increment step counter
    Step         = Step + 1;
end

%% Adjust BETA
% Trim un-used BETA
if (size(beta,2) > Step)
    beta(:, (Step+1):end) = [];
end

% Un-normalized BETA
xmean       = stand.xmean;
xstd        = stand.xstd;
ymean       = stand.ymean;
ystd        = stand.ystd;

% Adjust to avoid double shrinkage (Non-Naïve Elastic Net)
beta        = (1 + delta) * beta;

% Remove standardization from BETA
beta        = [ymean - xmean*(beta*ystd./xstd') ; beta*ystd./xstd'];

end
