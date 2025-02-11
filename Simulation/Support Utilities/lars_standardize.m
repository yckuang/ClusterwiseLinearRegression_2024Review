% Title: standardize
% Author: Sanush K Abeysekera
% Date: 19/08/2022
% Description: Standardize data for lars.
% Ref: B. Efron, T. Hastie, I. Johnstone, and R. Tibshirani. Least angle 
%      regression. Annals of Statistics. 32(2):407-499, 2004.
% Ref: H. Zou and T. Hastie. Regularization and variable selection via the
%      elastic net. J. Royal Stat. Soc. B. 67(2):301-320, 2005.

% Notes: Weighted OLS using normalized X and y
%        beta = inv(X'*diag(w)*X)*X'*diag(w)*y

function [X,y,w,stand] = lars_standardize(X,y,w,stdtype)
% X     -> Predictor covariates [n,p]
% y     -> Response [n,1]
% w     -> Weight vector (0~1) [n,1]
% stdtype  -> Standardization type
%          0: No Standardization, unweighted
%          1: Mean 0, standard deviation 1, weighted

% X     -> Standardized predictor covariates [n,p]
% y     -> Standardized response [n,1]
% w     -> Standardized weight (0~1) [n,1]
% stand -> Standardization structure [xmean,xstd,ymean,ystd]

% Dataset dimensions
[n,p]   = size(X);

% Default standardization type
if (nargin < 4 || isempty(stdtype))
    stdtype = 1;
end

% Default weight vectors
if (nargin < 3 || isempty(w))
    w  = ones(n,1)./n;
else
    w  = w./sum(w);
end

% Normalization type
switch stdtype
    case 0  % No normalization
        stand.xmean = zeros(1,p);
        stand.xstd  = ones(1,p);
        stand.ymean = 0;
        stand.ystd  = 1;

    case 1  % Standard normalization
        % Weighted data
        stand.xmean = mean(X,1);
        stand.xstd  = sqrt(mean((X-stand.xmean).^2,1));
        stand.ymean = mean(y,1);
        stand.ystd  = sqrt(mean((y-stand.ymean).^2,1));

        % Normalized data
        X = (X - stand.xmean) ./ stand.xstd;
        y = (y - stand.ymean) ./ stand.ystd;
        
        % Weighting
        X = sqrt(w) .* X;
        y = sqrt(w) .* y;
end

end