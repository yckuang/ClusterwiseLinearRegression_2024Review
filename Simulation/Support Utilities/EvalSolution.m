function [yhat,varargout] = EvalSolution(Solution,X,varargin)

bErr = false;
if nargin > 2
    y = varargin{1};
    bErr = true;
end

yhat    = X*Solution.beta;

if bErr 
    err	= yhat - y;
    varargout{1} = err;
end

