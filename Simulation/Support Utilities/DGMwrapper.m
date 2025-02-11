function Solution = DGMwrapper(b,y,X1,OptDGM)
evalcount = 0;

[N,pp1] = size(X1);
numelb	= numel(b);
if mod(numelb,pp1) ~= 0
    error('incompatible size between X and beta')
end

k = numelb/(pp1);
b = b(:)';
w = zeros(N,k);


lambda = OptDGM.lambda;

fun = @(bb) DGM_EvalCost(bb,y,X1,pp1,k); 
[b,info] = DiscreteGradientOptimization(fun,b,lambda,OptDGM);
evalcount	= evalcount	+ info.evalcount;

b           = reshape(b,[pp1,k]);
err         = y-X1*b;                                                       evalcount	= evalcount	+ 1;
[~, ptk]    = min(abs(err),[],2);
w(sub2ind(size(w),[1:N]',ptk)) = 1;
pr          = mean(w./sum(w,2),1);
r2          = 1 - sum(err.^2)/sum((y-mean(y)).^2);

Solution.beta   = b;
Solution.w      = w;
Solution.pr     = pr;
Solution.r2     = r2;
Solution.iterations = [evalcount,0];



