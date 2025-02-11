function Solution = CorrectPermutation(Solution, beta)
% Corrects the position of coloumns of the calculated coefficient matrix to
% match that of the true coefficient matrix
% 25/11/2022

betahat     = zeros(size(beta));
what        = zeros(size(Solution.w));
K           = size(betahat,2);

if K>1
    idx     = zeros(K^2,2);
    dist    = zeros(K^2,1);
    for i = 1 : K
        ptnow       = (i-1)*K+1 : i*K;
        dist(ptnow)	= sum(abs(beta(:,i)-Solution.beta),1);
        idx(ptnow,1)= i;
        idx(ptnow,2)= 1:K;
    end
    
    pair = zeros(K,2);
    for i = 1 : K-1
        [~,ptmin]       = min(dist);
        pair(i,:)       = idx(ptmin,:);
        bRemove         = idx(:,1)==pair(i,1) | idx(:,2)==pair(i,2);
        dist(bRemove)   = [];
        idx(bRemove,:)  = [];
        betahat(:,pair(i,1))= Solution.beta(:,pair(i,2));
        what(:,pair(i,1))   = Solution.w(:,pair(i,2));
    end
    betahat(:,idx(1,1)) = Solution.beta(:,idx(1,2));
    what(:,idx(1,1))    = Solution.w(:,idx(1,2));
    pr = mean(what,1);
    
    Solution.beta   = betahat;
    Solution.w      = what;
    Solution.pr     = pr;
end