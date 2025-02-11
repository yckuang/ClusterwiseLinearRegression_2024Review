function Solution = IncrDGM(X,y,Solution,Option)
% Bagirov 2013, Nonsmooth nonconvex optimization approach to clusterwise
% linear regression problems, European Journal of Operational Research


[N,p]   = size(X);
K       = Option.NumMix;

if ~isfield(Option,'Incr') || ~isfield(Option.Incr,'g1')
    if N <= 200
        Option.Incr.g1	= 0.30;
    elseif N <= 1000
        Option.Incr.g1   = 0.50;
    else
        Option.Incr.g1   = 0.95;
    end
end
if ~isfield(Option,'Incr') || ~isfield(Option.Incr,'g21')    
    Option.Incr.g21  = 0.95;
end
if ~isfield(Option,'Incr') || ~isfield(Option.Incr,'g22')    
    Option.Incr.g22  = 10;
end
if ~isfield(Option,'Incr') || ~isfield(Option.Incr,'g3')    
    Option.Incr.g3   = 10;
end
if ~isfield(Option,'Incr') || ~isfield(Option.Incr,'epsilon')
    Option.Incr.epsilon   = 1e-4;
end
if ~isfield(Option,'Incr') || ~isfield(Option.Incr,'maxiter')
    Option.Incr.maxiter   = 1000;
end

EvalCount = 0;
RegrCount = 0;

% Find the first component
X1      = [ones(N,1),X];
% beta	= X1\y .* (1 + 0.2*randn(p+1,1));                                   RegrCount = RegrCount+1;
beta	= X1\y;                                                             RegrCount = RegrCount+1;

Solution.w	= ones(N,1);

% Add components incrementally
for k = 2 : K
    Option.NumMix   = k;
    [A3,aerr,EvalCount,RegrCount]   = Algorithm3(beta,X1,y,Option,EvalCount,RegrCount);
    [A5,EvalCount,RegrCount]        = Algorithm4(A3,aerr,X1,y,Option,EvalCount,RegrCount);
    
    NumTrial    = size(A5,2);
    f           = inf;
    ptBest      = 1;
    for index = 1 : NumTrial
        b               = [beta,A5(:,index)];
        % DGM optimization
        Solution        = DGMwrapper(b,y,X1,Option.DGM);
        
        % Retain best solution
        r               = min(abs(y - X1*Solution.beta),[],2);
        newf            = rms(r);
        if newf < f
            f       = newf;
            ptBest  = index;
        end
        EvalCount = EvalCount + Solution.iterations(1);
        RegrCount = RegrCount + Solution.iterations(2);
    end
    beta   = [beta,A5(:,ptBest)]; %#ok<AGROW>
end
err             = y - X1*beta;                                         EvalCount = EvalCount + 1;
[~,ptk]         = min(abs(err),[],2);
w               = zeros(N,k);
w(sub2ind(size(w),[1:N]',ptk)) = 1;
% Run EM optimization
Solution.beta   = beta;
Solution.w      = w;
Solution.pr     = mean(w./sum(w,2),1);

Solution.iterations = [EvalCount,RegrCount];

end


%%
% Algorithm 3
function [A3,aerr,EvalCount,RegrCount] = Algorithm3(b,X1,y,Option,EvalCount,RegrCount)
    
    g1  = Option.Incr.g1;
    g21	= Option.Incr.g21;
    g22 = Option.Incr.g22;

    [N,pp1]     = size(X1);
    [~,km1]     = size(b);
    p           = pp1-1;
    
    X           = X1(:,2:end);
    u           = X*b(2:end,:);                                   EvalCount = EvalCount + 1;
    err         = (y-b(1,:))-u;
    [aerr,ptk]  = min(abs(err),[],2);
    u           = u(sub2ind(size(u),[1:N]',ptk)); % f_k-1

    % set A1 equ (8)
    bb      = y-u;
    drop    = zeros(1,N);
    parfor index = 1 : N    % use for loop to avoid creating N^2 matrix
        aerr1       = min(abs(y-X*b(2:end,ptk(index))-bb(index)),[],2); %abs(bb-bb(index));
        d           = aerr-aerr1;
        drop(index) = sum(d(d>0));
    end
    maxdrop = max(drop);
    A1   = find(drop > g1*maxdrop);
    
%     % (Additional filter not in the original paper). Truncate A1 by removing similar pointsets. 
%     % This will reduce the number of regression operations in dense dataset
%     A1          = [];
%     d           = aerr-abs(bb-bb(nList)');
%     bSel        = d > 0;
%     TotalSel    = sum(bSel,1);
%     covmat      = (bSel'*bSel)./sqrt(TotalSel'*TotalSel);
%     bHighCorr   = covmat > g21;             % Filter off any suggestion that are highly correlated
%     while ~isempty(nList)                   % Only retain the proposal with the largest drop in every group of highly correlated point set
%         pta1	= find(bHighCorr(1,:));
%         a1      = nList(pta1);
%         [~,jj]  = max(drop(a1));
%         A1 = [A1,a1(jj)];        
%         nList(pta1)         = [];
%         bHighCorr(pta1,:)	= [];
%         bHighCorr(:,pta1)	= [];
%     end

    % set A2
    A2  = zeros(p+1,length(A1));
    f   = zeros(1,length(A1));                                              RegrCount = RegrCount + length(A1);
    bbA1    = bb(A1);
    aerrA1  = aerr(A1);
    for index = 1 : length(A1)
        aerr1   = abs(bbA1-bbA1(index));
        d       = aerrA1-aerr1;
        % set B equ (9)
        bSel    = d > 0;
        A2(:,index)	= X1(A1(bSel),:)\y(A1(bSel));
        err2        = y-X1*A2(:,index);
        aerr2       = abs(err2);
        f(index)    = sum(min([aerr,aerr2],[],2));
    end
    minf = min(f);
    
    % set A3 equ (12)
    A3 = A2(:,f < g22*minf);

end

%% Algorithm 4, EM refinement of A3
function [A5,EvalCount,RegrCount] = Algorithm4(A3,aerr,X1,y,Option,EvalCount,RegrCount)
    
    g3      = Option.Incr.g3;
    epsilon = Option.Incr.epsilon;
    maxiter = Option.Incr.maxiter;
    
    NumTrial = size(A3,2);
    
    f = zeros(1,NumTrial);
    for index = 1 : NumTrial
        b = A3(:,index);
        rmsb = rms(b);
        loopcounter = 0;
        bContinue = true;
        while bContinue
            loopcounter = loopcounter + 1;
            aerr1 = abs(y-X1*b);                                            EvalCount = EvalCount + 1;
            bSel = aerr1 < aerr;
            bb	= X1(bSel,:)\y(bSel);                                       RegrCount = RegrCount + 1;
            if rms(bb-b)/rmsb < epsilon || loopcounter > maxiter
                bContinue = false;
            else
                b = bb;
            end
        end
        A3(:,index) = b;
        
        % Remove redundant b (converging to the same local minima)
        diff = sqrt(sum(((b-A3(:,1:index-1)).^2),1))/norm(b);
        if any(diff < epsilon)
            f(index) = inf;
        else            
            f(index) = sum(min([aerr,aerr1],[],2));
            A3(:,index) = b;
        end
    end
    minf = min(f);
    
    A5 = A3(:,f <= g3*minf);
    
end