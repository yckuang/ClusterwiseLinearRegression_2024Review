function Solution = EM(X,y,Solution,Option)

%% Center the data
[N,p]   = size(X);
            
X1  = [ones(N,1),X];
SSy = sum((y-mean(y)).^2);

%%
MaxIter                 = Option.MaxIter;
K                       = Option.NumMix;
prctScoreTh             = Option.prctScoreTh;
EliteLen                = Option.EliteLen;
Momentum                = Option.Momentum;

RegressWeight.Name  = '';
ScaleEstWeight.Name = '';
if isfield(Option,'Initialization')         % Initialization method
    Initialization	= Option.Initialization;
else
    Initialization.Name	= 'random';
end
if isfield(Option,'bCEM')                   % Crisp cluster membership
    bCEM	= Option.bCEM;
else
    bCEM	= true;
end
if isfield(Option,'DistType')               % distribution function that generates likelihood 
    DistType	= Option.DistType;
    if strcmpi(DistType.Name,'t')
        nu = DistType.nu;
        % Default weighting for t-robust modelling. May be overiden by
        % explicit spcification by ScaleEstWeight and RegressWeight
        ScaleEstWeight.Name= 't-robust';
        RegressWeight.Name = 't-robust';        
    end
else
    DistType.Name= 'normal';
end
if isfield(Option,'RegressWeight')           % Outlier weight-out method to achieve "robust" regression
    RegressWeight	= Option.RegressWeight;
end
if isfield(Option,'ScaleEstWeight')           % Outlier weight-out method to achieve "robust" regression
    ScaleEstWeight	= Option.ScaleEstWeight;
end
if isfield(Option,'RegressMethod')          % Regression method
    RegressMethod   = Option.RegressMethod;
else
    RegressMethod.Name= 'ls';
end

Initialization.Name     = lower(Initialization.Name);
DistType.Name           = lower(DistType.Name);
RegressWeight.Name      = lower(RegressWeight.Name);
ScaleEstWeight.Name     = lower(ScaleEstWeight.Name);
RegressMethod.Name      = lower(RegressMethod.Name);

%%
EvalCount = 0;
RegrCount = 0;

% Initialize solution
if isempty(Solution)
    switch Initialization.Name
        case {'kmean','kmean1'}    % Kmean clustering
            if strcmpi(Option.Initialization.Name,'kmean1')
                ptGroup = kmeans([y,X],K, 'MaxIter', 1000);
            else
                ptGroup = kmeans(X,K, 'MaxIter', 1000);
            end
            
            beta    = zeros(p+1,K);
            w       = ones(N,1);
            for k = 1 : K
                bSelect = ptGroup == k;
                x = X1(bSelect,:);
                beta(:,k) = x\y(bSelect);                                   RegrCount = RegrCount + 1;
            end
            Solution.beta   = beta;
            
            [~,err]	= EvalSolution(Solution,X1,y);                          EvalCount = EvalCount + 1;
            errsq   = err.^2;
            sigmasq	= sum(errsq,1)./sum(w,1);
            phi     = exp(-(err.^2)./(2*sigmasq));
            
            w       = phi./sum(phi,2);
            pr      = mean(w,1);
            pr      = pr./sum(pr);
            
            Solution.w      = w;
            Solution.pr     = pr;
            Solution.r2     = 0;
            Solution.Score  = 0;
            
        case {'tensor-chaganty2013','tensor-sedghi2016'} % Third order tensor factorization
            ymu     = mean(y);
            Xmu     = mean(X,1);
            v       = y-ymu;
            U       = X-Xmu;
            if strcmpi(Option.Initialization.FactorizationMethod, 'tensor-chaganty2013')
                Solution	= Tensor_Chaganty2013(v,U,Option);              % Tensor factorization is way more complex than what can be captured by EvalCount and RegrCount, omit
            elseif strcmpi(Option.Initialization.FactorizationMethod, 'tensor-sedghi2016')
                Solution	= Tensor_Sedghi2016(v,U,Option);                % Tensor factorization is way more complex than what can be captured by EvalCount and RegrCount, omit
            end
            
            err     = v - U*Solution.beta;
            res     = abs(err);
            beta0   = zeros(1,K);
            for k = 1 : K
                drp = zeros(N,1);
                e   = err(:,k);
                r   = res(:,k);
                s   = sum(r);
                parfor index = 1 : N
                    news = sum(min([abs(e-e(index)),r],[],2));
                    drp(index) = s-news;
                end
                [~,pt] = max(drp);
                beta0(k) = e(pt) + ymu - Xmu*Solution.beta(:,k);
            end
            Solution.beta = [beta0;Solution.beta];
            
            % initialize sample weight w
            w       = ones(N,1);
            [~,err]	= EvalSolution(Solution,X1,y);                          EvalCount = EvalCount + 1;
            errsq   = err.^2;
            sigmasq	= sum(errsq,1)./sum(w,1);
            phi     = exp(-(err.^2)./(2*sigmasq));
            w       = phi./sum(phi,2);
            pr      = mean(w,1);
            pr      = pr./sum(pr);
            beta    = Solution.beta;

            Solution.w      = w;
            Solution.pr     = pr;
            Solution.r2     = 0;
            Solution.Score  = 0;            
           
        otherwise       % Random initialization
            
            w       = rand(N,K);
            w       = w./sum(w,2);
            pr      = mean(w,1);
            pr      = pr./sum(pr);
            beta    = [zeros(1,K);randn(p,K)];
            
            Solution.pr     = pr;
            Solution.beta   = beta;
            Solution.w      = w;
            Solution.r2     = 0;
            Solution.Score  = 0;
            
    end
else
    w       = Solution.w;
    beta    = Solution.beta;
    pr      = Solution.pr;
    Solution.r2     = 0;
    Solution.Score  = 0;
end
[~,err]     = EvalSolution(Solution,X1,y);                                  EvalCount = EvalCount + 1;
errsq       = err.^2;
sigmasq     = sum(w.*errsq,1)./sum(w,1);

Solution.sigmasq = sigmasq;

%% Main EM loop

% Record of good solutions
for index = 1 : EliteLen
    BestSolution(index) = Solution;
end
ScoreBoard      = zeros(1,EliteLen);
ptScoreHist     = 0;
ScoreHistory    = zeros(1,MaxIter);


bConverge	= false;
bContinue	= true;
c = 0;
while bContinue
    
    c = c + 1;
    
    % \\------------------ Evaluate model ------------------// 
    
    [yhat,err]      = EvalSolution(Solution,X1,y);                          EvalCount = EvalCount + 1;
    errsq           = err.^2;
    r2              = 1 - sum((y-sum(w.*yhat,2)).^2)/SSy;
    Score           = r2;
    
    Solution.r2     = r2;
    Solution.Score  = Score;    
    
    % Record and report progress
    ptScoreHist     = ptScoreHist+1;
    ScoreHistory(:,ptScoreHist) = Score;
    % Record if found a good solution
    ptScore = find(Score > ScoreBoard, 1, 'first');
    if ~isempty(ptScore)
        BestSolution(ptScore+1:end) = BestSolution(ptScore:end-1);
        ScoreBoard(ptScore+1:end)	= ScoreBoard(ptScore:end-1);
        BestSolution(ptScore)       = Solution;
        ScoreBoard(ptScore)         = Score;
    end
    if Option.bReportStats
        ReportStatistics(Solution, ScoreHistory(1:ptScoreHist), X1, y, Option);
    end
    
    % Check termination criteria    
    Last5Score  = ScoreHistory(max([1,ptScoreHist-10]):ptScoreHist); %Change Minus For Deeper Tests
    if ptScoreHist > 10
        if min(Last5Score(1,:)) > 0
            bCondition1 = mean(abs(diff(Last5Score(1,:))))/max(Last5Score(1,:)) < prctScoreTh;
            bConverge   = bCondition1;
        end
    end
    bContinue = c < MaxIter & ~bConverge;
    
    
    if bContinue % ------------ EM core ------------
        
    % Scale estimation
    z2	= errsq./(sigmasq+eps);
    switch ScaleEstWeight.Name
        case 'tukey'
            v2      = z2/ScaleEstWeight.csq;
            u       = 3-3*v2+v2.^2;
            bOut    = z2 > 1;
            u(bOut) = 1./v2(bOut);
            u       = u/(ScaleEstWeight.csq*ScaleEstWeight.d);
        case 'huber'
            v2      = z2/ScaleEstWeight.csq;
            u       = ones(size(v2));
            bOut    = v2 > 1;
            u(bOut) = (2*sqrt(v2(bOut))-1)./v2(bOut);
            u       = u/(ScaleEstWeight.csq*ScaleEstWeight.d);
        case 'welsch'
            v2      = z2/ScaleEstWeight.csq;
            u       = zeros(size(v2));
            bZero   = v2 < 0.01;
            t       = v2(bZero);
            u(bZero)= 1 - t/2 + t.^2/6 - t.^3/24 + t.^4/120;        % Taylor series expansion near 0
            t       = v2(~bZero);
            u(~bZero)=(1-exp(-t))./t;                               % Standard expression
            u       = u/(ScaleEstWeight.csq*ScaleEstWeight.d);
        case 'cauchy'
            v2      = z2/ScaleEstWeight.csq;	v2(v2>1e6)=1e6;
            u       = zeros(size(v2));
            bZero   = v2 < 0.01;
            t       = v2(bZero);
            u(bZero)= 1 - t/2 + t.^2/3 - t.^3/4 + t.^4/5;           % Taylor series expansion near 0
            t       = v2(~bZero);
            u(~bZero)=log(1+t)./t;                                  % Standard expression
            u       = u/(ScaleEstWeight.csq*ScaleEstWeight.d);
        case 't-robust'
            u       = (1+nu)./(nu + errsq/sigmasq);
        case ''
            u       = ones(size(w));
    end
    sigmasq = sum(w.*u.*errsq,1)./(sum(w,1)+eps);
        
%     % Trim high-leverage outliers
%     switch TrimMethod.Name
%         case 'mcd'
%             % (to be implemented) % find ptTrim
%             effw(ptTrim,:) = 0;
%     end
        
    % \\----------------------- E-step -----------------------//
    
    % Evaluate likelihood function
    switch DistType.Name
        case 't'        % fit generalized student's t-distribution
            % re-estimate nu by maximizing loglikelihood function
            nuList = DistType.nuList;
            phi_test = -inf(length(nuList),1);
            for index = 1 : length(nuList)
                v	= nuList(index);
                a   = 0.5*(v+1);
                f_test          = (gamma(a)./(sqrt(v*sigmasq)*gamma(v/2)))./(1+errsq./(sigmasq*v)).^a;
                phi_test(index)	= sum(log(f_test*pr'),1);
            end
            [~,ptmax] = max(phi_test);
            nu  = nuList(ptmax);
            a   = 0.5*(nu+1);
            phi = (gamma(a)./sqrt(sigmasq))./(sqrt(nu)*gamma(nu/2)*(1+errsq./(sigmasq*nu)).^a);
            
        case 'normal'   % normal distribution
            phi	= exp(-errsq./(2*sigmasq));
    end
    
    % Update class weights
    postllhood	= pr.*phi;
    w           = postllhood./sum(postllhood,2);
    if bCEM % Force classification
        [~,ptmax] = max(w,[],2);
        pt  = sub2ind([N,K],[1:N]',ptmax);
        w   = zeros(size(w));
        w(pt) = 1;
    end    
    Solution.w	= w;
    
    % \\----------------------- M-Step -----------------------//
    
    % Update prior likelihood
    pr = mean(w,1);
    pr = pr./sum(pr);
    
    % Re-weight samples for regression
    switch RegressWeight.Name
        case 'tukey'
            z2      = errsq./(sigmasq*RegressWeight.csq+eps);
            u       = 1 - z2;
            u(u<0)  = 0;
            u       = u.^2;
        case 'huber'
            z2      = errsq./(sigmasq*RegressWeight.csq+eps);
            u       = ones(size(z2));
            u(z2>1) = 1./sqrt(z2(z2>1));
        case 'cauchy'
            z2      = errsq./(sigmasq*RegressWeight.csq+eps);
            u       = 1./(1+z2);
        case 'welsch'
            z2      = errsq./(sigmasq*RegressWeight.csq+eps);
            u       = exp(-z2);
        case 't-robust'
            u       = (nu+1)./(nu + errsq/(sigmasq+eps));
        case ''
            u       = ones(size(w));
    end
    effw	= w.*u;
    
    % Re-estimate beta using the new weight
    newbeta = zeros(size(beta));
    switch RegressMethod.Name
        case 'pls'  % PLS
            optPLS.XThreshold   = RegressMethod.XThreshold;
            optPLS.YThreshold   = RegressMethod.YThreshold;
            optPLS.NumComp      = RegressMethod.NumComp;
            for k = 1 : K
                if sum(effw(:,k)) > eps
                    ymu     = mean(effw(:,k).*y);
                    Xmu     = mean(effw(:,k).*X,1);
                    newbeta(2:end,k) = SIMPLS(X-Xmu,y-ymu,effw(:,k),optPLS);                RegrCount = RegrCount + 1;
                    newbeta(1,k)     = (ymu - Xmu*newbeta(2:end,k))/mean(effw(:,k));
%                 else    % No data point selected, try to partition from the largest cluster
%                     [newbeta(:,k),w] = CarveFromMaxCluster(X1,y,beta,pr,w,k);               EvalCount = EvalCount + 1;
%                     pr = mean(w,1);
%                     pr = pr./sum(pr);
                end
            end
            
        case 'larsen'
            optLARSEN	= RegressMethod;
            for k = 1 : K
                if sum(effw(:,k)) > eps
                    larsenbeta   = larsen(X,y,effw(:,k),optLARSEN);
                    newbeta(:,k) = larsenbeta(:,end);                                       RegrCount = RegrCount + 1;
%                 else    % No data point selected, try to partition from the largest cluster
%                     [newbeta(:,k),w] = CarveFromMaxCluster(X1,y,beta,pr,w,k);               EvalCount = EvalCount + 1;
%                     pr = mean(w,1);
%                     pr = pr./sum(pr);
                end
            end
            
        case 'ls'   % weighted least square solution
            for k = 1 : K
                ymu     = mean(effw(:,k).*y);
                Xmu     = mean(effw(:,k).*X,1);
                v       = y - ymu;
                U       = X - Xmu;
                UU      = (U'*(effw(:,k).*U));
                if 1/cond(UU) > 1e-10
                    newbeta(2:end,k)= UU \ (U'*(effw(:,k).*v));
                    newbeta(1,k)    = (ymu-Xmu*newbeta(2:end,k))/mean(effw(:,k));           RegrCount = RegrCount + 1;
%                 else    % No data point selected, try to partition from the largest cluster
%                     [newbeta(:,k),w] = CarveFromMaxCluster(X1,y,beta,pr,w,k);               EvalCount = EvalCount + 1;
%                     pr = mean(w,1);
%                     pr = pr./sum(pr);
                end
            end
    end
    % momentum stabilized update
    beta	= (1-Momentum)*newbeta + Momentum*(beta);
    
    % Update Solution
    Solution.beta   = beta;
    Solution.pr     = pr;
    end % ------------------- end of EM core -------------------
    
end
ScoreHistory(:,ptScoreHist+1:end) = [];
Solution = BestSolution(1); % Report the highest-score solution
if Option.bReportStats
    ReportStatistics(Solution, ScoreHistory, X1, y, Option);
end

Solution.iterations = [EvalCount,RegrCount]; 
Solution.bConverge = bConverge;

end

%% \\\ ********* ********* Local function space ********* ********* ///

function ReportStatistics(Solution, ScoreHist, X1, y, Option)
    global rgbtriplets 
    if isfield(Option,'ptplot') && length(ScoreHist)>1
        if ~isempty(Option.ptplot)
            ptplot  = Option.ptplot;
            
            [N,~]	= size(X1);
            
            beta    = Solution.beta;
            [pp1,K]	= size(beta);
            p       = pp1-1;
            
            w       = Solution.w;
            r2      = Solution.r2;
            
            [yhat,err]	= EvalSolution(Solution,X1,y);
            
            hProgress = figure(ptplot);    clf;
            set(hProgress,'Name','Progress Statistics');
            
            NumRow = ceil( (4+floor((K+1)/2))/2 );
            
            % Display Score
            subplot(NumRow,2,1:2)
            plot(ScoreHist(1,:),'Linewidth',2)  % plot only r^2
            xlabel('Iterations')
            ylabel('Score')
            title(strcat(['Score=',num2str(Solution.Score(1)), '    R^{2}=',num2str(r2)]), 'interpreter', 'tex')
            yaxishig = max(ScoreHist(1,:));
            yaxislow = max([0,min(ScoreHist(1,:))]);
            if yaxishig > yaxislow
                axis([1,length(ScoreHist(1,:)),yaxislow,yaxishig])
            end
            
            subplot(NumRow,2,3)
            title('Data weight')
            colormap('turbo')
            image(w','CDataMapping','scaled');
            set(gca,'CLim',[0,1]);
            colorbar
            
            for ptMix = 1 : K
                displace = mod(ptMix-1,2)/2;
                subplot(NumRow,2,3+ceil(ptMix/2))
                hold on
                stem([0:p]-displace,beta(:,ptMix),'Color',rgbtriplets(ptMix,:),'Marker','none')
                title('regression coefficients')
                hold off
            end
            
            subplot(NumRow,2,4+floor((K+1)/2))
            HistEdge= linspace(min(min(err)),max(max(err)),31)';
            HistBin = 0.5*(HistEdge(1:end-1)+HistEdge(2:end));
            Freq	= histcounts(err, HistEdge)';
            hb      = bar(HistBin,Freq,1,'EdgeColor','none');
            hb(1).FaceColor = rgbtriplets(1,:);
            title('Error distribution')
            drawnow
        end
    end
end

function [newbeta,w] = CarveFromMaxCluster(X1,y,beta,pr,w,k)
    
    [~,MaxCluster] = max(pr);
    [~,ptMaxCluster] = max(w,[],2);
    bSampleMaxCluster = (ptMaxCluster == MaxCluster);
    
    X1_     = X1(bSampleMaxCluster,:);
    y_      = y(bSampleMaxCluster);
    
    yhat    = X1_*beta(:,MaxCluster);
    err     = yhat - y_;
    berrpos = err>0;
    
    if mean(err) > 0
        bSelectSide = berrpos;
        aerr = err(bSelectSide);
    else
        bSelectSide = ~berrpos;
        aerr = -err(bSelectSide);
    end
    prc = prctile(aerr,[80,95],'method','approximate');
    bSelect = aerr > prc(1) & aerr <= prc(2);
    
    A = X1_(bSelectSide,:);
    A = A(bSelect,:);
    B = y_(bSelectSide);
    B = B(bSelect);
    newbeta     = A\B;
    
    N = length(y);
    pt = 1:N;
    pt = pt(bSampleMaxCluster);
    pt = pt(bSelectSide);
    pt = pt(bSelect)';
    
    temp = w(pt,k);
    w(pt,k) = w(pt,MaxCluster);
    w(pt,MaxCluster) = temp;
    
end
