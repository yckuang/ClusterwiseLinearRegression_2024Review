clc
clear

strcolor = 'brcmgk';
dir = cd;   cd('..');   parentdir = cd; cd(dir)
addpath(fullfile(dir,'Support Utilities'))
filepath = dir;


RegrMethod = {'ls', 'larsen', 'pls'};
InitType = {'random', 'kmean', 'incr01', 'tensor-sedghi2016'};
ScaleEsWeight = {'', 'huber', 'tukey', 'welsch', 'cauchy'};
RegrWeight = {'', 'huber', 'tukey', 'welsch', 'cauchy'};

column_names = {'RegrMethod', 'InitType', 'DistrType', 'ScaleEstWeight', 'RegrWeight', ...
                'K', 'p', 'nonzerodimratio', 'N/K', 'dotprod', 'Outlier', 'Noise', 'Covariance', 'Offset', ...
                'Acc_CLR', 'Acc_MLR', 'R_theory', 'R_em_CLR', 'R_em_MLR', 'RMSE_CLR', 'RMSE_MLR', ...
                'time_CLR', 'time_MLR', 'IterEval_CLR', 'IterRegr_CLR', 'IterEval_MLR', 'IterRegr_MLR'};
filename = join(['CEM_individual','.xlsx']);
filename = fullfile(filepath, filename);
writecell(column_names, filename);
fclose('all');


for loop = 1 : 1e5
    
    K = 1+randi(3);
    p = K+randi(40);
    
    Option.MaxIter          = 1000;
    Option.Momentum         = 0.2;
    Option.prctScoreTh      = 1e-2;
    Option.EliteLen         = round(2*K);
    Option.prElite          = 0.1;
    Option.ptplot           = 314;
    Option.bReportStats     = false;
    Option.SmallGroupTh     = 0.1;
    Option.PerturbCount     = 3;%2*Option.EliteLen;
    Option.MaxEliteCombi    = 100;
    
    % ResolveHyperplane parameters
    Option.ResHyp.DeflateThreshold      = 1e-4;
    Option.ResHyp.EigspcOp2Threshold    = 3;
    Option.ResHyp.ComplClusThreshold    = 4;
    Option.ResHyp.QThreshold            = 0.8;
    
    % Incremental initialization parameters
    Option.Incr.g1      = 0.5;
    Option.Incr.g21     = 0.8;
    Option.Incr.g22     = 10;
    Option.Incr.g3      = 10;
    Option.Incr.epsilon = 1e-3;
    Option.Incr.maxiter	= 100;
    Option.Incr.surprisethreshold = 6;
    
    % % DGM optimization parameters
    % Option.DGM.lambda   = 0.1353353;
    % Option.DGM.delta	= 1e-4;
    % Option.DGM.reltol   = 1e-4;
    % Option.DGM.abstol	= 1e-6;
    % Option.DGM.c        = 0.3;      % range (0,1)
    % Option.DGM.maxiter  = 10000;
    % Option.DGM.DGopt    = [];
    % Option.DGM.beta     = 0.8;      % range (0,1)
    
    tolerance   = 0.1;
    
    Method = [];
    for RegrMethodNum = 1 : 1 % ls larsen pls
        
        ScaleEsWeightNum    = 1     % '', 'huber', 'tukey', 'welsch', 'cauchy'
        RegrWeightNum       = 3     % '', 'huber', 'tukey', 'welsch', 'cauchy';
        
        
        
        % Excution for the combination of parameter
        SamplSz_Grp = (randi(45)+5)*p;
        dotprod     = 0.8*rand;     % dot products of beta
        noise       = 0.4*rand;
        offset      = 0;%8*rand;
        N           = K*SamplSz_Grp;
        
        covariance      = 0;
        nonzerodimratio = 1;
        outlier         = 0;
        offsety         = 0;
        NumCorupt       = round(outlier*N);
        
        Input   = SetInput('p', p, 'nonzerodimratio', nonzerodimratio, 'N', N,  'K', K, 'outlier', outlier, ...
                        'dotprod', dotprod, 'covariance', covariance, 'offset', offset, 'noise', noise, ...
                        'offsety', offsety, 'tolerance', tolerance);
        
        [X, y, beta, R_theory] = GenData(Input);
        
        if R_theory > 0.5 
            w                   = rand(N,K);
            w                   = w./sum(w,2);
            pr                  = mean(w,1);
            pr                  = pr./sum(pr);
            IniSolution.pr      = pr;
            IniSolution.beta    = [zeros(1,K);randn(p,K)];
            IniSolution.w       = w;
            IniSolution.r2      = 0;
            IniSolution.Score   = 0;

            % CLR method ==================================================
            Method = Initialize('RegrMethod', RegrMethod{RegrMethodNum}, ...
                                'InitType', 'random',...
                                'ScaleEstWeight', ScaleEsWeight{ScaleEsWeightNum},...
                                'RegrWeight', RegrWeight{RegrWeightNum}, ...
                                'CEM',true, 'DistrType', 't', 'p', p);
            Option = SetOption(Option, Method, K);
            
            tic %robust method        
            Solution    = EM(X,y,IniSolution,Option);        
            elaptime_CLR = toc;

            Solution    = CorrectPermutation(Solution, beta);
            yhat        = sum(Solution.w.*(Solution.beta(1,:) + X*Solution.beta(2:end,:)),2);
            betaerr     = Solution.beta-beta;

            accuracy_CLR= mean(max([zeros(1,K);1-sqrt(sum(betaerr.^2,1)./sum(beta.^2,1))],[],1));
            rmse_CLR    = sqrt(mean((y(NumCorupt+1:end)-yhat(NumCorupt+1:end)).^2));     % Exclude outliers in the calculations
            iter_CLR    = Solution.iterations;

            % Calculate R from EM Solution
            s           = sqrt(Solution.sigmasq);
            yhat        = [ones(N,1),X]*Solution.beta;
            yhat1       = yhat./(s);
            yhat2       = yhat1./(s);
            a2          = 1/sum(1./(s.^2));
            a           = sqrt(a2);
            b           = a2*sum(yhat2,2);
            b2          = b.^2;
            c2          = a2*sum(yhat1.^2,2);
            G           = sqrt(K)*a/prod(s.^(1/K))*exp( -0.5/a2*(c2-b2) );
            Q           = mean(G);
            R_em_CLR	= 1-Q;

            % Reference leastsquare method ===================================
            Method = Initialize('RegrMethod', RegrMethod{RegrMethodNum}, ...
                                'InitType', 'random',...
                                'ScaleEstWeight', ScaleEsWeight{ScaleEsWeightNum},...
                                'RegrWeight', RegrWeight{RegrWeightNum}, ...
                                'CEM',false, 'DistrType', 't', 'p', p);
            Option = SetOption(Option, Method, K);

            tic %reference method        
            Solution    = EM(X,y,IniSolution,Option);
            elaptime_MLR = toc;

            Solution    = CorrectPermutation(Solution, beta);
            yhat        = sum(Solution.w.*(Solution.beta(1,:) + X*Solution.beta(2:end,:)),2);
            betaerr     = Solution.beta-beta;

            accuracy_MLR= mean(max([zeros(1,K);1-sqrt(sum(betaerr.^2,1)./sum(beta.^2,1))],[],1));
            rmse_MLR    = sqrt(mean((y(NumCorupt+1:end)-yhat(NumCorupt+1:end)).^2));     % Exclude outliers in the calculations
            iter_MLR    = Solution.iterations;
            
            % Calculate R from EM Solution
            s           = sqrt(Solution.sigmasq);
            yhat        = [ones(N,1),X]*Solution.beta;
            yhat1       = yhat./(s);
            yhat2       = yhat1./(s);
            a2          = 1/sum(1./(s.^2));
            a           = sqrt(a2);
            b           = a2*sum(yhat2,2);
            b2          = b.^2;
            c2          = a2*sum(yhat1.^2,2);
            G           = sqrt(K)*a/prod(s.^(1/K))*exp( -0.5/a2*(c2-b2) );
            Q           = mean(G);
            R_em_MLR	= 1-Q;

            % Save results ===================================================
            alldata = { upper(Method.regr.Name) 'random' upper(Method.distr.Name) upper(Method.scaleweight.Name) upper(Method.regrweight.Name) ...
                        K p nonzerodimratio N/K dotprod outlier noise covariance offset ...
                        accuracy_CLR accuracy_MLR R_theory R_em_CLR R_em_MLR rmse_CLR rmse_MLR ...
                        elaptime_CLR elaptime_MLR iter_CLR iter_MLR};
            
            SafeAppend(alldata,filename, 'Sheet', 1,'WriteMode','append');
            
        end
        
    end    
end