clc
clear

strcolor = 'brcmgk';
dir = cd;   cd('..');   parentdir = cd; cd(dir)
addpath(fullfile(dir,'Support Utilities'))
filepath = dir;


RegrMethod = {'ls', 'larsen', 'pls'};
InitType = {'random', 'kmean', 'incr01', 'tensor-sedghi2016'};
ScaleEsWeight = {'', 'huber', 'tukey', 'welsch', 'cauchy', 't-robust'};
RegrWeight = {'', 'huber', 'tukey', 'welsch', 'cauchy', 't-robust'};
DistrType = {'normal', 't'};

column_names = {'RegrMethod', 'InitType', 'DistrType', 'ScaleEstWeight', 'RegrWeight', ...
                'K', 'p', 'nonzerodimratio', 'N/K', 'dotprod', 'Outlier', 'Noise', 'Covariance', 'Offset', ...
                'AccGain', 'RMSE Red', 'TimeRatio', 'IterRed_Eval', 'IterRed_Regr', 'R'};
filename = join(['Robust_idv','.xlsx']);
filename = fullfile(filepath, filename);
writecell(column_names, filename);
fclose('all');


for loop = 1 : 1e5
    
    K = 1+randi(3);
    p = max([1+randi(39),K+1]);
    
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
   
    tolerance   = 0.1;
    Method = [];    
    
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
    
    % ---------------- Reference leastsquare method ----------------
    Method = Initialize('RegrMethod', 'ls', 'InitType', 'random', 'CEM',true, 'p', p, ...
                        'ScaleEstWeight', '','RegrWeight', '', 'DistrType', 'normal');
    Option = SetOption(Option, Method, K);
    
    tic %reference method
    Solution    = EM(X,y,IniSolution,Option);
    elaptime_Ref = toc;
    
    Solution    = CorrectPermutation(Solution, beta);
    yhat        = sum(Solution.w.*(Solution.beta(1,:) + X*Solution.beta(2:end,:)),2);
    betaerr     = Solution.beta-beta;
    
    accuracy_Ref= mean(max([zeros(1,K);1-sqrt(sum(betaerr.^2,1)./sum(beta.^2,1))],[],1));
    rmse_Ref    = sqrt(mean((y(NumCorupt+1:end)-yhat(NumCorupt+1:end)).^2));     % Exclude outliers in the calculations
    iter_Ref    = Solution.iterations;

    for ScaleEsWeightNum = 1 : 6     %{'', 'huber', 'tukey', 'welsch', 'cauchy', 't-robust'};
        for RegrWeightNum = 1 : 6     % {'', 'huber', 'tukey', 'welsch', 'cauchy', 't-robust'};
            for DistrTypeNum = 1:2     % 'normal', 't'
                if ~( ((ScaleEsWeightNum==6 || RegrWeightNum==6) && DistrTypeNum==1) || (ScaleEsWeightNum==1 && RegrWeightNum==1 && DistrTypeNum==1) )

                    Method = Initialize('RegrMethod', 'ls', 'InitType', 'random','CEM',true, 'p', p, ...
                                        'ScaleEstWeight', ScaleEsWeight{ScaleEsWeightNum},...
                                        'RegrWeight', RegrWeight{RegrWeightNum}, ...
                                        'DistrType', DistrType{DistrTypeNum});
                    Option = SetOption(Option, Method, K);
                    
                    tic %robust method
                    Solution    = EM(X,y,IniSolution,Option);
                    elaptime = toc;
                    
                    Solution    = CorrectPermutation(Solution, beta);
                    yhat        = sum(Solution.w.*(Solution.beta(1,:) + X*Solution.beta(2:end,:)),2);
                    betaerr     = Solution.beta-beta;
                    
                    accuracy= mean(max([zeros(1,K);1-sqrt(sum(betaerr.^2,1)./sum(beta.^2,1))],[],1));
                    rmse    = sqrt(mean((y(NumCorupt+1:end)-yhat(NumCorupt+1:end)).^2));     % Exclude outliers in the calculations
                    iter    = Solution.iterations;
                    
                    AccGain     = accuracy - accuracy_Ref;
                    rmseRed     = rmse_Ref - rmse;
                    TimeRatio   = elaptime/elaptime_Ref;
                    iterRed     = iter_Ref - iter ;
                    
                    % Save results ===================================================
                    
                    alldata = { 'LS' 'RANDOM' upper(Method.distr.Name) upper(Method.scaleweight.Name) upper(Method.regrweight.Name), ...
                                K p nonzerodimratio N/K dotprod outlier noise covariance offset, ...
                                AccGain, rmseRed TimeRatio iterRed R_theory};
                    
                    SafeAppend(alldata, filename,'WriteMode','append');
                    
                end
            end
        end
    end
    
end