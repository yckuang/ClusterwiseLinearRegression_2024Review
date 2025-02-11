clc
clear

strcolor = 'brcmgk';
dir = cd;   cd('..');   parentdir = cd; cd(dir)
addpath(fullfile(dir,'Support Utilities'))
filepath = dir;


% Set the number of times to execute the file
numExecutions = 500

RegrMethod = {'ls', 'larsen', 'pls'};
InitType = {'random', 'kmean', 'kmean1', 'incr01', 'tensor-sedghi2016'};
ScaleEsWeight = {'', 'huber', 'tukey', 'welsch', 'cauchy', 't-robust'};
RegrWeight = {'', 'huber', 'tukey', 'welsch', 'cauchy', 't-robust'};
DistrType = {'normal', 't'};

%N range
SamplSz_Grp = 500%[100, 500, 1000, 2000];    % sample size per 
NZDR        = [1.0:-0.2:0.4];
DOTPROD     = [0]%[0:0.2:0.8];      % dot products of beta
OUTLIER     = 0%[0.0:0.1:0.2];
NOISE       = [0.1:0.2:0.5];
COVARIANCE	= 0;                % Covariance of X
OFFSET      = [0,2,4]*2;
offsety     = 0;

%K = 2;
% p = 30
for K = 2 : 2

Option.MaxIter          = 1000;
Option.Momentum         = 0.2;
Option.prctScoreTh      = 1e-2;
Option.EliteLen         = round(2*K);
Option.prElite          = 0.1;
Option.ptplot           = 314;
Option.bReportStats     = false;
Option.SmallGroupTh     = 0.1;
Option.PerturbCount     = 3%2*Option.EliteLen;
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

% DGM optimization parameters
Option.DGM.lambda   = 0.1353353;
Option.DGM.delta	= 1e-4;
Option.DGM.reltol   = 1e-4;
Option.DGM.abstol	= 1e-6;
Option.DGM.c        = 0.3;      % range (0,1)
Option.DGM.maxiter  = 10000;
Option.DGM.DGopt    = [];
Option.DGM.beta     = 0.8;      % range (0,1)

bCEM    = true;

SampleSize  = K*SamplSz_Grp;
tolerance   = 0.1;

column_names = {'RegrMethod', 'InitType', 'DistrType', 'ScaleEstWeight', 'RegrWeight', ...
                'K', 'p', 'nonzerodimratio', 'N/K', 'dotprod', 'Outlier', 'Noise', 'Covariance', 'Offset', ...
                'Acc', 'stdAcc', 'convRate', 'tolerance', 'RMSE', 'kurtRMSE', ...
                'IterEval', 'stdIterEval', 'IterRegr', 'stdIterRegr', 'time', 'stdtime'};
            
filename = join(['K',num2str(K),'sparse.xlsx']);
filename = fullfile(filepath, filename);
writecell(column_names, filename, 'Sheet',strcat('K=',num2str(K)));
fclose('all');

for p = 5 : 5 : 80
    
Method = [];
for RegrMethodNum = 1 : 2 % ls larsen pls
    for InitTypeNum = 2 : 2 % 'random', 'kmean', 'kmean1', 'incr01', 'tensor-sedghi2016'
        for ScaleEsWeightNum = 2:2 % '', 'huber', 'tukey', 'welsch', 'cauchy', 't-robust'
            for RegrWeightNum = 2:2     % '', 'huber', 'tukey', 'welsch', 'cauchy', 't-robust';
                for DistrTypeNum = 1:1     % 'normal', 't'
                    
                    Method = Initialize('RegrMethod', RegrMethod{RegrMethodNum}, ...
                                        'InitType', InitType{InitTypeNum},...
                                        'ScaleEstWeight', ScaleEsWeight{ScaleEsWeightNum},...
                                        'RegrWeight', RegrWeight{RegrWeightNum}, ...
                                        'CEM',bCEM, 'DistrType', DistrType{DistrTypeNum}, 'p', p);


                    % Excution for the combination of parameter
                    for ptNZDR = 1: length(NZDR)
                        nonzerodimratio = NZDR(ptNZDR);

                        for ptSize = 1 : length(SampleSize)
                            N = SampleSize(ptSize); 
                            for ptDOTPROD = 1 : length(DOTPROD)
                                dotprod = DOTPROD(ptDOTPROD);

                                for ptOUTLIER = 1: length(OUTLIER)
                                    outlier = OUTLIER(ptOUTLIER);
                                    NumCorupt = round(outlier*N);

                                    for ptNOISE = 1 : length(NOISE)
                                        noise = NOISE(ptNOISE);


                                        for ptCOVARIANCE = 1 : length(COVARIANCE)
                                            covariance = COVARIANCE(ptCOVARIANCE);

                                            for ptOFFSET = 1 : length(OFFSET)
                                                offset = OFFSET(ptOFFSET);                                        

                                                Input = SetInput('p', p, 'nonzerodimratio', nonzerodimratio, 'N', N,  'K', K, 'outlier', outlier, ...
                                                                    'dotprod', dotprod, 'covariance', covariance, 'offset', offset, 'noise', noise, 'offsety', offsety, ...
                                                                    'tolerance', tolerance);

                                                Output = struct('iterations', [], 'accuracy', [], 'converged', []);
                                                Option = SetOption(Option, Method, K);

                                                elaptime    = inf(1,numExecutions);
                                                iterations  = zeros(numExecutions,2);
                                                accuracy	= zeros(numExecutions,1);
                                                converged   = zeros(numExecutions,1);
                                                rmse        = zeros(numExecutions,1);

                                                initmethod = Method.init.Name;

                                                parfor exec = 1: numExecutions
%                                                  for exec = 1: numExecutions
%                                                  disp('parfor is off')

                                                    [X, y, beta]    = GenData(Input);

                                                    tic %start timer
                                                    switch initmethod
                                                        case 'incr'
                                                            Solution    = IncrEM(X,y,[],Option);
                                                        case 'incr01'
                                                            %Solution    = Incr01EM(X,y,[],Option);
                                                            Solution    = Incr01EM01(X,y,[],Option);
                                                        otherwise
                                                            Solution    = EM(X,y,[],Option);
                                                            %Solution    = EM01(X,y,[],Option);
                                                            %Solution    = EM02(X,y,[],Option);
                                                    end
                                                    elaptime(exec) = toc;

                                                    Solution        = CorrectPermutation(Solution, beta);
                                                    yhat            = sum(Solution.w.*(Solution.beta(1,:) + X*Solution.beta(2:end,:)),2);
                                                    betaerr         = Solution.beta-beta;
                                                    acc             = mean(max([zeros(1,K);1-sqrt(sum(betaerr.^2,1)./sum(beta.^2,1))],[],1));
                                                    rmserr          = sqrt(mean((y(NumCorupt+1:end)-yhat(NumCorupt+1:end)).^2));     % Exclude outliers in the calculations

                                                    accuracy(exec)  = acc;
                                                    rmse(exec)      = rmserr;
                                                    iterations(exec,:)= Solution.iterations;
                                                    converged(exec) = 1-acc < Input.tolerance;

                                                end
                                                Output.iterations   = iterations;
                                                Output.accuracy     = accuracy;
                                                Output.converged    = converged;
                                                Output.rmse         = rmse;
                                                % Save results
                                                PutToExcel(Input,Method,Output,filename,elaptime);
                                                
                                            end 
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
            end
        end
    end
end

end
end