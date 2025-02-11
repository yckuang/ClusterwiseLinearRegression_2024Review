function PutToExcel(Input, Method, Output, filename, elaptime)
% Output data in the correct form to excel.
% 07/09/2023 Ye Chow

initType        = upper(Method.init.Name); 
regrMethod      = upper(Method.regr.Name); 
scaleEstWeight  = upper(Method.scaleweight.Name); 
regrWeight      = upper(Method.regrweight.Name); 
distrType       = upper(Method.distr.Name); 

p = Input.p;
K = Input.K;
N = Input.N;
nonzerodimratio = Input.nonzerodimratio;
dotprod         = Input.dotprod;
outlier         = Input.outlier;
noise           = Input.noise;
covariance      = Input.covariance;
offset          = Input.offset;
tolerance       = Input.tolerance;


meanIterations  = mean(Output.iterations,1);
stdIterations   = std(Output.iterations,[],1);
meanAccuracy    = mean(Output.accuracy);
stdAccuracy     = std(Output.accuracy);
convRate        = mean(Output.converged);
meanRMSE        = mean(Output.rmse);
kurtRMSE        = kurtosis(Output.rmse);
meanTime        = mean(elaptime);
stdTime         = std(elaptime);

%printData = {'Iterations', 'std', 'Accuracy', 'std', 'ConvRate'; ...
%            meanIterations, stdIterations, meanAccuracy, stdAccuracy, convRate}

alldata = { regrMethod initType distrType scaleEstWeight regrWeight ...
            K p nonzerodimratio N/K dotprod outlier noise covariance offset ...
            meanAccuracy stdAccuracy convRate tolerance meanRMSE kurtRMSE ...
            meanIterations(1) stdIterations(1) meanIterations(2) stdIterations(2) meanTime stdTime};

SafeAppend(alldata,filename,'Sheet',strcat('K=',num2str(K)),'WriteMode','append');

end