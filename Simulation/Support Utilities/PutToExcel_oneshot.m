function PutToExcel_oneshot(Input, Method, Output, filename, elaptime)
% Output data in the correct form to excel.
% 28/11/2023 Ye Chow

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


Iterations      = Output.iterations;
Accuracy        = Output.accuracy;
convRate        = Output.converged;
RMSE            = Output.rmse;
R_theory        = Output.R_theory;
R_em            = Output.R_em;
Time            = elaptime;

alldata = { regrMethod initType distrType scaleEstWeight regrWeight ...
            K p nonzerodimratio N/K dotprod outlier noise covariance offset ...
            Accuracy R_theory R_em convRate tolerance RMSE ...
            Iterations(1) Iterations(2) Time};

%add another row
bLoop = true;
loopcounter = 0;
while bLoop
    [fileID,errmsg] = fopen(filename);
    if isempty(errmsg)
        bLoop = false;
        writecell(alldata, filename, 'Sheet',strcat('K=',num2str(K)),'WriteMode','append');
        fclose(fileID);
    else
        loopcounter = loopcounter + 1;
        pause(0.5);
        if loopcounter > 100
            bLoop = false;
            fprintf('Too many failed attempts to open file. Abandon write operation\n')
        end
    end
end

end