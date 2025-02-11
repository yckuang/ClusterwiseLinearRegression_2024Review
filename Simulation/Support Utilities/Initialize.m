function Method = Initialize(varargin)
% Initilializes data given initilialisation type, trim type, regression
% method, scale estimate & regression weight, and likelihood distribution
% functions.
% 21/11/2022

% InitType: 'random','kmean','tensor-chaganty2013'
% TrimType: '','MCD'
% RegrMethod: 'ls','larsen','pls'
% ScaleEstWeight: '','huber','tukey','welsch','cauchy','t-robust'
% RegrWeight: '','huber','tukey','welsch','cauchy','t-robust'
% DistrType: 'normal', 't'

%%

Method.init.Name = 'kmean'; 
Method.trim.Name = ''; 
Method.regr.Name = 'ls'; 
Method.scaleweight.Name = ''; 
Method.regrweight.Name = ''; 
Method.distr.Name = 'normal'; 

for k = 1:length(varargin)
    if strcmpi(varargin{k},'InitType')
        InitType = varargin{k+1}; 
        Method.init.Name = InitType; 
    elseif strcmpi(varargin{k},'RegrMethod')
        RegrMethod = varargin{k+1}; 
        Method.regr.Name = RegrMethod; 
    elseif strcmpi(varargin{k},'ScaleEstWeight')
        ScaleEstWeight = varargin{k+1}; 
        Method.scaleweight.Name = ScaleEstWeight; 
    elseif strcmpi(varargin{k},'RegrWeight')
        RegrWeight = varargin{k+1}; 
        Method.regrweight.Name = RegrWeight; 
    elseif strcmpi(varargin{k},'DistrType')
        DistrType = varargin{k+1}; 
        Method.distr.Name = DistrType; 
    elseif strcmpi(varargin{k},'p')
        p = varargin{k+1};
	elseif strcmpi(varargin{k},'cem')
        Method.bCEM = varargin{k+1};
    end        
end


switch InitType
    case 'tensor-sedghi2016'
%         Method.init.FactorizationMethod.Name    = 'rpsf';
        Method.init.FactorizationMethod.Name	= 'tpm';
    case 'tensor-chaganty2013'
%         Method.init.FactorizationMethod.Name    = 'rpsf';
        Method.init.FactorizationMethod.Name	= 'tpm';
        Method.RegressMethod.Name               = 'ls';
    case 'kmean'
  
    case 'random'
    
    case 'incr'

end

switch DistrType
    case 't'
        Method.distr.nuList         = logspace(0,log10(20),30); % parameter of the t-distribution
        Method.distr.nu             = 10;                       % Initial value of nu  
    case 'normal'

end

switch ScaleEstWeight
    case 'huber'
        Method.scaleweight.csq	= 0.7979^2;
        Method.scaleweight.d	= 1;
    case 'tukey'
        Method.scaleweight.csq	= 1.5476^2;
        Method.scaleweight.d	= 0.5;
    case 'welsch'
        Method.scaleweight.csq	= 1.6035^2;
        Method.scaleweight.d	= 0.25;
    case 'cauchy'
        Method.scaleweight.csq	= 2.4027^2;
        Method.scaleweight.d	= 1/7;
end


switch RegrWeight
    case ''

    case 'huber'
        Method.regrweight.csq	= 4.6852^2;
    case 'tukey'
        Method.regrweight.csq	= 1.3449^2;
    case 'welsch'
        Method.regrweight.csq	= 2.9843^2;
    case 'cauchy'
        Method.regrweight.csq	= 2.3848^2;
end
 


switch RegrMethod 
    case 'ls'

    case 'larsen'
        Method.regr.singularthresh = 1e-3;
        Method.regr.stdtype   = 1;
        Method.regr.delta     = 0;  
        Method.regr.stop      = 0;
    case 'pls'
        Method.regr.XThreshold    = 0.99;
        Method.regr.YThreshold    = 0.99;
        Method.regr.NumComp       = p;  
end

end