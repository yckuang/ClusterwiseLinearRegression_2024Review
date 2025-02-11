function Input   = SetInput(varargin)
% Sets input variables to be used in data generation.
% 25/11/2022
%%
for k = 1:length(varargin)
    if strcmpi(varargin{k},'N')
        Input.N = varargin{k+1}; 
    elseif strcmpi(varargin{k},'p')
        Input.p = varargin{k+1}; 
    elseif strcmpi(varargin{k},'K')
        Input.K = varargin{k+1}; 
	elseif strcmpi(varargin{k},'dotprod')
        Input.dotprod = varargin{k+1};        
    elseif strcmpi(varargin{k},'outlier')
        Input.outlier = varargin{k+1}; 
	elseif strcmpi(varargin{k},'noise')
        Input.noise = varargin{k+1};        
    elseif strcmpi(varargin{k},'covariance')
        Input.covariance = varargin{k+1}; 
    elseif strcmpi(varargin{k},'offset')
        Input.offset = varargin{k+1}; 
    elseif strcmpi(varargin{k},'nonzerodim')
        Input.nonzerodim = varargin{k+1};
    elseif strcmpi(varargin{k},'tolerance')
        Input.tolerance = varargin{k+1};    
    elseif strcmpi(varargin{k},'nonzerodimratio')
        Input.nonzerodimratio = varargin{k+1};
    elseif strcmpi(varargin{k},'offsety')
        Input.offsety = varargin{k+1};
    elseif strcmpi(varargin{k},'repitions')
        Input.repitions = varargin{k+1};
    end
end
end