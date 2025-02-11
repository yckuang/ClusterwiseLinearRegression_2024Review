function [errmsg] = SafeAppend(data,filename, varargin)


bLoop = true;
loopcounter = 0;
while bLoop
    [fileID,errmsg] = fopen(filename);
    if isempty(errmsg)
        bLoop = false;
        writecell(data,filename,varargin{:});
    else
        loopcounter = loopcounter + 1;
        pause(0.5);
        if loopcounter > 100
            bLoop = false;
            fprintf('Too many failed attempts to open file. Abandon write operation\n')
        end
    end
    fclose('all');    
end