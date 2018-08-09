function savelog = ParRefresh(savelog,varargin)

if ~isempty(savelog)
    dbstop in ParRefresh at 6 if ~isempty(savelog.Error);
    
    ImaginaryCheck = 0;
    
    if length(varargin) == 1
        type = varargin{1};
    else
       error('Need a Type') 
    end
    
    if strcmpi(type,'Read')
        finishedjobs = find([savelog.Read] == 1);
        delete(savelog(finishedjobs));
        savelog = savelog(~finishedjobs);
    elseif strcmpi(type,'Finished')
        finishedjobs = find(ismember({savelog.State},'finished'));
        delete(savelog(finishedjobs));
        savelog = savelog(~finishedjobs);
    else
        error('WRONG Type!')
    end
else
    savelog = parfeval(@sum,1,1,1);
    fetchNext(savelog);
end

end