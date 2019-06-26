function struct = Exist(struct,field,varargin)

if ~isfield(struct, field)
    binary = 0;
elseif isempty(struct.(field))
    binary = 0;
else
    binary = 1;
end

if binary == 0 && isempty(varargin) 
    error(['Struct Contains no field ' field]);
elseif binary == 0 && ~isempty(varargin) 
    struct.(field) = varargin{1};
elseif binary == 1 && ~isempty(varargin) 
    
end

end