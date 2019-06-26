function Clean = CleanDir(Directory,varargin)

Unclean = dir(Directory);

Clean = {Unclean(~ismember({Unclean.name},[{'.','..','.DS_Store','._.DS_Store'} varargin])).name};

Clean = Clean(~cellfun(@(c) contains(c,'._'),Clean));

end