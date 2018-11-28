function DatScale(folder,scale)

exp_files = dir([folder filesep '*.vhdr']);
mkdir([folder '_new'])

for k = 1:length(exp_files)

hdr   = ft_read_header([folder filesep exp_files(k).name]);
event = ft_read_event([folder filesep exp_files(k).name]); %Mk<Marker number>=<Type>,<Description>,<Position in data points>,
dat   = ft_read_data([folder filesep exp_files(k).name]);

hdr.orig.MarkerFile = [folder '_new' filesep exp_files(k).name(1:end-5) '.vmrk'];
hdr.orig.DataFile = [folder '_new' filesep exp_files(k).name(1:end-5) '.eeg'];
dat = dat * scale;

ft_write_data([folder '_new' filesep  exp_files(k).name(1:end-5) ],dat,'header',hdr,'dataformat','brainvision_eeg','event',event)
end

end