function dc_select_labels(popup,event,label,varargin)

pos      = get(0,'DefaultFigurePosition');

 
extra = 50;
pos(3:4) = [290+extra+extra 300];
dlg      = dialog('Name', 'Select Labels', 'Position', pos,'WindowStyle','normal');
set(gca, 'Visible', 'off'); % explicitly turn the axis off, as it sometimes appears

userdata.label    = label;

if ~isempty(popup.UserData)
    select            = ismember(label,popup.UserData(:)');     % ensure that it is a row array
else
    select = ~logical(1:length(label));
end

userdata.select   = select;
userdata.unselect = ~select;
set(dlg, 'userdata', userdata); 

if ~isempty(varargin)
    title = varargin{1};
else
    title = {'unselected','selected'};
end


uicontrol(dlg, 'style', 'text',       'position', [ 10 240+20 80+extra+extra  20], 'string', title{1});
uicontrol(dlg, 'style', 'text',       'position', [200+extra 240+20 80+extra  20], 'string', title{2});
uicontrol(dlg, 'style', 'listbox',    'position', [ 10  40+20 80+extra 200], 'min', 0, 'max', 2, 'tag', 'lbunsel') 
uicontrol(dlg, 'style', 'listbox',    'position', [200+extra  40+20 80+extra 200], 'min', 0, 'max', 2, 'tag', 'lbsel') 
uicontrol(dlg, 'style', 'pushbutton', 'position', [105+extra 175+20 80  20], 'string', 'add all >'   , 'callback', @label_addall);
uicontrol(dlg, 'style', 'pushbutton', 'position', [105+extra 145+20 80  20], 'string', 'add >'       , 'callback', @label_add);
uicontrol(dlg, 'style', 'pushbutton', 'position', [105+extra 115+20 80  20], 'string', '< remove'    , 'callback', @label_remove);
uicontrol(dlg, 'style', 'pushbutton', 'position', [105+extra  85+20 80  20], 'string', '< remove all', 'callback', @label_removeall);
uicontrol(dlg, 'style', 'pushbutton', 'position', [ 55+extra  10    80  20], 'string', 'Cancel',       'callback', 'close');
uicontrol(dlg, 'style', 'pushbutton', 'position', [155+extra  10    80  20], 'string', 'OK',           'callback', {@label_submit,dlg,popup});
label_redraw(dlg);

function label_submit(h,events,dlg,pop)
h = get(h, 'parent');
userdata = get(h, 'userdata');
pop.UserData =  userdata.label(userdata.select);
close(dlg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label_redraw(h)
userdata = get(h, 'userdata');
set(findobj(h, 'tag', 'lbsel'  ), 'string', userdata.label(userdata.select));
set(findobj(h, 'tag', 'lbunsel'), 'string', userdata.label(userdata.unselect));
% set the active element in the select listbox, based on the previous active element
tmp = min(get(findobj(h, 'tag', 'lbsel'), 'value'));
tmp = min(tmp, length(get(findobj(h, 'tag', 'lbsel'), 'string')));
if isempty(tmp) | tmp==0
  tmp = 1;
end
set(findobj(h, 'tag', 'lbsel'  ), 'value', tmp);
% set the active element in the unselect listbox, based on the previous active element
tmp = min(get(findobj(h, 'tag', 'lbunsel'), 'value'));
tmp = min(tmp, length(get(findobj(h, 'tag', 'lbunsel'), 'string')));
if isempty(tmp) | tmp==0
  tmp = 1;
end
set(findobj(h, 'tag', 'lbunsel'  ), 'value', tmp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label_addall(h, eventdata, handles, varargin)
h = get(h, 'parent');
userdata = get(h, 'userdata');
userdata.select   = logical(1:length(userdata.label));
userdata.unselect = ~logical(1:length(userdata.label));
set(findobj(h, 'tag', 'lbunsel'  ), 'value', 1);
set(h, 'userdata', userdata);
label_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label_removeall(h, eventdata, handles, varargin)
h = get(h, 'parent');
userdata = get(h, 'userdata');
userdata.unselect = logical(1:length(userdata.label));
userdata.select   = ~logical(1:length(userdata.label));
set(findobj(h, 'tag', 'lbsel'  ), 'value', 1);
set(h, 'userdata', userdata);
label_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label_add(h, eventdata, handles, varargin)
h = get(h, 'parent');
userdata = get(h, 'userdata');
if any(userdata.unselect)
  listbox = get(findobj(h, 'tag', 'lbunsel'  ));
  add = ismember(userdata.label,listbox.String(listbox.Value));
  userdata.select(add)   = 1;
  userdata.unselect(add) = 0;
  set(h, 'userdata', userdata);
  label_redraw(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label_remove(h, eventdata, handles, varargin); 
h = get(h, 'parent');
userdata = get(h, 'userdata');
if any(userdata.select)
  listbox = get(findobj(h, 'tag', 'lbsel'  ));
  remove = ismember(userdata.label,listbox.String(listbox.Value));
  userdata.select(remove)   = 0;
  userdata.unselect(remove) = 1;
  set(h, 'userdata', userdata);
  label_redraw(h);
end

