function dc_pmask(mainfig)

ButtonH=uicontrol('Parent', mainfig,'Style','pushbutton','String','p Mask','Position',[10 10 100 25],'Visible','on','Callback',@pmask);
ButtonH.UserData = @ones;

    function pmask(PushButton, EventData)
        
        Axes = PushButton.Parent.Children.findobj('Type','Axes');
        Axes = Axes(arrayfun(@(c) ~isempty(c.String), [Axes.Title]));
        
        for a = 1:length(Axes)
            
            imag = Axes(a).Children.findobj('Type','Image');
            
            if isempty(imag)
                imag =  Axes(a).Children.findobj('Type','Surface');
                
            end
            
            if ~isempty(imag)
                if isempty(imag.UserData)
                    imag.UserData = logical(~isnan(imag.CData));
                end
                
                placeholder = imag.UserData;
                imag.UserData = imag.AlphaData;
                imag.AlphaData = placeholder;
            else
                imag =  Axes(a).Children.findobj('Type','contour');
                placeholder = imag.UserData;
                imag.UserData = imag.ZData;
                imag.ZData = placeholder;
            end
        end
    end
end