function figure_makeLegendNice(leg, obj, ax, settings)
    leg.Box             = 'off';
    leg.FontName        = 'Arial';
    leg.FontWeight      = 'normal';
    leg.FontSize        = settings.AxisFontSize;
    leg.TextColor       = settings.AxisFontColor;
    leg.Location        = 'southeast';
    leg.EdgeColor       = 'none';
    leg.Color           = 'none';
    leg.Units           = 'centimeters';
    
    
    
    %set legend line width
    %lineWidthFrac       = settings.legendLineWidth;

    lineh = findobj(obj,'type','line');
    %lineW = (ax.XLim(2) - ax.XLim(1)) * lineWidthFrac;
    lineW = settings.legendLineWidth;
    
    try
        for ll = 1:length(lineh)/2
            curLine             = (ll-1)*2 + 1;
            curX                = lineh(curLine).XData;
            newX                = [curX(2) - lineW, curX(2)];
            
            

            lineh(curLine).XData     = newX;
            lineh(curLine+1).XData   = curX(2) - lineW/2;
            
            dX = newX(1);
        end
    catch
        for ll = 1:length(lineh)/2
            curLine             = (ll-1)*2 + 2;
            curX                = lineh(curLine).XData;
            newX                = [curX(2) - lineW, curX(2)];
            
            lineh(curLine).XData     = newX;
            lineh(curLine-1).XData   = curX(2) - lineW/2;
            dX = newX(1);
        end
    end
        
%set legend position
    pos                 = settings.legendLocation;
    posAx               = ax.Position;
    posLeg              = leg.Position;
    if strcmp(pos,'lowerright')
        left = posAx(1) + posAx(3) - posLeg(3);
        bot  = posAx(2);
    elseif strcmp(pos,'upperright')
        left = posAx(1) + posAx(3) - posLeg(3);
        bot  = posAx(2) + posAx(4) - posLeg(4);
    elseif strcmp(pos,'upperleft')
        left = posAx(1);
        bot  = posAx(2) + posAx(4) - posLeg(4);
    elseif strcmp(pos,'lowerleft')
        left = posAx(1);
        bot  = posAx(2);
    end
    leg.Position        = [left bot posLeg(3:4)];
