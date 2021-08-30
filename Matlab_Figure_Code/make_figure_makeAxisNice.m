
function figure_makeAxisNice(ax, imPos, settings)

ax.Box                         = 'off';
ax.XMinorTick                  = 'off';
ax.YMinorTick                  = 'off';
ax.FontName                    = 'Arial';
ax.FontWeight                  = 'normal';
ax.FontSize                    = settings.AxisFontSize;

ax.Units                       = 'centimeter';
ax.Position                    = imPos;

ax.XLabel.Units                = 'centimeter';
ax.XLabel.FontName             = 'Arial';
ax.XLabel.FontWeight           = 'normal';
ax.XLabel.FontSize             = settings.LabelFontSize;
ax.YLabel.Units                = 'centimeter';
ax.YLabel.FontName             = 'Arial';
ax.YLabel.FontWeight           = 'normal';
ax.YLabel.FontSize             = settings.LabelFontSize;

end
