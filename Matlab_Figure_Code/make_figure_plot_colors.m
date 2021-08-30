function colors = make_figure_plot_colors()
    colors = struct();
    
    colors.g    = [0 1 0]; %0 .6 .2
    colors.r    = [1 0 0];
    colors.DP   = [.5 .5 .9];
    colors.DT   = [.9 .9 0]*0.8;
    colors.b    = [0 0 1];
    colors.k    = [0 0 0];
    
    colors.medianCol = [255 38 0]/255;


    %Colorblind friendly 3 color sheme
    colors.green    = [102  194 165]/255;
    colors.red      = [252  141 98]/255;
    colors.blue     = [141  160 203]/255;
    colors.grey     = [1 1 1]*0.25;
    
    
    %Colorblind friendly 3 color sheme - Vibrant
    colors.greenV    = [27  158 119]/255;
    colors.redV      = [217  95 2]/255;
    colors.blueV     = [117  112 179]/255;
    
    %sequential blue
    colors.blue1     = [189,201,225]/255;
    colors.blue2     = [116,169,207]/255;
    colors.blue3     = [5,112,176]/255;

    colors.MagentaDark    = [208,28,139]/255;
    colors.MagentaLight   = [233,163,201]/255;
    colors.GreenLight     = [161,215,106]/255;
    colors.GreenDark      = [77,172,38]/255;
    colors.Orange        = [253,174,97]/255;
    colors.Blue           = [44,123,182]/255;
    colors.Grey           = [0.3, 0.3, 0.3];

  
end