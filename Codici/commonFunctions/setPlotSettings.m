function[] = setPlotSettings(tileTitle)
    %{
        sets the shared settings of all the plots
    %}   
    
    % activate grid and axis
    grid on;
    
    set(gca,'Color','w'); 
    set(gcf,'Color','w');

    % axis colors
    ax = gca;
    % axis color BluePoli
    ax.XColor = '#000000'; ax.YColor = '#000000'; ax.ZColor = '#000000'; 
    
    xaxis= get(gca, 'XAxis');
    xaxis.TickLabelInterpreter = 'latex'; % latex for x-axis
    yaxis= get(gca, 'YAxis');
    
    for i = 1:length(yaxis)
        yaxis(i).TickLabelInterpreter = 'latex';   % tex for y-axis
        yaxis(i).Exponent = 0;
       
        if i ==2
            yaxis(i).Color = '#ca7c6a';
        else
            yaxis(i).Color = '#000000';
        end
    end
    
    zaxis= get(gca, 'ZAxis');
    zaxis.TickLabelInterpreter = 'latex'; % latex for z-axis

    % title of each tile
    tileTitle.Color = '#000000';
    tileTitle.Interpreter = 'latex';
    tileTitle.FontSize = 16;
    tileTitle.FontWeight = "bold";

end
