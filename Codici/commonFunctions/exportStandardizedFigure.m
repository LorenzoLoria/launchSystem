function exportStandardizedFigure(fig, figureName, percTextwidth, lineWidth, varargin)
% exportStandardizedFigure - this function standardizes plot graphics and exports the
%                   figure as a .pdf file
%
% INPUTS:
%            fig - figure variable
%     figureName - figure name, string, (note: do not add '.pdf')
%  percTextwidth - percentage of the linewitdth as indicated in LaTeX,
%                  (70% -> 0.7), double
%      lineWidth - choose line width plot (0 is function decided, any other number is the actual lineWidth)    
%       varargin - optional inputs:
% 		        addMarkers -  logic, if true every lines has a different
%                                    marker, default: true
% 	          changeColors -  logic, if true the colors of the lines will
%                                    be changed, default: true
% 	       changeLineStyle -  logic, if both 'changeLineStyle' and
%                                    'changeColors' are true the lines
%                                    with the same color will have two
%                                    different line styles, default: false
% 			       WHratio - double, width/height ratio; if 'WHratio' is 0
%                                    the ratio will not change,
%                                    default: 0 (current ratio)
%            forcedMarkers - double, number of markers for each line. If it
%                                    is set to 0, all points will have a
%                                    marker; default: 0 (all points)
%                     grid -  logic, true to show grid, false otherwise,
%                                    default: true
%           legendLocation - string, location of the legend; 
%                                    default: 'southoutside'
%        legendOrientation - string, orientation of the legend; 
%                                    default: 'horizontal'
%                exportPDF -  logic, true to export pdf file, default: true
%                exportFIG -  logic, true to export fig file, defalut: false
%       satelliteMapColors -  logic, if true, color palette for satellite
%                                    map is chosen, default: false
%               figurePath - string, path in which the pdf file and/or the
%                                    fig file will be saved, default: ''
%          overwriteFigure -  logic, if true it will overwrite the .pdf
%                                    file and the .fig file with the same 
%                                    name, default: false
%
% --------------------------------EXAMPLE--------------------------------
% >> exportStandardizedFigure(gcf,'nameFig',0.67, 'forcedMarkers', 6, ...
%        'WHratio', 1)
% Figure saved successfully!
%
% text to copy:
% \includegraphics[width=0.67\textwidth]{<add figure path>\nameFig.pdf}
%------------------------------------------------------------------------


%% check input validity
if percTextwidth>1
    error('figure width is larger than the page!')
end

%% Parse input
p = inputParser;

addParameter(p, 'addMarkers', true, @islogical);
addParameter(p, 'forcedMarkers', 0);
addParameter(p, 'changeColors', true, @islogical);
addParameter(p, 'changeLineStyle', false, @islogical);
addParameter(p, 'WHratio', 0);
addParameter(p, 'grid', true, @islogical);
addParameter(p, 'exportPDF', true, @islogical);
addParameter(p, 'exportFIG', false, @islogical);
addParameter(p, 'satelliteMapColors', false, @islogical);
addParameter(p, 'figurePath', '');
addParameter(p, 'overwriteFigure', false, @islogical);
addParameter(p, 'legendLocation', 'southoutside');
addParameter(p, 'legendOrientation', 'horizontal');

parse(p, varargin{:});

%% Recall data
%%% diplay plot config
addMarkers       = p.Results.addMarkers;         % add different markers to lines
changeColors     = p.Results.changeColors;       % change lines colors
changeLineStyle  = p.Results.changeLineStyle;    % to diversify lines of the same color;
                                                 % change colors has to be true
legendLocation   = p.Results.legendLocation;
legendOrientation = p.Results.legendOrientation;

% width\height ratio
if p.Results.WHratio == 0
    changeWHratio = false;  % true = change width/height ratio; false = keep the same ratio
else
    changeWHratio = true;
    WHratio = p.Results.WHratio;
end
% forced markers
if p.Results.forcedMarkers == 0
    forcedMarkers = false;
    nForcedMarkers = nan;
else
    forcedMarkers = true;
    nForcedMarkers = p.Results.forcedMarkers*ones(1,length(fig.Children));
end

if (~addMarkers) && (forcedMarkers)
    warning('addMarker is false and forcedMarkers is true: markers are not displayed')
end

figPath = p.Results.figurePath;
if ~strcmp(figPath,'') && ~(strcmp(figPath(end),'\') || strcmp(figPath(end),'/'))
    figPath = strcat(figPath,'\');
end

%%% dimensions setup
mult = 1.5;             % multiplier for .fig file view
textwidthCm = 16.54;    % \textwidth length in centimeters

fontsize = 9;           % reference font size
fontsizeLegend = fontsize*0.8;

%% colors and options lists
if not(p.Results.satelliteMapColors)
    colors = {
        '#4658A9'               % SkywardBlue
        '#D3212D'               % AmaranthRed
        '#00A86B'               % Jade
        '#FDAF00'               % yellowHoney
        '#E94196'               % Pink
        '#FF5A36'               % PortlandOrange
        '#6CC407'               % AppleGreen
        '#B100FF'               % violet
        '#88ACE0'               % LightCobalt
        };
else
    colors = {
        '#D3212D'               % AmaranthRed
        '#FDAF00'               % yellowHoney
        '#E94196'               % Pink
        '#6CC407'               % AppleGreen
        '#B100FF'               % violet
        '#FF5A36'               % PortlandOrange
        '#4658D9'               % SkywardBlue (more blue)
        '#04B96D'               % Jade (lighter)
        '#80A9DD'               % LightCobalt (darker)
        };
end

markers = {'o'
    '<'
    'square'
    'p'
    '+'
    '>'
    'h'
    'diamond'
    'v'
    '+'
    '*'
    '^'
    'x'};

linestylesList = {'-'
    '--';
    ':';
    '-.'};

nColors = length(colors);
nMarkers = length(markers);
nStyles = length(linestylesList);

%% figure
f = fig;                                             % figure

widthPos = textwidthCm*percTextwidth*mult;
if ~changeWHratio
    WHratio = f.Position(3)/f.Position(4); % retrieve current WHratio
end
heightPos = widthPos/WHratio;

f.Units = "centimeters";
f.Position(3:4) = [widthPos heightPos];

%%
for k = 1:length(f.Children)                         % for each subfigure
    
    if strcmp(f.Children(k).Type, 'uicontextmenu')  
        continue;   % To avoid errors on ContextMenu
    end

    ax = f.Children(k);                              % axes

    % no minor grid
    if isfield(ax,'MinorGridLineStyle')
        ax.MinorGridLineStyle = 'none';
    end

    % interpreter
    listFieldnames = fieldnames(ax);
    indexInterpreter = find(contains(listFieldnames,'Interpreter'));
    
    if ~isempty(indexInterpreter)
        for i = 1:length(indexInterpreter)
            ax.(listFieldnames{indexInterpreter(i)}) = 'latex';
        end
    end

    for i = 1:length(listFieldnames)
        try
            subfieldNames = fieldnames(ax.(listFieldnames{i})); 
            indexSubInterpreter = find(strcmp(subfieldNames, 'Interpreter')); 
            if ~isempty(indexSubInterpreter)
                for j = 1:length(indexSubInterpreter)
                    ax.(listFieldnames{i}).(subfieldNames{indexSubInterpreter(j)}) = 'latex'; 
                end
            end
        catch
        end
    end


    % fontName
    indexFontName = find(contains(listFieldnames,'FontName'));
    if ~isempty(indexFontName)
        for i = 1:length(indexFontName)
            ax.(listFieldnames{indexFontName(i)}) = 'Palatino Linotype';
        end
    end

    % fontSize
    indexFontSize = find(contains(listFieldnames,'FontSize'));
    if ~isempty(indexFontSize)
        removeInd = [find(contains(listFieldnames,'FontSizeMode')); ...
            find(contains(listFieldnames, 'FontSizeMultiplier'))];
        for i = 1:length(indexFontSize)
            if sum(removeInd == indexFontSize(i))==0
                ax.(listFieldnames{indexFontSize(i)}) = fontsize*mult;
            end
        end
    end

    ax.LineWidth = 0.5;
    if isa(ax,'matlab.graphics.illustration.Legend') % check if axes is a legend
        leg = ax;
        leg.Location = legendLocation;
        leg.Orientation = legendOrientation;
        leg.FontSize = fontsizeLegend*mult;
        while leg.Position(3)>0.8 && leg.NumColumns>1
            leg.NumColumns = leg.NumColumns-1;
        end
    elseif isa(ax, 'matlab.graphics.axis.Axes') || isa(ax, 'matlab.graphics.axis.GeographicAxes')
        % grid
        if isa(ax, 'matlab.graphics.axis.Axes') && p.Results.grid 
            ax.XGrid = "on";
            ax.YGrid = "on";
            ax.ZGrid = "on";
        elseif isa(ax, 'matlab.graphics.axis.Axes')
            ax.XGrid = "off";
            ax.YGrid = "off";
            ax.ZGrid = "off";     
        end

        jLines = 0;
        jHistograms = 0;
        for j = 1:length(ax.Children)
            if isa(ax.Children(j),'matlab.graphics.chart.primitive.Line')
                jLines = jLines+1;
                lines(jLines) = ax.Children(j);
            elseif isa(ax.Children(j),'matlab.graphics.chart.decoration.ConstantLine')
                ax.Children(j).LineWidth = 1*mult;
            elseif changeColors && isa(ax.Children(j),'matlab.graphics.chart.primitive.Histogram')
                jHistograms = jHistograms+1;
                histgrams(jHistograms) = ax.Children(j);
            end
        end
        existLegend = isfield(ax,'Legend');
        if existLegend
            leg = ax.Legend;
            leg.Location = legendLocation;
            leg.Orientation = legendOrientation;
            leg.FontSize = fontsizeLegend*mult;
            while leg.Position(3)>0.8
                leg.NumColumns = leg.NumColumns-1;
            end
        end

        %% lines
        if exist('lines','var')
            N = length(lines);
            for i = 1:N
                % changeLineStyle
                if changeLineStyle
                    index = mod(i-1, nStyles) + 1;
                    lines(N+1-i).LineStyle = linestylesList{index};
                end

                % changeColors
                if changeColors
                    index = mod(i-1, nColors) + 1;
                    lines(N+1-i).Color = colors{index};
                end

                % addMarkers
                if addMarkers
                    index = mod(i-1, nMarkers) + 1;
                    lines(N+1-i).Marker = markers{index};
                    lines(N+1-i).MarkerSize = mult*3;
                end

                % forcedMarkers
                nElements = length(lines(N+1-i).MarkerIndices);
                if forcedMarkers && nElements>nForcedMarkers(k)
                    markerIndices = linspace(1,nElements,nForcedMarkers(k));
                    lines(N+1-i).MarkerIndices = round(markerIndices,0);
                end
                if lineWidth == 0
                lines(N+1-i).LineWidth =1.5*mult;              % general settings
                else
                lines(N+1-i).LineWidth = lineWidth;
                end
            end
        end

        if changeColors && exist('histgrams','var')
            N = length(histgrams);
            for i = 1:N
                index = mod(i-1, nColors) + 1;
                histgrams(N+1-i).FaceColor = colors{index};
            end
        end
    end

    clear("lines"); 
end

%% export
if p.Results.exportPDF
    s1 = strcat(figPath, figureName,'.pdf');
    
    if exist(s1,"file")==0 || p.Results.overwriteFigure
        exportgraphics(f,s1,"ContentType","vector")
    else
        time = string(datetime('now','Format','yyyy-MM-dd-HHmmss'));
        s2 = strcat(figPath, figureName, time, '.pdf');
        exportgraphics(f,s2,"ContentType","vector")
        warning('figure named ''%s'' instead of ''%s''',s2,s1)
    end

    fprintf('Figure saved successfully!\n\ntext to copy:\n')
    latexstr = strcat('\includegraphics[width=',num2str(percTextwidth,2),...
        '\textwidth]{<add figure path>\',figureName,'.pdf}');
    disp(latexstr)
end

if p.Results.exportFIG
    s1 = strcat(figPath, figureName,'.jpg');
    if p.Results.overwriteFigure || exist(s1,"file")==0
         exportgraphics(f,s1,"ContentType","auto","Resolution",350)
    else
        time = string(datetime('now','Format','yyyy-MM-dd-HHmmss'));
        s2 = strcat(figPath, figureName, time, '.jpg');
        exportgraphics(f,s2,"ContentType","auto","Resolution",350)
        warning('figure saved as ''%s'' instead of ''%s''',s2,s1)
    end
end
