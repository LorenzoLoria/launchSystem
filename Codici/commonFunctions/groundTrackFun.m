function [lon,lat] = groundTrackFun(tSpan,rVec,wP,thetaG,plots,square,titleString,settings)

    % DESCRIPTION ---------------------------------------------------------
    % The function calculates the ground track of a satellite given its position 
    % vector over time, angular velocity of the Earth's rotation, and the initial 
    % and latitude over a specified time span. Optionally, the function can also plot 
    % the ground track on a 2D map of the Earth.
    %
    % PROTOTYPE
    % [lon,lat] = groundTrackFun(tSpan,rVec,wP,thetaG,plots,square,titleString)
    %
    % INPUT:
    %   tSpan         Vector of time values for which the 
    %                    satellite's position is evaluated.             [s]
    %   rVec          Matrix with satellite's position vector in the   
    %   wP            Angular velocity of the Earth's rotation      [rad/s]
    %   thetaG        Initial Greenwich sidereal angle                [rad]
    %   plots         Indicating whether to plot the ground track 
    %                     - 1 plot
    %                     - 0 no plot
    %   square        Indicating whether to add a zoomed-in subplot
    %                     - 1 zoom
    %                     - 0 no zoom
    %
    % CONTRIBUTORS:
    % Lorenzo Loria
    %
    % ---------------------------------------------------------------------

    
dec = asin(rVec(:,3)./vecnorm(rVec')');
rightAsc=zeros(1,length(tSpan));

        for i = 1:size(rVec,1)
            if rVec(i,2)/norm(rVec(i,:)) > 0
                rightAsc(i) = acos(rVec(i,1)/(norm(rVec(i,:))*cos(dec(i))));
            else
                rightAsc(i) = 2*pi-acos(rVec(i,1)/(norm(rVec(i,:))*cos(dec(i))));
            end
        end
    
    
    th = thetaG+ wP*tSpan;   
    
    lon = rad2deg(rightAsc-th);
    lon = mod(lon + 180, 360) - 180;
    lat = rad2deg(dec)';
    for i=1:length(lon)-1
        
        if abs(lon(i)-lon(i+1))>=180
        lon=[lon(1:i),nan,lon(i+1:end)];
        lat =[lat(1:i),nan,lat(i+1:end)];
        else 
        end
    end

if plots == 1
    
    Earth2D
    setPlotSettings(title(titleString))
    hold on
    plot(lon(1:100),lat(1:100),'linewidth',2,'Color',settings.color.terracotta)
    plot(lon(101:200),lat(101:200),'linewidth',2,'Color',settings.color.orange)
    plot(lon(201:300),lat(201:300),'linewidth',2,'Color',settings.color.gray)
    %scatter(lon(1),lat(1),"LineWidth",1.5,Marker="o",MarkerEdgeColor='#D62728')
    %scatter(lon(end),lat(end),"LineWidth",1.5,Marker="diamond",MarkerEdgeColor='#D62728')
    
    xlim([-180,180])
    ylim([-90,90])
    legend(["Ground Track" "Starting Point" "Final Point"],'BackgroundAlpha',0.5)
    if square == 1
    axes('position',[.70 .65 .1 .1])
    box on 
    indexOfInterest = (lon > 0.5*(lon(1)+lon(end))-0.4) & (lon < 0.5*(lon(1)+lon(end))+0.4); % range of t near perturbation
    Earth2DZoom
    ax = gca;
    ax.YColor = [1 1 1];
    ax.XColor = [1 1 1];
    ax.XTick = [];
    ax.YTick = [];

    hold on
    plot(lon,lat,"LineWidth",2,'Color','#D62728') % plot on new axes
    scatter(lon(1),lat(1),"LineWidth",1.5,Marker="o",MarkerEdgeColor='#D62728')
    scatter(lon(end),lat(end),"LineWidth",1.5,Marker="diamond",MarkerEdgeColor='#D62728')
    if lon(1)<lon(end)
        
        lim1 = lon(1);
        lim2=lon(end);
    else
        lim1 =lon(end); 
        lim2= lon(1);
    end
    
    if lat(1)<lat(end)
        lim3 = lat(1);
        lim4 = lat(end);
    else
        lim3 = lat(end);
        lim4 = lat(1);
    end
    
    xlim([lim1-0.05,lim2+0.05])
    ylim([lim3-0.005,lim4+0.005])
    hold off
    else
    end
    hold off

else
end