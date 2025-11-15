
function EarthPlot(r)
    % EarthPlot - plot of the Earth sphere with specific radius and texture
    % Input:
    %   r - Earth radius [km]

    textureImage = imread('EarthTexture.jpg');

    % Genera la superficie sferica
    [theta, phi] = meshgrid(linspace(0, 2*pi, size(textureImage, 2)), linspace(0, pi, size(textureImage, 1)));
    X = r * sin(phi) .* cos(theta);
    Y = r * sin(phi) .* sin(theta);
    Z = r * cos(phi);
    
    surface(X, Y, Z, 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'CData', textureImage);
    axis equal;
    grid on
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
   
end
