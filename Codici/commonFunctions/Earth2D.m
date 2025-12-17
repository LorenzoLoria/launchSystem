function Earth2D()
% Function: Earth2D
% Author: Lorenzo Loria
%
% Description:
% This function creates a 2D visualization of the Earth using a texture image 
% of the planet. The image is displayed with longitude and latitude axes, making 
% it suitable as a background for plotting ground tracks or other geographical data.
%
% Inputs:
% - None
%
% Outputs:
% - None
%
% Functionality:
% 1. Loads an image file ('Earth.jpg') representing a 2D map of the Earth.
% 2. Displays the image in a geographic coordinate system, with longitude on the 
%    x-axis and latitude on the y-axis.
% 3. Ensures the image is correctly oriented with latitude increasing from bottom 
%    to top by flipping the image vertically and setting the y-axis direction to normal.
% 4. Adds appropriate axis labels and a title for clarity.
% 5. Configures the plot to maintain equal scaling for longitude and latitude axes.
%
% Usage:
% Earth2D()
%
% Notes:
% - Ensure the file 'Earth.jpg' is available in the current directory or specified path.
% - The function is primarily designed for use as a base map for other plots.
 
    textureImage = imread('Earth_NASA.jpg');
    
   
    xRange = [-180, 180]; 
    yRange = [-90, 90];   

    imagesc(xRange, yRange, flipud(textureImage)); 
    set(gca, 'YDir', 'normal'); 

    xlabel('Longitude [deg]');
    ylabel('Latitude [deg]');
    axis equal
    
    axis on;
    
end
