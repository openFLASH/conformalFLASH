%% makeConvGrid
% Define the X and Y axis coordinates for the cartesian grid used to define the shape of the CEF spike
% The same grid is used for the convolution involgin the propagation of the fluence through the system
%
% |pts_Spot| defines the coordinate of pixels, centered on a PBS spots. The number of pixels defining a spot is |NpxlSpot| and |SpotSize| define the size of the square pixel grid.
% |X_grid| and |Y_grid| define a grid that is large enough to include all the CEF spikes AND a circle with a radius 2*sigma so as to accomodate
% the PBS spot located at the border of the field
%
%% Syntax
% |[X_grid , Y_grid , NpxlSpot , pts_Spot , SpotSize] =  makeConvGrid(spotSigma , step , CEF_x , CEF_y)|
%
%
%% Description
% |[X_grid , Y_grid , NpxlSpot , pts_Spot , SpotSize] =  makeConvGrid(spotSigma , step , CEF_x , CEF_y)| Description
%
%
%% Input arguments
% |spotSigma| -_SCALAR_- Sigma (mm) of the lateral Gaussian dose distribution in a PBS spot (assuming circular spot)
%
% |step| -_SCALAR_- Pixel size (mm)
%
% |CEF_x| -_SCALAR VECTOR_- X coordinate (IEC gantry, mm) of the centre of the spikes
%
% |CEF_y| -_SCALAR VECTOR_- Y coordinate (IEC gantry, mm) of the centre of the spikes
%
%% Output arguments
%
% |X_grid| -_SCALAR VECTOR_- |Xmeas(x)| Cartesian X coordinate (mm) of pixel (x,y) in the grid
%
% |Y_grid| -_SCALAR VECTOR_- |Ymeas(y)| Cartesian Y coordinate (mm) of pixel (x,y)  in the grid
%
% |NpxlSpot| -_SCALAR_- Width (in pixel) of the grid that accomodate the PBS spot
%
% |pts_Spot| -_SCALAR VECTOR_- |pts_Spot(i,:)=[x,y]| Coordinates of the pixels defining one PBS spot
%
% |SpotSize| -_SCALAR MATRIX_- Number of pixels [Nx,Ny] in a spot (i.e. nb of points defined in |pts_Spot|)
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [X_grid , Y_grid , NpxlSpot , pts_Spot , SpotSize] =  makeConvGrid(spotSigma , step , CEF_x , CEF_y)

  radius0 = 2 .* spotSigma; % 95% of fluence

  %make a cartesian grid around a spike
  x_grid = -radius0 + step/2 : step : radius0 + step/2; %Define the coordinate of the pixel at the centre of the pixel -> add step/2
  y_grid = -radius0 + step/2 : step : radius0 + step/2;
  NpxlSpot = numel(x_grid); %Number of pixels along the edge of the cartesian grid
  pts_Spot  = makeListOfCoordinates(x_grid , y_grid);
  SpotSize = [numel(x_grid) , numel(y_grid)]; %Number of pixels (x,y) at each spike

  %coordinate grid where the fluence is to be computed
  % Add a margin
  X_grid = CEF_x(1) + x_grid(1) - step : step : CEF_x(end) + x_grid(end) + step;
  Y_grid = CEF_y(1) + y_grid(1) - step : step : CEF_y(end) + y_grid(end) + step;

end
