%% drawPolygon
% Draw an polygon with a height |h_step| on a square grid of pixels |Zgrid|
% The pixels outside the polygon are left to the initial value defined in |Zgrid|
% See [1] for some basic math information about polygon
%
%% Syntax
% |[Zgrid, GridMask] = drawPolygon(nrSides , pts_Spike , SpikeSize , apothem , h_step , param)|
%
%
%% Description
% |[Zgrid, GridMask] = drawPolygon(nrSides , pts_Spike , SpikeSize , apothem , h_step , param)| Description
%
%
%% Input arguments
% |nrSides| -_INTEGER_- Number of sides to the polygon.
%          * |nrSides = -1| : multiple elliptical spikes per PBS spot
%          * |nrSides = 1| : elliptical spike
%          * |nrSides = 4| : square spike
%          * |nrSides = 6| : hexagonal spike
%
% |pts_Spike| -_SCALAR MATRIX_- |pts_Spike(:,i)= [x,y]| Coordinate (mm) of the i-th pixel in |Zgrid(sub2ind(SpikeSize,i))|
%
% |SpikeSize| -_SCALAR MATRIX_- Number of pixels [Nx,Ny] in |Zgrid|
%
% |apothem| -_SCALAR_- Distance (mm) between the centre of the polygon and one edge of the polygon
%           -_SCALAR VECTOR_- In the case of a rectangle, |apothem(1)| and |apothem(2)| are the half length of each side ofthe rectangle
%
% |h_step| -_SCALAR_- The value of the pixels of |Zgrid| which are inside thehexagon will be replaced by |h_step|
%
% |param| -_STRUCTURE_- [OPTIONAL. Only needed if |nrSides <=1|]
%   * |param.Sx_Sy| -_SCALAR_- Ratio aX / aY of the ellipse moin axes
%   * |param.radius| -_SCALAR_- Radius of the circle in which all the spikes are to be placed
%
%
%% Output arguments
%
% |Zgrid| -_SCALAR MATRIX_- |Zgrid(x,y)| Height (mm) at the pixel (x,y). The matrix has dimension [Nx,Ny]
%
% |GridMask| -_SCALAR MATRIX_- |GridMask(x,y)| = 1 if the pixel at (x,y) is inside the polygon. =0 otherwise. The matrix has dimension [Nx,Ny]
%
%
%% Contributors
% Authors : L. Hotoiu (open.reggui@gmail.com)
%
%
% REFERENCE
% [1] https://www.wikihow.com/Calculate-the-Area-of-a-Hexagon

function [Zgrid, GridMask] = drawPolygon(nrSides , pts_Spike , SpikeSize , apothem , h_step , param)

  switch nrSides

  case -1
        %Draw 4 cylinders
        % GridMask =   GridMask .* 0; %Put all to zero and then draw circles
        Nb = 3;
        dR = 2.* param.radius ./ Nb; %Distance between 2 cylinders
        R = apothem ./ Nb; %Radius of cylinder
        Ori = -floor(Nb./2) .* dR; %Bottom right corner

        for ix = 0:Nb-1
          for iy = 0:Nb-1
            dX = Ori + ix .* dR;
            dY = Ori + iy .* dR;
            GridMask = ( ( (pts_Spike(1,:) - dX ).^2 + (pts_Spike(2,:) - dY ).^2 ) <= R.^2 );
          end
        end

  case 1
        % Draw ellipse
        %Remove all voxels outside of polygon
       [Sx , Sy] = getEllipseAxes(apothem , param.Sx_Sy);
       R = rot(-rad2deg(param.ang) , [0,0,0]); %Rotation matrix around Z
       M = R(1:2,1:2);
       pts_Spike = M * pts_Spike; %Rotate the ellipse
       GridMask = (( (pts_Spike(1,:)./Sx).^2 + (pts_Spike(2,:)./Sy).^2 ) <= 1) ;

  case 4
        % Draw square
        %Set to zero outside the polygon
        GridMask = (pts_Spike(1,:) <= apothem )   .* (pts_Spike(1,:) >= -apothem)    .* (pts_Spike(2,:) >= -apothem)    .* (pts_Spike(2,:) <= apothem );

  case 41
        % Draw rectangle
        %Set to zero outside the polygon
        GridMask = (pts_Spike(1,:) <= apothem(1)) .* (pts_Spike(1,:) >= -apothem(1)) .* (pts_Spike(2,:) >= -apothem(2)) .* (pts_Spike(2,:) <= apothem(2));

  case 6
        % Draw hexagonS
        side = 2.*apothem./sqrt(3); % radius of circle circumscribed to hex = side; needed for diagonal coordinates
        %Remove all voxels outside of polygon
        GridMask = ~( (pts_Spike(1,:) > apothem) + (pts_Spike(1,:) < -apothem) + (pts_Spike(1,:) > (tan(deg2rad(-360/nrSides)) .* (pts_Spike(2,:) - side)) ) + ...
                  (pts_Spike(1,:) < (tan(deg2rad(+360/nrSides)) .* (pts_Spike(2,:) - side)) ) + (pts_Spike(1,:) < (tan(deg2rad(-360/nrSides)) .* (pts_Spike(2,:) + side)) ) + ...
                  (pts_Spike(1,:) > (tan(deg2rad(+360/nrSides)) .* (pts_Spike(2,:) + side)) ) );

  otherwise
        nrSides
        error('Unknown polygon')
  end

    % Set height of grid
    Zgrid = GridMask .* h_step;
end
