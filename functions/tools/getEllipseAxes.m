%% getEllipseAxes
% Compute the length of the ellipse axes
%
%% Syntax
% |[Ax , Ay] = getEllipseAxes(apothem , Sx_Sy) |
%
%
%% Description
% |[Ax , Ay] = getEllipseAxes(apothem , Sx_Sy) | Description
%
%
%% Input arguments
% |apothem| -_SCALAR_- Distance (mm) between the centre of the polygon and one edge of the polygon
%
% |Sx_Sy| -_SCALAR_- Ratio aX / aY of the ellipse moin axes
%
%
%% Output arguments
%
% |Ax| - _SCALAR_ - Length (mm) of the main axis along the X axis
%
% |Ay| - _SCALAR_ - Length (mm) of the main axis along the Y axis
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [Ax , Ay] = getEllipseAxes(apothem , Sx_Sy)

    if (Sx_Sy > 1)
      %Sx is larger. It must be the apothem
      Ax = apothem;
      Ay = apothem ./ Sx_Sy;
    else
      %Sy is larger. It must be the apothem
      Ax = apothem .* Sx_Sy;
      Ay = apothem;
    end
end
