%% CourantSnyder
% Compute the spot sigma of the lateral spread of a PBS spot using the Courant Snyder equation
% The Courant Snyder function is:
% sigma1X^2 =  SpotSize1x^2 -  2 * Correlation1x .*  SpotSize1x .* Divergence1x * Z + Divergence1x^2 * Z^2
%
% The same equation is valid for X and Y for the two bi-normal functions.
%
%
%% Syntax
% |sigX = CourantSnyder(sigX0 , sigTh0 , rho0 , z)|
%
%
%% Description
% |sigX = CourantSnyder(sigX0 , sigTh0 , rho0 , z)| Description
%
%
%% Input arguments
% |sigX0| -_SCALAR VECTOR_- Spot sigma (mm) at the nozzle exit
%
% |sigTh0| -_SCALAR VECTOR_- The beam divergence (radian)
%
% |rho0| -_SCALAR VECTOR_- The correlation coefficient between |sigX0| and |sigTh0|
%
% |z| -_SCALAR VECTOR_- The distance (mm) between the nozzle frame and the plane in which the PBS spot is to be drawn
%
%
%% Output arguments
%
% |sigX| -_SCALAR VECTOR_- |sigX(i)| Spot sigma (mm) in a plane at distance |z| from nozzle frame
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)
%
% REFERENCE
% [1] http://www.openmcsquare.org/documentation_commissioning.html

function sigX = CourantSnyder(sigX0 , sigTh0 , rho0 , z)
  sigX = sqrt(sigX0.^2 - 2 .* rho0 .* sigX0 .* sigTh0 .* z + sigTh0.^2 .* z.^2 );
end
