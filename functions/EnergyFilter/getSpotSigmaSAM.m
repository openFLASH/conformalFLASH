%% getSpotSigmaSAM
% Compute the sigma of the lateral spot profile afte going through some thickness of CEF
% The function uses the semi-analytical model
%
%% Syntax
% |sigmaLat = getSpotSigmaSAM(E0 , thickness , distance, matl)|
%
%
%% Description
% |sigmaLat = getSpotSigmaSAM(E0 , thickness , distance, matl)| Description
%
%
%% Input arguments
% |thickness| - _SCALAR VECTOR_ - |thickness(i)| Thickness (mm) of the i-th step of the CEF
%
%
%% Output arguments
%
% |sigmaLat| - _SCALAR VECTOR_ - |sigmaLat(i)| Spot sigma (mm) of the i-th step of the CEF
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function sigmaLat = getSpotSigmaSAM(E0 , thickness )

  pS  = getSAmodelParam();
  sigmaLat = polyval(pS , thickness);

end
