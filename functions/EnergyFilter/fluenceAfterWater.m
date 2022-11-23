%% fluenceAfterWater
% Compute the fluence map after a water tank. Take into acount the lateral scattrering of water
%
%% Syntax
% |Fluence = fluenceAfterWater(Fluence , E0 , WaterThickeness , pxlSize)|
%
%
%% Description
% |Fluence = fluenceAfterWater(Fluence , E0 , WaterThickeness , pxlSize)| Description
%
%
%% Input arguments
% |Fluence(x,y,E)| -_SCALAR MATRIX_-  |Fluence(x,y,E)|  fluence at position (x,y) for proton of energy |E|
%
% |E0| -_SCALAR_- Energy (MeV) of the incoming beam
%
% |WaterThickeness| -_SCALAR_- Thickness (mm) of the water tank
%
% |pxlSize| -_SCALAR_- Distance (mm) between points in the fluence map. The fluence map must be equally spaced in x and y
%
% |sigmaWater| -_SCALAR VECTOR_- [OPTIONAL. If absent, the sigma are computd using Moliere angle] |sigmas_far(E)|  Sigma (mm) of the Gaussian kernel for lateral scatter for E-th energy layer
%
%% Output arguments
%
% |Fluence(x,y,E)| -_SCALAR MATRIX_-  |Fluence(x,y,E)|  fluence at position (x,y) for proton of energy |E|
%
% |sigmaWater| -_SCALAR VECTOR_- |sigmas_far(E)|  Sigma (mm) of the Gaussian kernel for lateral scatter for E-th energy layer
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [Fluence , sigmaWater] = fluenceAfterWater(Fluence, E0, WaterThickeness, pxlSize, sigmaWater)

    Wthick = ones(1,size(Fluence,3)) .* WaterThickeness; %Same thickness of water for all energy layer
    distance = 0;
    diverg = 0;
    materiel = 'water';
    if isempty(sigmaWater)
        [~, ~ , thtM] = getSpotSigmaMoliere(E0 , Wthick , distance, materiel); %Get the beam divergence from Moliere scattering
        %Add the variance of each process to obtain the sigma of the Gaussian kernel for the convolution
        % Moliere scattering + intrinsic beam divergence at the base of the CEF + add broadening due to transmission through |distance|
        sigmaWater = sqrt( (thtM .* Wthick).^2 + (thtM .* distance ).^2 + (diverg .* Wthick).^2 + (diverg .* distance ).^2); %sigma in mm
        sigmaWater = round(sigmaWater ./ pxlSize); %sigma in pixels
    end 
    [Fluence , sigmaWater]  = convScatter(Fluence, sigmaWater); %Compute fluence in measurement plan at distance 0mm after water tank

end
