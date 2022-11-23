%% getSpotFluence
% Compute the fluence map of a PBS spot from the MCsquare beam data library
%
%% Syntax
% |[FluenceMap , Fmax , SpotData] = getSpotFluence(E0 , Diso2Meas , BDL)|
%
%
%% Description
% |[FluenceMap , Fmax , SpotData] = getSpotFluence(E0 , Diso2Meas , BDL)| Description
%
%
%% Input arguments
% |PtsG| - _SCALAR MATRIX_ -  |PtsG(i,:)=[x,y]| Coordinates (as defined by |x_axis| and |y_axis|) of the i-th point in the grid
%
% |PosCent| -_SCALAR VECTOR_- |PosCent = [x,y]| Position (mm) of the centre of the PBS spot
%
% |SpotData| - _STRUCTURE_ - The information about the PBS spot. See function |getSpotFromBDL|
%
% |Zg| -_SCALAR_- Z coordinate (mm) (in the IEC gantry CS) of the plane in which the spot is to be drwan
%                 Z is positive between isocentre and nozzle. Z is negative downstream to isocentre
%
%% Output arguments
%
% |FluenceMap| - _SCALAR VECTOR_ - |FluenceMap(i)| is the Fluence at position |PtsG(i,:)|
%
% |Fmax| -_SCALAR_- MAximum of the fluence map
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [FluenceMap , Fmax ] = getSpotFluence(PtsG , PosCent ,  SpotData , Zg)

  sigX1 = CourantSnyder(SpotData.SpotSize1x , SpotData.Divergence1x , SpotData.Correlation1x , SpotData.iso2Nozzle - Zg);
  sigY1 = CourantSnyder(SpotData.SpotSize1y , SpotData.Divergence1y , SpotData.Correlation1y , SpotData.iso2Nozzle - Zg);
  sigX2 = CourantSnyder(SpotData.SpotSize2x , SpotData.Divergence2x , SpotData.Correlation2x , SpotData.iso2Nozzle - Zg);
  sigY2 = CourantSnyder(SpotData.SpotSize2y , SpotData.Divergence2y , SpotData.Correlation2y , SpotData.iso2Nozzle - Zg);

  %Compute spot profile from the sum of the 2 Gaussian defined in BDL
  FluenceMap = biNorm(PtsG , SpotData.Weight1 .* 2 .* pi .* sigX1 .* sigY1,  PosCent , sigX1 , sigY1 , 0 , SpotData.SpotTilt) ...
          + biNorm(PtsG , SpotData.Weight2 .* 2 .* pi .* sigX2 .* sigY2,  PosCent , sigX2 , sigY2 , 0 , SpotData.SpotTilt);
  Fmax  = max(FluenceMap,[],'all'); %Find the highest dose
end
