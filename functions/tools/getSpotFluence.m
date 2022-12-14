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
% |SpotDataNOZL| - _STRUCTURE_ - The information about the PBS spot as saved in the BDL. See function |getSpotFromBDL|. These are the Courant Snyder beam parameter in the plane of the nozzle exit.
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

function [FluenceMap , Fmax ] = getSpotFluence(PtsG , PosCent ,  SpotDataNOZL , Zg)

  %In the BDL, the Courant Snyder parameters are sepcified in the plane of nozzle exit
  % Recompute the parameters in the plane of isocenter
  SpotDataISO = CS_spot_noz2iso(SpotDataNOZL);

  %Use Courant Snyder equation to compute beam parameter in plane Zg. The spot parameter are specified in isocenter plane.
  sigX1 = CourantSnyder(SpotDataISO.SpotSize1x , SpotDataISO.Divergence1x , SpotDataISO.Correlation1x , Zg);
  sigY1 = CourantSnyder(SpotDataISO.SpotSize1y , SpotDataISO.Divergence1y , SpotDataISO.Correlation1y , Zg);
  sigX2 = CourantSnyder(SpotDataISO.SpotSize2x , SpotDataISO.Divergence2x , SpotDataISO.Correlation2x , Zg);
  sigY2 = CourantSnyder(SpotDataISO.SpotSize2y , SpotDataISO.Divergence2y , SpotDataISO.Correlation2y , Zg);

  %Compute spot profile from the sum of the 2 Gaussian defined in BDL
  FluenceMap = biNorm(PtsG , SpotDataISO.Weight1 .* 2 .* pi .* sigX1 .* sigY1,  PosCent , sigX1 , sigY1 , 0 , SpotDataISO.SpotTilt) ...
          + biNorm(PtsG , SpotDataISO.Weight2 .* 2 .* pi .* sigX2 .* sigY2,  PosCent , sigX2 , sigY2 , 0 , SpotDataISO.SpotTilt);
  Fmax  = max(FluenceMap,[],'all'); %Find the highest dose
end
