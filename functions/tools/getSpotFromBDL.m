%% getSpotFromBDL
% Get information about PBS spot properties at a given energy from the Beam Data Library of MCsquare
% The function computes the Courrant Snyder parameters, energy spread and protons per MU at for all
% the specified energies.
% The X and Y axes of the ellipses are rotated by an angle |SpotData.SpotTilt| with respect to the IEC gantry CS
%
%% Syntax
% |SpotData = getSpotFromBDL(BDL_file , E)|
%
%
%% Description
% |SpotData = getSpotFromBDL(BDL_file , E)| Description
%
%
%% Input arguments
% |BDL_file| -_STRING_- path to the BDL file
%
% |En| -_SCALAR VECTOR_- List of energies (MEV) at which the spot properties are requested
%
%
%% Output arguments
%
% |SpotData| - _STRUCTURE_ - The infroamtion about the PBS spot at the specified energy
%   * |SpotData.NominalEnergy| -_SCALAR VECTOR_- the nominal energy (in MeV) of the i-th energy.
%   * |SpotData.iso2Nozzle| -_SCALAR_- Distance (mm) from isocentre to nozzle exit. This is the plane in which the spot has a sigma|SpotData.SpotSizeNx|
%   * |SpotData.EnergySpread(i)| -_SCALAR VECTOR_- the standard deviation of the Gaussian distribution used to model the energy spectrum at the exit of the nozzle (in % of the nominal energy).
%   * |SpotData.ProtonsMU(i)| -_SCALAR VECTOR_- the number of protons delivered per MU at this specific nominal energy.
%   * |SpotData.Weight1(i)| -_SCALAR VECTOR_- the weight of the first Gaussian for the double Gaussian optical model. The sum of Weight1 and Weight2 should be 1.0.
%   * |SpotData.SpotSize1x(i)| -_SCALAR VECTOR_- the in air spot size (in mm) along the x direction at the nozzle exit (for the first 2D Gaussian).
%   * |SpotData.Divergence1x(i)| -_SCALAR VECTOR_- the spot divergence (in rad) along the x direction (for the first 2D Gaussian).
%   * |SpotData.Correlation1x(i)| -_SCALAR VECTOR_- the correlation between the spot size and the divergence along the x direction (for the first 2D Gaussian). Its value must be between -0.99 and 0.99.
%   * |SpotData.SpotSize1y(i)| -_SCALAR VECTOR_- the in air spot size (in mm) along the y direction at the nozzle exit (for the first 2D Gaussian).
%   * |SpotData.Divergence1y(i)| -_SCALAR VECTOR_- the spot divergence (in rad) along the y direction (for the first 2D Gaussian).
%   * |SpotData.Correlation1y(i)| -_SCALAR VECTOR_- the correlation between the spot size and the divergence along the y direction (for the first 2D Gaussian). Its value must be between -0.99 and 0.99.
%   * |SpotData.Weight2(i)| -_SCALAR VECTOR_- the weight of the second Gaussian for the double Gaussian optical model. The sum of Weight1 and Weight2 should be 1.0.
%   * |SpotData.SpotSize2x(i)| -_SCALAR VECTOR_- the in air spot size (in mm) along the x direction at the nozzle exit (for the second 2D Gaussian).
%   * |SpotData.Divergence2x(i)| -_SCALAR VECTOR_- the spot divergence (in rad) along the x direction (for the second 2D Gaussian).
%   * |SpotData.Correlation2x(i)| -_SCALAR VECTOR_- the correlation between the spot size and the divergence along the x direction (for the second 2D Gaussian). Its value must be between -0.99 and 0.99.
%   * |SpotData.SpotSize2y(i)| -_SCALAR VECTOR_- the in air spot size (in mm) along the y direction at the nozzle exit (for the second 2D Gaussian).
%   * |SpotData.Divergence2y(i)| -_SCALAR VECTOR_- the spot divergence (in rad) along the y direction (for the second 2D Gaussian).
%   * |SpotData.Correlation2y(i)| -_SCALAR VECTOR_- the correlation between the spot size and the divergence along the x direction (for the second 2D Gaussian). Its value must be between -0.99 and 0.99.
%   * |SpotData.SpotTilt| -_SCALAR_- Rotation angle (radian) around Ziec axis to go from IEC gantry to the ellipses axes
%
%
%% REFERENCE
% [1] https://gitlab.com/openmcsquare/commissioning/-/blob/master/Tools/load_BDL.m
% [2] http://www.openmcsquare.org/documentation_commissioning.html
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function SpotData = getSpotFromBDL(BDL_file , E)

  BDL = load_BDL(BDL_file);
  SpotData = interpolateAllFields(BDL , E); %Interpolate / extrapolate the value of all BDL parameters at the specified energies
  SpotData = setfield(SpotData,'Energy',E);
  SpotData.iso2Nozzle = getNozzle2IsoDistanceFromBDL(BDL_file);

end

%============================
% Interpolate the value of all the fields present in the BDL structure at the specified energy
%============================
function SpotData = interpolateAllFields(BDL , E)

  SpotData = struct; %The structure with the interpolated values

  NominalEnergy = BDL.NominalEnergy; %Get the nominal energy defined in the BDL
  BDL = rmfield(BDL , 'NominalEnergy'); %Remove nominal intergy so that we do not interpolate on this field

  fieldNames = fieldnames(BDL);
  for f = 1:length(fieldNames)
    %Interpolate the value of each field at the specified energy
    F = getfield(BDL,fieldNames{f});
    Vq = interp1(NominalEnergy,F,E , 'linear','extrap'); %Linear interpolaotion and extrapolation
    SpotData = setfield(SpotData,fieldNames{f},Vq) ;
  end

  if ~isfield(SpotData,'SpotTilt')
    SpotData.SpotTilt = zeros(numel(E),1);
  end


end
