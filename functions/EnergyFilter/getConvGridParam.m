%% getConvGridParam
% Get the dimensions of the convolution filter from the spike information
%
%% Syntax
% |[radius0 , step , nrSides] = getConvGridParam(Plan , b)|
%
%
%% Description
% |[radius0 , step , nrSides] = getConvGridParam(Plan , b)| Description
%
%
%% Input arguments
% |Plan| - _struct_ - MIROpt structure with updated information:
%   * |Plan.Spike.intrpCTpxlSize| -_SCALAR_- Lateral spatial resolution (mm) of the 3D printer
%   * |Plan.Spike.SpikeType| -_STRING_- Type of spike to be designed. The centre of the spike corresponds to the BP with smaller range ('up') or the largest range ('down') or randomise pixel column ('random'), apply Gaussian filter to "smear"('smooth'), draw elliptical spike ('ellipse')
%   * |Plan.Beams(b).GridLayout| -_STRING_- Layout of the PBS spots on the grid. Options: HEXAGONAL (default), SQUARE
%   * |Plan.Beams(b).RidgeFilter(k).apothem| -_SCALAR_- Apothem (mm) of the base of the spike
%
% |b| -_SCALAR_- Index of the beam in |Plan.Beams(b)| for which the fluence is to be computed
%
%% Output arguments
%
% |radius0| - _SCALAR_ - Radius of the circle circumscribing the base of the spike
%
% |step| - SCALAR_ - Distance (mm) between 2 pixels of the fluence map
%
% |nrSides| - _INTEGER_ - Number of side of the polygon defining the base of the spike
%
% |latticePeriod| -_SCALAR_- Lateral spacing between spots in mm
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [radius0 , step , nrSides , GridLayout , latticePeriod] = getConvGridParam(Plan , b)

  apothem_base = max([Plan.Beams(b).RidgeFilter(:).a_max],[],'all');

  if (isfield(Plan.Beams(b),'GridLayout'))
    GridLayout =  Plan.Beams(b).GridLayout;
  else
    GridLayout =  'HEXAGONAL'; %default value for MIROPT
  end

  latticePeriod = Plan.Beams(b).SpotSpacing;

  switch GridLayout
    case 'HEXAGONAL'
      nrSides = 6;
      radius0 = max(2.*apothem_base./sqrt(3)); %max Radius of the spike
    case 'SQUARE'
      nrSides = 4;
      radius0 = max(sqrt(2).*apothem_base); %max Radius of the spike
  end

  switch Plan.Spike.SpikeType
  case {'ellipse','smooth'}
      nrSides = 1;
      radius0 = max(apothem_base);
  end

  step = Plan.Spike.intrpCTpxlSize; %Lateral resolution of the fluence map

end
