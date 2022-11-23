%% getCEFmaxThickness
% Compute the maximum thickness of the CEF from the min and max energy layers
%
%% Syntax
% |CEFmaxThick = getCEFmaxThickness(Layers , Spike , ScannerDirectory)|
%
%
%% Description
% |CEFmaxThick = getCEFmaxThickness(Layers , Spike , ScannerDirectory)| Description
%
%
%% Input arguments
%   * |Layers| -_STRUCTURE_- Definition of the position of the PBS spots;
%      * |Layers(L).Energy| -_SCALAR_- Energy (MeV) of the L-th layer
%
% |Spike| -_STRUCTURE_- Description of the propoerties of the CEF spikes
%   * |Spike.MaterialID| - _STRING_ - Name of the CEF material, as defined in the file "plugins\openMCsquare\lib\Materials\list.dat"
%   * |Spike.MinThickness| -_SCALAR_- Thickness (mm) in of the base on which the spikes are built.
%
% |ScannerDirectory| - _STRING_ - Name of the folder containing the definition of the CT scanner properties in MCsquare in folder "plugins\openMCsquare\lib\Scanners"
%
%% Output arguments
%
% |CEFmaxThick| - _SCALAR_ - MAximum thickness (mm) of the CEF
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function CEFmaxThick = getCEFmaxThickness(Layers , Spike , ScannerDirectory)

  N_layer = length([Layers(:).Energy]);
  T_max = max([Layers(:).Energy]);

  water = materialDescription('water');
  [~, ~ , SPRcef] =  getMaterialSPR(Spike.MaterialID , ScannerDirectory); %Relative stopping power of the CEF material
  R_max = energy2range(T_max, water.alpha,water.p) + Spike.MinThickness .* SPRcef ./10; %Add the minimum thickness to the max range

  for k_layer = 1:N_layer
      R = energy2range(Layers(k_layer).Energy, water.alpha,water.p); %alpha and p for WATER
      t_RF(k_layer) = 10 .* (R_max-R) ./ SPRcef;  %Height (mm) of the hexagon/square for k_layer.
  end

  CEFmaxThick = max(t_RF); %Maximum WET of the rnage compensator

end
