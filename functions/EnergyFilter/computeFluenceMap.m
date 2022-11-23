%% computeFluenceMap
% Compute nominal fluence map
% Ignore the beam divergence. The spot are placed at the (X,Y) coordinates in the isocentre plane.
% The coordinates of the centre of the spots is the same as the coordinates of the centres of the lattice cell
% Compute lateral scatter through water with WET of the base of CEM + range shifter + patient
%
% The different fluence maps are ordered in different layers with the smae order as |E_flu |
%
%% Syntax
% |Fluence = computeFluenceMap(X_flu , Y_flu , E_flu , Beams , sigmaAtCEF , WaterThickeness , sigmas)|
%
%
%% Description
% |Fluence = computeFluenceMap(X_flu , Y_flu , E_flu , Beams , sigmaAtCEF , WaterThickeness , sigmas)| Description
%
%
%% Input arguments
%
% |X_flu| -_SCALAR VECTOR_- |X_flu(x)| Cartesian X coordinate (mm) of pixel (x,y) in the fluence map
%
% |Y_flu| -_SCALAR VECTOR_- |Y_flu(y)| Cartesian Y coordinate (mm) of pixel (x,y)  in the fluence map
%
% |E_flu| -_SCALAR MATRIX_- |E_flu(E)| is the Energy (MeV) for the fluence map |Fluence(:,:,E)|
%
% |Beams| -_STRUCTURE_-
%   * |Beams.Layers(L).Energy| -_SCALAR_- Energy (MeV) of the L-th layer
%   * |Beams.RidgeFilter| -_VECTOR of STRUCT_- Structure describing the shape of the Conformal Energy Filter filer
%     * |Beams.RidgeFilter(k).x_centre| -_SCALAR_- X coordinate (mm in IEC Gantry CS) of central axis of the k-th spike
%     * |Beams.RidgeFilter(k).y_centre| -_SCALAR_- Y coordinate (mm in IEC Gantry CS) of central axis of the k-th spike
%     * |Beams.RidgeFilter(k).w_step(L)| -_SCALAR VECTOR_-  Weight of the L-th energy layer for the  at location k, expressed in unit proportional to number of protons at **maximum** energy
%     * |Beams.RidgeFilter(k).Energy| -_SCALAR VECTOR_- |RF(k).Energy(L)| Energy (MeV) of the protons coming out of the L-th ring from the base
%
%
% |sigmas| -_STRUCTURE_- [OTPIONAL: if absent, the sigma are computed by the function]. sigma (mm) of the lateral spread of the beam at |distance| of the scatterer
%       * |sigmas.water| -_SCALAR_- Lateral spread introduced by the water
%
%
%% Output arguments
%
% |Fluence(x,y,E)| -_SCALAR MATRIX_-  |Fluence(x,y,E)|  fluence at position (x,y) for proton of energy |E_flu(E)|
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function Fluence = computeFluenceMap(X_flu , Y_flu , E_flu , Beams , sigmas )

  sigmaCEF = Beams.sigmaAtCEF; %mm PBS spot before range shifter and CEF
  CEF_Zg = Beams.RangeModulator.IsocenterToRangeModulatorDistance;

  %Pts = makeListOfCoordinates(X_flu , Y_flu); %List of coordinates wher to compute the fluence
  [pts , Xmap , Ymap] = makeListOfCoordinates(X_flu , Y_flu); %List of coordinates where to compute the fluence
  E0 = max([Beams.Layers(:).Energy],[],'all'); %The maximum energy of incoming protons

  x_spot = [Beams.RidgeFilter(:).x_centre]; %Position of the centre of the spot. The spots centered on the lattice cell, so we can use the cell coordiantes
  y_spot = [Beams.RidgeFilter(:).y_centre];
  Nb_Spots = numel(x_spot); %Number of spikes

  Fluence = zeros(numel(X_flu) , numel(Y_flu) , numel(E_flu));

  %Add the fluence from each spot
  % Return Fluence(x,y,E), the fluence map at each energy layer
  for SpotIdx = 1:Nb_Spots
    w = Beams.RidgeFilter(SpotIdx).w .* Beams.RidgeFilter(SpotIdx).w_step; %Weight of spot wrt to neigbourgh * weight of E layer within the spot
    Fluence = fluenceSngSpot(Fluence , E_flu , w , Beams.RidgeFilter(SpotIdx).Energy , [x_spot(SpotIdx) , y_spot(SpotIdx)] , size(Xmap) , sigmaCEF  , pts , CEF_Zg);
  end

  %Add scatter in water and range shifter to the measurement plan
  pxlSize = min([min(diff(X_flu)),min(diff(Y_flu))]);
  Fluence = fluenceAfterWater(Fluence , E0 , Beams.CEFbaseWET , pxlSize , []); %Scatter in CEM base
  WaterThickeness = Beams.Iso2Skin + Beams.RSinfo.RangeShifterWET;   %WET of the base of CEM + range shifter + patient
  Fluence = fluenceAfterWater(Fluence , E0 , WaterThickeness , pxlSize , []); %Scatter in RS + water

end

%=========================
% Compute the fluence map for all the energy layers of one single spot
%
% INPUT
% |Fluence(x,y,E)| -_SCALAR MATRIX_-  |Fluence(x,y,E)|  fluence at position (x,y) for proton of energy |E_flu(E)|
%
% |E_flu| -_SCALAR MATRIX_- |E_flu(E)| is the Energy (MeV) for the fluence map |Fluence(:,:,E)|
%
% |wSptLayer| -_SCALAR MATRIX- wSptLayer(L) weight of the Bragg peak at position |pos| for energy layer L.
%
% |E_step| -_SCALAR VECTOR_- |E_step(L)| Energy (MeV) of the protons of the L-th energy layer
%
% |pos| -_SCALAR VECTOR_- (x,y) coordinates (mm) of the centre of the spots
%
% |sFmap| -_SCALAR VECTOR_ Size (pixel) of the fluence map in (x,y) dimensions
%
% |sigmaCEF.BDL| -_STRUCTURE_- Spot information from the Beam Data Library. See |getSpotFromBDL| for more information
%
% |PtsG| - _SCALAR MATRIX_ -  |PtsG(i,:)=[x,y]| Coordinates (as defined by |x_axis| and |y_axis|) of the i-th point in the grid
%
% |Zg| -_SCALAR_- Z coordinate (mm) (in the IEC gantry CS) of the plane in which the spot is to be drwan
%=========================
function Fluence = fluenceSngSpot(Fluence, E_flu , w_spot , E_step , pos , sFmap , sigmaCEF   , pts , CEF_Zg)

  for E = 1:numel(E_step)
      G = getSpotFluence(pts , pos , sigmaCEF.BDL , CEF_Zg);
      G = reshape(G,sFmap(1),sFmap(2));
      ind_E = find(E_flu == E_step(E)); %Find the index in the flmuence matrix where to store the fluence map for energy layer E_step(E)
      Fluence(:,:,ind_E) = Fluence(:,:,ind_E) + G .* w_spot(E);
  end
end
