%% CEFmultiStair
% Compute the energy filter modulator parameter in the case of a multiple staircase per spot
%
%% Syntax
% |RF = CEFmultiStair(RF , k_RF , kspt , SpotPositions , t_RF , tw_RF , w_RF , w_RF0 , Beams  , Layers , Spike , nrSides , mag , apothem , BaseArea)|
%
%
%% Description
% |RF = CEFmultiStair(RF , k_RF , kspt , SpotPositions , t_RF , tw_RF , w_RF , w_RF0 , Beams  , Layers , Spike , nrSides , mag , apothem , BaseArea)| Description
%
%
%% Input arguments
% |RF| -_VECTOR of STRUCT_- Structure describing the shape of the Conformal Energy Filter filer
%
% |k_RF| -_SCALAR_- Index in |RF| to which the spike information is to be added
%
% |kspt| -_SCALAR_- Index of the PBS spot to which the |k_RF| spike corresponds
%
% |spotPosition| - SCALAR MATRIX_ - |spotPosition(i,:) = [x,y]| The i-th spot is located at position [x,y] in the BEV coordinates
%
% |t_RF| -_SCALAR VECTOR_- Thickness (mm) of the different steps of the spikes, corresponding to each energy layer
%
% |tw_RF| -_SCALAR VECTOR_- |tw_RF(i)| Weight of the SOBP allocated to the i-th spike
%
% |w_RF| -_SCALAR MATRIX- w_RF(i,L) weight of the hexagon/square at SpotPositions(i,:) for energy layer L. The weight is rescaled for proton of the maximum evergy
%
% |w_RF0| -_SCALAR MATRIX- w_RF(i,L) weight of the hexagon/square at SpotPositions(i,:) for energy layer L. The weight is assuming the ernergy of layer l
%
% |Beam| -_STRUCTURE_- Information about the beam
%
% |Layers| -_STRUCTURE_- Information from  |Beam.Layers| with the weight renormalised at max energy.
%
% |Spike| -_STRUCTURE_- Description of the propoerties of the CEF spikes
%
% |nrSides| - _INTEGER_ - Number of side of the polygon defining the base of the spike
%
% |mag| - _SCALAR VECTOR_ - Magnification factor [Mx,My] to project from the isocentre plane to the projection plane
%
% |apothem| -_SCALAR_- Radius (mm) of the base of the spikes
%
% |BaseArea| -_SCALAR_- Area (mm2) of the base of the spike
%
%% Output arguments
%
% |RF| -_VECTOR of STRUCT_- Structure describing the shape of the Conformal Energy Filter filer with the |k_RF| updated
%     * |RF.h_grid| -_SCALAR VECTOR_- fraction of area allocated to that step h
%     * |RF.h_step| -_SCALAR VECTOR_- Height of the step h in the staircase
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function RF = CEFmultiStair(RF , k_RF , kspt , SpotPositions , t_RF , tw_RF , w_RF , w_RF0 , Beams  , Layers , Spike , nrSides , mag , apothem , BaseArea)

  N_layer = length([Beams.Layers(:).Energy]);

  RF(k_RF).latticePeriod = Beams.SpotSpacing .* min(mag); %The lattice period projected into the plane of the CEF. Take the smallest magnification
  RF(k_RF).x_centre = SpotPositions(kspt,1) .* mag(1); %Position of the spike in the CEF plane
  RF(k_RF).y_centre = SpotPositions(kspt,2) .* mag(2);
  RF(k_RF).x_spot = SpotPositions(kspt,1); %Position of the PBS spot in the isocentre plane
  RF(k_RF).y_spot = SpotPositions(kspt,2);
  RF(k_RF).apothem = apothem ; %Apothem (mm) of the base of the spikes
  RF(k_RF).w = tw_RF(kspt); %The total weight of the SOBP: sum the charges at position (x,y) of the BP from all layers
  RF(k_RF).h_outside = Spike.MinThickness;
  RF(k_RF).NbColumns = Spike.NbColumns;


  RF(k_RF).Energy = [Layers(:).Energy]; %MeV Energy at exit of the step of the CEF

  RF(k_RF).a_max = apothem;
  RF(k_RF).a_min = 0;

  RF(k_RF).t_RF = t_RF; %mm Thickness corrsponding to each energy
  RF(k_RF).w_step = w_RF(kspt,:); %Record the weight of the spot, expressed in unit proportional to number of protons  at **maximum** energy
  RF(k_RF).w0_step = w_RF0(kspt,:); %Record the weight of the spot, expressed in unit proportional to number of protons  at **maximum** energy
  RF(k_RF).SpikeType = Spike.SpikeType; %Record the spike type


  RF(k_RF).h_grid = w_RF(kspt,:);%fraction of area allocated to that step h
  RF(k_RF).h_grid = RF(k_RF).h_grid ./ sum(RF(k_RF).h_grid); %Normalise to 1
  RF(k_RF).h_step = t_RF; %The height of the corrsponding step


end
