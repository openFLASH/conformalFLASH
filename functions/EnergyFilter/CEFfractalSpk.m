%% CEFfractalSpk
% Compute the energy filter modulator parameter in the case of a multiple suqre spikes per spot
%
%% Syntax
% |RF = CEFfractalSpk(RF , k_RF , kspt , SpotPositions , t_RF , tw_RF , w_RF , w_RF0 , Beams  , Layers , Spike , nrSides , mag , apothem , BaseArea)|
%
%
%% Description
% |RF = CEFfractalSpk(RF , k_RF , kspt , SpotPositions , t_RF , tw_RF , w_RF , w_RF0 , Beams  , Layers , Spike , nrSides , mag , apothem , BaseArea)| Description
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
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function RF = CEFfractalSpk(RF , k_RF , kspt , SpotPositions , t_RF , tw_RF , w_RF , w_RF0 , Beams  , Layers , Spike , nrSides , mag , apothem , BaseArea)

  N_layer = length([Beams.Layers(:).Energy]);

  RF(k_RF).latticePeriod = Beams.SpotSpacing .* min(mag); %The lattice period projected into the plane of the CEF. Take the smallest magnification
  RF(k_RF).x_centre = SpotPositions(kspt,1) .* mag(1); %Position of the spike in the CEF plane
  RF(k_RF).y_centre = SpotPositions(kspt,2) .* mag(2);
  RF(k_RF).x_spot = SpotPositions(kspt,1); %Position of the PBS spot in the isocentre plane
  RF(k_RF).y_spot = SpotPositions(kspt,2);
  RF(k_RF).w = tw_RF(kspt); %The total weight of the SOBP: sum the charges at position (x,y) of the BP from all layers
  RF(k_RF).h_outside = Spike.MinThickness;
  RF(k_RF).NbColumns = Spike.NbColumns;

  RF(k_RF).h_step = t_RF; %Add the height of the available layers
  RF(k_RF).Energy = [Layers(:).Energy]; %MeV Energy at exit of the step of the CEF

  RF(k_RF).a_max = apothem;
  RF(k_RF).a_min = 0;

  RF(k_RF).t_RF = t_RF; %mm Thickness corrsponding to each energy
  RF(k_RF).w_step = w_RF(kspt,:); %Record the weight of the different energy layers of the spot, expressed in unit proportional to number of protons  at **maximum** energy
  RF(k_RF).w0_step = w_RF0(kspt,:); %Record the weight of the different energy layers of the spot, expressed in unit proportional to number of protons  at **maximum** energy
  RF(k_RF).SpikeType = Spike.SpikeType; %Record the spike type

  RF(k_RF).h_grid = zeros(Spike.NbColumns , Spike.NbColumns); %Prepare an empty grid of spikes

  AngStp = 2.*pi / (4 .* Spike.NbColumns) ; %angular step = 2 pi / perimeter
  angles = 0 : AngStp : 2.*pi; %The sequence of angles

  %At the center of the square, place the tower with the highest weight
  [~ , heavisetStep] = max(w_RF(kspt,:));
  RF(k_RF).h_grid(round(Spike.NbColumns ./ 2),round(Spike.NbColumns ./ 2)) = t_RF(heavisetStep);

  %Around each circle, add the proper number of towers to meet the proportions defined in w_RF(kspt,:)
  % This create a tower distirbution with circular symetry
  % Around each circle, there is the correct proportion of spike height
  for r = 1:round(sqrt(2) .* Spike.NbColumns ./ 2)
    %Loop for every radius from center to edge of cell
    Xs = round (r .* cos(angles) + Spike.NbColumns./2); % Coordinates of all pixels on a circle at radius r. Pixel indices are round number
    Ys = round (r .* sin(angles) + Spike.NbColumns./2);
    outOfBox = logical((Xs > Spike.NbColumns) + (Ys > Spike.NbColumns) + (Xs < 1) + (Ys < 1)); % These points are outside of the square
    Xs(outOfBox) = []; %REmove the points outside of the square
    Ys(outOfBox) = [];

    DistMat = interSpotDistance([Xs',Ys']);
    KeepPts = ones(1,numel(Xs));
    [i,j] = find(DistMat ==0); %find the pair of points at zero distance

    %Make sure that the list contains only unique points. Remove duplicates
    for di = 1:numel(Xs)
      maskDoublon = logical((i == di) .* (j > i)); %These points are duplicates of the i-th point
      KeepPts(j(maskDoublon)) = 0;
    end
    pts = find(KeepPts); %Indices of the pixels where to add a tower

    NbPts = numel(pts); %Nb of points on this ring
    NbSpk_h_step = round(w_RF(kspt,:) .* NbPts); %Number of spikes of the height h_step. Assume UNIFORM fluence throughout the cell
    NbSpk_h_step = assignResid2smallest(NbSpk_h_step , NbPts);  %Adjust the smallest weight to get the correct total number of spikes
    T = [];
    for ti = 1:numel(NbSpk_h_step)
      if(NbSpk_h_step(ti) > 0)
        T = [T , repmat(ti,1,NbSpk_h_step(ti))]; %List of the point indices that are located on the circle of radius r
      end
    end

    idxSpk = randperm(NbPts); %Randomly place the spikes on the circle
    %idxSpk = 1:NbPts; %TODO !!!!!!!!!!!!!!!!! For timing test, make sure to disable the random permutation to reproduce laways the same computation
    idx = sub2ind(size(RF(k_RF).h_grid),  Xs(pts(idxSpk)),Ys(pts(idxSpk))); %Index of the points on the circle
    RF(k_RF).h_grid(idx) = t_RF(T); %Add the towers on the points located on the circle at radius r
  end %for r

[~ , idxStep] = find(w_RF(kspt,:)); %List of the steps with a weight
RF(k_RF).h_grid(RF(k_RF).h_grid==0) = t_RF(min(idxStep)); %If there is one square not assigned to a tower, assign it to the shortest tower


end
