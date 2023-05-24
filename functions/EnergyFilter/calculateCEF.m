%% calculateCEF
% Compute the shape and position of the spikes of the Conformal Energy Filter
% so as to create the same SOBP with Conformal Energy Filter than what is specified in the
% weight of the spots from different energy layers.
% There is one Conformal Energy Filter spike at each (x,y) location of a PBS spot.
% The spike is a stack of hexagons or squares with decreasing apothem.
% The crosssection area of the hexagon/square is proportional to the weight of the spot in the layer L
% The height from the base to the top of the hexagon/square L is proportional to the range of the BP for energy L
%
% The position of the spike is defined for the isocentre to range
%% Syntax
% |[RF , R_max , t_RF , IsocenterToRangeModulatorDistance , mag] = calculateCEF(Beams, apothem, Spike , SpotPositions , GridLayout , BDL , ScannerDirectory , nrSides)|
%
%
%% Description
% |[RF , R_max , t_RF , IsocenterToRangeModulatorDistance , mag] = calculateCEF(Beams, apothem, Spike , SpotPositions , GridLayout , BDL , ScannerDirectory , nrSides)| Description
%
%
%% Input arguments
%
% |Beam| -_STRUCTURE_- Information about the beam
%   * |Beams.sigmaAtCEF| -_STRUCTURE_- Sigma of the lateral Gaussian dose distribution in a PBS spot at the base of the CEF
%       * |Beams.sigmaAtCEF.Sx| -_SCALAR_- Sigma (mm) of the lateral Gaussian dose distribution in a PBS spot at the base of the CEF along the X axis
%       * |Beams.sigmaAtCEF.Sy| -_SCALAR_- Sigma (mm) of the lateral Gaussian dose distribution in a PBS spot at the base of the CEF along the Y axis
%       * |Beams.sigmaAtCEF.r| -_SCALAR_- Correlation between X and Y
%       * |Beams.sigmaAtCEF.sigma| -_SCALAR_- Sigma (mm) of the lateral Gaussian dose distribution in a PBS spot at the base of the CEF (assuming circular spot)
%   * |Beams.Layers| -_STRUCTURE_- Definition of the position of the PBS spots;
%      * |Layers(L).Energy| -_SCALAR_- Energy (MeV) of the L-th layer
%      * |Layers(L).SpotPositions(k)| -_SCALAR VECTOR_- =[x,y] position (mm) of the k-th spot in the layer
%      * |Layers(L).SpotWeights(k)| -_SCALAR VECTOR_-  Weight of the k-th spot in the layer
%
% |spotSigma| -_SCALAR_- Sigma (mm) of the lateral Gaussian dose distribution in a PBS spot at the base of the CEF
%
% |apothem| -_SCALAR_- Radius (mm) of the base of the spikes
%
% |Spike| -_STRUCTURE_- Description of the propoerties of the CEF spikes
%   * |Spike.MaterialID| - _STRING_ - Name of the CEF material, as defined in the file "plugins\openMCsquare\lib\Materials\list.dat"
%   * |Spike.MinThickness| -_SCALAR_- Thickness (mm) in of the base on which the spikes are built.
%   * |Spike.intrpCTpxlSize| -_SCALAR_- Resolution (mm) in the thickness of the deosited layer by the 3D printer
%   * |Spike.PreOptimization| -_STRING_- Algorithm to pre optimise spike shape: uniform , Gaussian
%
% |spotPosition| - SCALAR MATRIX_ - |spotPosition(i,:) = [x,y]| The i-th spot is located at position [x,y] in the BEV coordinates
%
% |GridLayout| -_STRING_- Layout of the PBS spots on the grid. Options: HEXAGONAL (default), SQUARE
%
%  |BDL| -_STRING_- Full path to the Beam data library
%
% |ScannerDirectory| - _STRING_ - Name of the folder containing the definition of the CT scanner properties in MCsquare in folder "plugins\openMCsquare\lib\Scanners"
%
% |nrSides| - _INTEGER_ - Number of side of the polygon defining the base of the spike
%
%
%% Output arguments
%
% |RF| -_VECTOR of STRUCT_- Structure describing the shape of the Conformal Energy Filter filer
% |RF(k).GridLayout| -_STRING_- Lattice layout of the spikes 'HEXAGONAL' or 'SQUARE'
% |RF(k).x_centre| -_SCALAR_- X coordinate (mm in IEC Gantry CS) of central axis of the k-th spike projected in the plane of the range modulator
% |RF(k).y_centre| -_SCALAR_- Y coordinate (mm in IEC Gantry CS) of central axis of the k-th spike projected in the plane of the range modulator
% |RF(k).x_spot| -_SCALAR_- X coordinate (mm in IEC Gantry CS) of PBS spot of the k-th spike in the isocentre plane
% |RF(k).y_spot| -_SCALAR_- Y coordinate (mm in IEC Gantry CS) of PBS spot of the k-th spike in the isocentre plane
% |RF(k).w| -_SCALAR_- Weight of the PBS spot of **maximum** energy to deliver at the base of spike = weight of the SOBP. Proportional to the number of protons of maximum energy to deliver
% |RF(k).h_outside| -_SCALAR_- Height (mm) of the base of the k-th spike
% |RF(k).a_max| -_SCALAR VECTOR_- |RF(k).a_max(L)| Max apothem (mm) of the L-th polygon (apothem of the L-th hexagon/square)
% |RF(k).a_min| -_SCALAR VECTOR_- |RF(k).a_min(L)| Min apothem (mm) of the L-th polygon (apothem of the (L+1)-th hexagon/square)
% |RF(k).h_step| -_SCALAR VECTOR_- |RF(k).h_step(L)| Height (mm) of the L-th polygon from the base
% |RF(k_RF).w_step(L)| -_SCALAR VECTOR_-  Weight of the L-th energy layer for the  at location k, expressed in unit proportional to number of protons at **maximum** energy
% |RF(k_RF).Energy| -_SCALAR VECTOR_- |RF(k).Energy(L)| Energy (MeV) of the protons coming out of the L-th ring from the base
% |RF(k_RF).latticePeriod| -_SCALAR_- Lateral spacing between spots (in mm) in the plane of the CEF. For anamorphic projection, takes the smallest period
%
% |R_max| -_SCALAR_- Range (cm) of the spots to deliver on the base of the ridge filter
%
% |t_RF| -_SCALAR VECTOR_- Thickness (mm) of the different steps of the spikes, corresponding to each energy layer
%
% |IsocenterToRangeModulatorDistance| -_SCALAR_- Distance (mm) from isocentre to the base of the CEF.
%
% |mag| - _SCALAR VECTOR_ - Magnification factor [Mx,My] to project from the isocentre plane to the projection plane
%
% REFERENCES
% [1] PBS 10s Concept Performance Assessment by Jonathan Hubeau MID 29367
%
%% Contributors
% Authors : Jonathan Hubeau (open.reggui@gmail.com); Rudi Labarbe (for wrapper function); Lucian Hotoiu (for geometrical variations)

function [RF , R_max , t_RF , IsocenterToRangeModulatorDistance , mag] = calculateCEF(Beams, apothem, Spike , SpotPositions , GridLayout , BDL , ScannerDirectory , nrSides)

%Find the overall dimension of the spot map, over all layers
Nspots = size(SpotPositions,1); %Number of SOBP
N_layer = length([Beams.Layers(:).Energy]);
T_max = max([Beams.Layers(:).Energy]);

water = materialDescription('water');
[~, ~ , SPRcef] =  getMaterialPropCT(Spike.MaterialID , ScannerDirectory); %Relative stopping power of the CEF material
R_max = energy2range(T_max, water.alpha,water.p); %Range in water of the maximum energy

%Put weight at right location  = > step width for given energy
w_RF = zeros(Nspots,N_layer);

%The spot weights are renormalised at max energy to match the number of protons required at 230MeV at the base of each step of the spike
%MCsquare compute the weight in MU. However, the link between MU and proton charge is a function of the proton energy
%We want to compute the step size based on number of protons
Plan.Beams.Layers = Beams.Layers; %Create a fake plan to call the function
Layers0 = Plan.Beams.Layers; %The original weight
Plan = weight2charge(Plan, BDL , T_max); %convert the weight from the energy of the l-th layer to the **maximum** energy
Layers = Plan.Beams.Layers; %The weight are renormalised at max energy.

for k_layer = 1:N_layer
  %Loop for each layer
    if (size(Layers(k_layer).SpotPositions,2)>1)
      %There are several spots in the layer
      NbSpInLayer = length(Layers(k_layer).SpotPositions(:,1));
    else
      %There is a single spot in the lyer
      NbSpInLayer = 1;
      SpotPositions = SpotPositions';
    end
    for k = 1:NbSpInLayer
        %Loop for each spot in layer
        [~,ix] = min(sum((SpotPositions - repmat(Layers(k_layer).SpotPositions(k,:),Nspots,1)).^2,2)); %Find the SOBP position (ix) that is closest to position of spot k
        w_RF(ix,k_layer) = Layers(k_layer).SpotWeights(k); %w_RF(i,L) weight of the hexagon/square at SpotPositions(i,:) for energy layer L
        w_RF0(ix,k_layer) = Layers0(k_layer).SpotWeights(k);
    end
end

%Calculate total weight for each spike (=SOBP)
tw_RF = sum(w_RF,2); %Sum over the layers

%Calculate steps width assuming same base for all spikes
switch nrSides
case 4 %'SQUARE'
    BaseArea = (2.*apothem).^2; % area of base for a square grid
  case 6 %'HEXAGONAL'
    %This is an hexagonal grid
    BaseArea = 2*sqrt(3) .* apothem.^2; %Area of an hexagon
  case 1 %Ellipse
    [Ax , Ay] = getEllipseAxes(apothem , Beams.sigmaAtCEF.Sx ./ Beams.sigmaAtCEF.Sy);
    BaseArea = pi .* Ax .* Ay; %Area of an ellipse
end %switch

%Renormalise the weights so that the sum of weight is equal to 1
% Indeed the integral of the Gaussian over all space =1
w_RF = w_RF ./ repmat(tw_RF,1, N_layer);
wNaN = find(isnan(w_RF));
w_RF(wNaN) = 0; %The layers with all weights =0

%Calculate required thickness for each energy
t_RF = zeros(N_layer,1);
for k_layer = 1:N_layer
    R = energy2range(Layers(k_layer).Energy, water.alpha,water.p); %alpha and p for WATER
    % (R_max-R) is the difference of energy before the range shifter. To compute the height of the CEM steps, it is this DELTA of energy
    % that matters. We do not need to subtract the range shifter shift from the two terms R_max and R.
    CEMthick = getRSThickness( T_max , Layers(k_layer).Energy , Spike.MaterialID);
    t_RF(k_layer) = rounding(Beams.CEFbaseThickness + CEMthick , Spike.intrpCTpxlSize); %Height (mm) of the hexagon/square for k_layer. Rounded to |Spike.intrpCTpxlSize|

    %t_RF(k_layer) = rounding(Beams.CEFbaseThickness + 10 * (R_max-R) ./ SPRcef , Spike.intrpCTpxlSize); %Height (mm) of the hexagon/square for k_layer. Rounded to |Spike.intrpCTpxlSize|
                %Beams.CEFbaseThickness : The base of the range shifter acts as a residual range shifter. Units are mm;
end

%Compute magnification for projection from isocentre plane to CEF base
IsocenterToRangeModulatorDistance = getIsocenterToRangeModulatorDistance(Beams, BDL);
mag = magnification(IsocenterToRangeModulatorDistance , BDL);
fprintf('Isocenter to Modulator Tray Distance : %f mm \n', IsocenterToRangeModulatorDistance)
fprintf('Magnification factor X : %f \n', mag(1))
fprintf('Magnification factor Y : %f \n', mag(2))

%Calculate geometry
RF = struct;
RF(1).GridLayout = GridLayout;
k_RF = 0;

%loop for every spike along the X and Y axis
for kspt = 1:Nspots

    if tw_RF(kspt) ~= 0
        %There is an SOBP at position (kspt). Add a spike
        k_RF = k_RF+1;

        switch Spike.SpikeType %All spike have the same shape. We cann switch based on type of 1st spike
          case 'fractal'
            RF = CEFfractalSpk(RF , k_RF , kspt , SpotPositions , t_RF , tw_RF , w_RF , w_RF0 , Beams  , Layers , Spike , nrSides , mag , apothem , BaseArea);
          otherwise
            RF = CEFsingleSpkPerSpot(RF , k_RF , kspt , SpotPositions , t_RF , tw_RF , w_RF , w_RF0 , Beams  , Layers , Spike , nrSides , mag , apothem , BaseArea);
        end

    else
        %There is no SOBP at this  position. Do not add a spike here
    end

  end

end
