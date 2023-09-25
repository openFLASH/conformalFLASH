%% optimizeTrajectory
% compute the optimum spot trajectory for each beam
% Order the spot in the multi layer PBS plan and then compute the optimum trajectory of the monolayer plan
%
%% Syntax
% |[SpotTrajectoryInfo , Plan] = optimizeTrajectory(Plan  , ROI , plots_BEV)|
%
%
%% Description
% |[SpotTrajectoryInfo , Plan] = optimizeTrajectory(Plan  , ROI , plots_BEV)| Description
%
%
%% Input arguments
% |Plan| - _STRUCTURE_ - The optimised treatment plan
%   * |Plan.Beams(b)| -_STRUCTURE_- Information about the b-th beam
%       * |Plan.Beams(b).spotSigma| -_SCALAR_- Spot lateral sigma (mm) at maximum of deepest Bragg peak along the optical axis
%       * |Plan.Beams(b).NbScarves| -_SCALAR_- Number of scarves to paint on the beam's eye view
%       * |Plan.Beams(b).BDL| -_STRING_- Beam data library. Name of the folder in REGGUI\plugins\openMCsquare\lib\BDL
%   *|Plan.scanAlgoGW| -_STRUCTURE_- [OPTIONAL. Only required to use a scanAlgo Gateway to get the timing] Information about the scanAlgo gateway
%       * |scanAlgoGW.scanalgoGateway_IP| -_STRING_- IP address, including port, to the scanAlkgo gatewat
%       * |scanAlgoGW.room_id| -_STRING_- Room ID as defined  inthe gateway
%       * |scanAlgoGW.snout_id|  -_STRING_- snout ID as defined in the gateway
%       * |scanAlgoGW.spot_id| -_STRING_- Spot tune ID as defined in the gateway
%
% |ROI| - _struct_ - MIROpt structure containing information about all volumes in the RTSTRUCT file. The following data must be present in the structure:
%     * |ROI(i).mask1D| - _array_ - Logical column vector storing a binary mask for ROI i (voxels inside the volume of interest are equal to 1, and those outside are equal to zero).
%
% |plotNb| -_SCALAR VECTOR_- [OTPIONAL] |plotNb(b)| Figure number where to plot the trajectory of beam b
%
%% Output arguments
%
% |SpotTrajectoryInfo| -_STRUCTURE_- Pre-computed spot trajectory. Used to accelerate computations.
%                         |SpotTrajectoryInfo{b}| Parameter for beam |b|
%   * |SpotTrajectoryInfo.beam{b}.sobpSequence| -_SCALAR VECTOR_-  Order of the indices of |spot| to sort the spots. |OrderedSOBP = spot(sobpSequence,:)|
%   * |SpotTrajectoryInfo.beam{b}.Nmaps| -_STRUCTURE_-  Topological information on the initial spot sequence
%       * |Nmaps.NeighbourghMap| -_SCALAR VECTOR_- |NeighbourghMap(d,i)| d=# of delivered spot; i=# of impacted spot; |NeighbourghMap(d,i)|  = 0 if there is an impact
%       * |Nmaps.NeighbourghWeightMap| -_SCALAR VECTOR_- |NeighbourghWeightMap(d,i)| fraction of the dose of spot d (# of delivered spot) that is also delivered at spot i (i=# of impacted spot)
%       * |Nmaps.NeighbourghTimeMap| -_SCALAR MATRIX_- |NeighbourghTimeMap(d,i)| Time (ms) required to go from spot d to spot i
%   * |SpotTrajectoryInfo.sobpPosition{b}| - CELL VECTOR_ - beamletPosition{b}(i,:) = [x,y] The i-th spot of the b-th beam is located at position [x,y] in the BEV coordinates
%   * |SpotTrajectoryInfo.weight2spot| - SCALAR MATRIX_ - |weight2spot(weightIdx,:) = [b,BeamLetNb,l]| The spot at index |weightIdx| is related to the b-th beam and the SOBP # |BeamLetNb| in layer # l
%   * |SpotTrajectoryInfo.scanAlgoGW| -_STRUCTURE_- [OPTIONAL. Present on if |Plan.scanAlgoGW|  is defined ] Information about the scanAlgo gateway
%           * |scanAlgoGW.scanalgoGateway_IP| -_STRING_- IP address, including port, to the scanAlkgo gatewat
%           * |scanAlgoGW.room_id| -_STRING_- Room ID as defined  inthe gateway
%           * |scanAlgoGW.snout_id|  -_STRING_- snout ID as defined in the gateway
%           * |scanAlgoGW.spot_id| -_STRING_- Spot tune ID as defined in the gateway
%
% |Plan| -_STRUCTURE_- Updated data about plan
%     * |Plan.bevPTV| - _CELL VECTOR OF SCALAR MATRIX_ - |Plan.bevPTV{b}(x,y)| Value of the paralell projection of the |target| on a a surface perpendicular to the |beam| axis. The surface is centered on |beam| axis. It is a square with size equal to |bev_size|*|bev_size|
%     * |Plan.bevOAR| - _CELL VECTOR OF SCALAR MATRIX_ - |Plan.bevOAR{b}(x,y)| Value of the paralell projection of the union of all OAR on a a surface perpendicular to the |beam| axis. The surface is centered on |beam| axis. It is a square with size equal to |bev_size|*|bev_size|
%     * |Plan.bev_x| -_CELL VECTOR OF SCLAR VECTOR_- |Plan.bev_x{b}(i)| X coordinate of Plan.bevPTV{b}(i,:) or Y coordinate of Plan.bevPTV{b}(:,i) for b-th beam

%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [SpotTrajectoryInfo , Plan] = optimizeTrajectory(Plan, ROI, plots_BEV, monoEnergyFlag)

  if (nargin < 3)
    for b = 1:size(Plan.Beams,2)
      plots_BEV{b} = [];
    end
  end

  if nargin < 4
      monoEnergyFlag = true;
  end

  if monoEnergyFlag
      %Convert the multi-layer plan into a single layer plan
      [sobpPosition , weight2spot , DRcritical , Plan] =  collectSpotsinBeamlets(Plan , ROI); %Gather the spots into beamlets for form SOBP
    
      SpotTrajectoryInfo.sobpPosition = sobpPosition;
      SpotTrajectoryInfo.weight2spot = weight2spot;
  else
      for b = 1:size(Plan.Beams,2)
          SpotTrajectoryInfo.sobpPosition = Plan.Beams(b).Layers.SpotPositions;
          SpotTrajectoryInfo.weight2spot = Plan.Beams(b).Layers.SpotWeights;
      end
  end


  %Define the way the spot timing is computed
  if isfield(Plan , 'scanAlgoGW')
    SpotTrajectoryInfo.TimingMode = 'scanAlgo';
  else
    SpotTrajectoryInfo.TimingMode = 'Model';
  end

  for b = 1: size(Plan.Beams,2) %Loop for each beam
          sigmaAtIso = max([Plan.Beams(b).spotSigma]); %mm Sigma of the Gaussian PBS spot at isocentre

          if ~isfield(Plan.Beams(b),'NbScarves')
            Plan.Beams(b).NbScarves = []; %Number of scaves is undefined. We will have to compute it
          end

          if min(size(sobpPosition{b}))==1
            %There is a single spot.Trajectory optimisation is simple :-)
            fprintf('One single spot \n')
            SpotTrajectoryInfo.beam{b}.sobpSequence  = 1;
            SpotTrajectoryInfo.beam{b}.Nmaps.NeighbourghMap = 1 ;% -_SCALAR VECTOR_- |NeighbourghMap(d,i)| d=# of delivered spot; i=# of impacted spot; |NeighbourghMap(d,i)|  = 0 if there is an impact
            SpotTrajectoryInfo.beam{b}.Nmaps.NeighbourghWeightMap = 1 ;% -_SCALAR VECTOR_- |NeighbourghWeightMap(d,i)| fraction of the dose of spot d (# of delivered spot) that is also delivered at spot i (i=# of impacted spot)
            SpotTrajectoryInfo.beam{b}.Nmaps.NeighbourghTimeMap = 0;% -_SCALAR MATRIX_- |NeighbourghTimeMap(d,i)| Time (ms) required to go from spot d to spot i

          else
            %There are several spot. Optimize trajectory
                switch SpotTrajectoryInfo.TimingMode
                  case 'scanAlgo'
                    scanAlgoGW = Plan.scanAlgoGW;
                    fprintf('----- Optimising trajectory with scanAlgo gateway: %s \n',Plan.scanAlgoGW.scanalgoGateway_IP)
                    SpotTrajectoryInfo.scanAlgoGW = Plan.scanAlgoGW;

                    if exist('flashAlgo_G4')
                      [sobpSequence , Nmaps] = flashAlgo_G4(sobpPosition{b}, sigmaAtIso , Plan.BDL ,  DRcritical{b} , plots_BEV{b} , Plan.scanAlgoGW); %Find the best scan trajectory for beam b.
                    else
                      [sobpSequence , Nmaps] = flashAlgo_G3(sobpPosition{b}, sigmaAtIso , Plan.BDL ,  DRcritical{b} , Plan.Beams(b).NbScarves , plots_BEV{b} , Plan.scanAlgoGW); %Find the best scan trajectory for beam b.
                    end

                  case 'Model'
                    fprintf('----- Optimising trajectory with simple model \n')
                    if exist('flashAlgo_G4')
                      [sobpSequence , Nmaps] = flashAlgo_G4(sobpPosition{b}, sigmaAtIso , Plan.BDL ,  DRcritical{b} , plots_BEV{b}); %Find the best scan trajectory for beam b.
                    else
                      [sobpSequence , Nmaps] = flashAlgo_G3(sobpPosition{b}, sigmaAtIso , Plan.BDL ,  DRcritical{b} , Plan.Beams(b).NbScarves , plots_BEV{b}); %Find the best scan trajectory for beam b.
                    end

                  otherwise
                    SpotTrajectoryInfo.TimingMode
                    error('Undefined timing mode')
                end

                SpotTrajectoryInfo.beam{b}.Nmaps = Nmaps;
                SpotTrajectoryInfo.beam{b}.sobpSequence  = sobpSequence;

          end
  end

end
