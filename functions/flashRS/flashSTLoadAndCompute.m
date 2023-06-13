%% flashSTLoadAndCompute
% Load a DICOM RT plan (shoot through plan) and:
%  * Compute the dose map
%  * Compute the dose rate
%
%
%% Syntax
% |[handles, Plan] = flashSTLoadAndCompute(planFileName, CTname , rtstructFileName , output_path , BeamProp , RTstruct, DoseRate , CEMprop)|
%
%
%% Description
% |[handles, Plan] = flashSTLoadAndCompute(planFileName, CTname , rtstructFileName , output_path , BeamProp , RTstruct, DoseRate , CEMprop)| Description
%
%
%% Input arguments
%
% |planFileName| -_STRING_- Full path and file name to the RT DICOM plan
%
% |CTname| -_STRING_- File name and full path to the CT scan
%
% |rtstructFileName| -_STRING_- File name and full path to the DICOM RT struct
%
% |output_path| -_STRING_- Folder in which the results are saved
%
% |BeamProp| -_STRUCTURE_- INformation about beam properties
%   *  |BeamProp.Inozzle| -_SCALAR_- Nozzle average current (nA)
%   *  |BeamProp.fractions| -_SCALAR_- Number of fractions for the treatment.
%   *  |BeamProp.BDL| -_STRING_- Beam data library. Name of the folder in REGGUI\plugins\openMCsquare\lib\BDL
%   *  |BeamProp.ScannerDirectory| - _STRING_ - Name of the folder containing the definition of the CT scanner properties in MCsquare in folder "plugins\openMCsquare\lib\Scanners"
%   *  |BeamProp.MachineType| - _STRING_ - Description of the treatment machine (PROTEUSone , PROTEUSplus)
%   *  |BeamProp.FLAGOptimiseSpotOrder| -_BOOL_- TRUE = Optimise the spot order in the plan. FALSE = use spot order form plan but check that this is the same spot order than MIROPT
%   *  |BeamProp.FLAGcheckSpotOrdering| -_BOOL_- TRUE = Check that spot ordering in plan matches scanAlgo output
%
% |RTstruct| -_STRUCTURE_- Information about the RT structs used in the TPS
%   * |RTstruct.selected_ROIs| -_CELL VECTOR_- |RTstruct.selected_ROIs{i}| is a stirng with the name of the i-th RT struct in which the dose rate is to be computed
%   * |RTstruct.ExternalROI| -_STRING_- Name of the RT struct with the body contour
%   * |RTstruct.TargetROI| -_STRING_- Name of the RT struct with the PTV
%
% |scanAlgoGW| -_STRUCTURE_- [OPTIONAL : default : no conexion to scanAlgo] Information about the scanAlgo gateway
%    * |scanAlgoGW.scanalgoGateway_IP| -_STRING_- IP address, including port, to the scanAlkgo gatewat
%    * |scanAlgoGW.room_id| -_STRING_- Room ID as defined  inthe gateway
%    * |scanAlgoGW.snout_id|  -_STRING_- snout ID as defined in the gateway
%    * |scanAlgoGW.spot_id| -_STRING_-  beam supply point as defined in the site config jar site.properties
%
% |spots| -_STRUCTURE_- [OPTIONAL. If absent, the spot information is obtained from plan]
%     * |spots.name| -_STRING_- Name of the beam to which the log record refers
%     * |spots.spots.xy| - _SCALAR VECTOR_ - Average spot position (x,y) at isocenter over the delivery of the s-th spot
%     * |spots.spots.weight| - _SCALAR_ - Monitor unit of the s-th spot
%     * |spots.spots.time| - _SCALAR_ - Time (in ms) at the begining of the delivery of the s-th spot
%     * |spots.spots.duration| - _SCALAR_ - Time (in ms) at the end of the delivery of the s-th spot of the l-th energy layer
%
%
%% Output arguments
%
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are stored.
%
% |handles| - _STRUCT_ - REGGUI data structure.
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [handles, Plan] = flashSTLoadAndCompute(planFileName, CTname , rtstructFileName , output_path , BeamProp , RTstruct, scanAlgoGW , spots)



if nargin < 7
  scanAlgoGW = [];
end

if nargin < 8
  spots = [];
end

warning('on') %Turn on warning messages


%Add MIROpt specific configuration
%-------------
CEMprop = struct;
Plan = configMiropt_RS(BeamProp, CEMprop, output_path);

%Load the CT scan and update the data
%-------------------------------------
[handles, Plan] =  loadRSctScan(CTname , [], Plan);

%Load the plan
%-------------
if ~isempty(scanAlgoGW)
  %Define the parameters for the conextion to scnaAlgo gateway
  Plan.scanAlgoGW.scanalgoGateway_IP = scanAlgoGW.scanalgoGateway_IP;
  Plan.scanAlgoGW.room_id = scanAlgoGW.room_id;
  Plan.scanAlgoGW.snout_id = scanAlgoGW.snout_id; %TODO get this infor from plan
  Plan.scanAlgoGW.spot_id = scanAlgoGW.spot_id; %TODO get this infor from plan
end


[handles, Plan] = parseFLASHplan(planFileName , Plan, handles);

%If records from logs are provided, then overwrite the spot info in plan
% with the log records
if ~isempty(spots)
  fprintf('Overwriting spot info with record from logs \n')
  Plan = overwriteWithLogs(Plan, spots);
end

if (~exist(fullfile(Plan.output_path,'Outputs'),'dir'))
  mkdir (fullfile(Plan.output_path,'Outputs'))
end
copyfile (planFileName, fullfile(Plan.output_path,'Outputs','Plan.dcm')); %Copy the RS plan into the output folder to be used when creating the dose maps


%Load the RT structures
%----------------------
[handles , Plan ] = loadEmptyStructs( rtstructFileName, handles, Plan, RTstruct );
ROI = [];
Plan  = updatePlanCTparam(handles , Plan  );

%Compute the optimum spot trajectory
%----------------------------------
if Plan.FLAGOptimiseSpotOrder
  %Optimise spot trajectory
  % If |Plan.scanAlgoGW| is defined, then use scanAlgo. In the output |SpotTrajectoryInfo.scanAlgoGW| will be defined
  % Otherwise, use simple model
  fprintf('Optimising spot trajectory \n')
  [SpotTrajectoryInfo , Plan] = optimizeTrajectory(Plan  , ROI);
  Plan.SpotTrajectoryInfo = SpotTrajectoryInfo;

else
  %Do not optimise spot trajectory. Use the one from the plan
  fprintf('Using spot trajectory from plan\n')

  [sobpPosition , weight2spot ] =  collectSpotsinBeamlets(Plan , []);
              %Gather the spots into beamlets for form SOBP
              %This is a mono layer plan, so we should get back the same spots
  Plan.SpotTrajectoryInfo.sobpPosition =  sobpPosition;
  Plan.SpotTrajectoryInfo.weight2spot = weight2spot;

  for b = 1:numel(Plan.Beams)
    Plan.SpotTrajectoryInfo.beam{b}.sobpSequence  = weight2spot(:,2);

    if isfield(Plan.Beams(b).Layers ,'time') && isfield(Plan.Beams(b).Layers ,'duration')
      %Spot timing provided from log. We will use the logs
      Plan.SpotTrajectoryInfo.beam{b}.TimePerSpot = Plan.Beams(b).Layers.duration'; %|TimePerSpot(s)| Duration (ms) of the s-th spot
      deltaT = diff(Plan.Beams(b).Layers.time); %time difference between spots
      deltaT = deltaT - Plan.Beams(b).Layers.duration(2:end);
      Plan.SpotTrajectoryInfo.beam{b}.dT = deltaT; %|dT(st)| Sweep (ms) to move from the s-1-th spot to the s-th spot. dT has one less element than |Plan.SpotTrajectoryInfo.beam{b}.TimePerSpot| because it does not have the sweep time of the first spot
      Plan.SpotTrajectoryInfo.TimingMode = 'Record'; %The spot timing is recovered from logs

    else
      %No spot timing provided. We will compute the trajectory using the simple model or scanAlgo
      Plan.SpotTrajectoryInfo.beam{b}.Nmaps = getTopologicalMaps(sobpPosition{b} , Plan.BDL , Plan.Beams(b).spotSigma(b) , scanAlgoGW); %Get the inital topological map;
      Plan.SpotTrajectoryInfo.TimingMode = 'Model'; %The spot timing is from a model
    end

    if Plan.showGraph
      figure(100+b)
      hold on
      plot(Plan.SpotTrajectoryInfo.sobpPosition{b}(:,1) , Plan.SpotTrajectoryInfo.sobpPosition{b}(:,2) , '-k')
    end
  end
end

PlanMono = Plan; %This is a Monolayer plan already

%Compute dose influence matix
fprintf('Inserting aperture in CT \n')
[Plan , handles ] = setApertureInCT(handles , Plan , Plan.CTname , Plan.CTname); %Add an apertrue block in the CT scan

for b = 1:numel(Plan.Beams)
  %TODO deal with plan contianing setup beams
  fprintf('Beam %d \n' , b)
  fprintf('Inserting range shifter in CT \n')
  [Plan , handles] = setRangeShifterinCT(handles , Plan , Plan.CTname);
  Plan.Beams(b).NumberOfRangeShifters = 0;  %Remove the range shifter from the MCsquare beam model. The range shifter is now inserted in the CT scan
end

%Update the size of the masks after expanding the CT scan
Plan  = updatePlanCTparam(handles, Plan);
nvoxels = prod(handles.size);

handles.size
Body = Get_reggui_data(handles , Plan.ExternalROI);

temp = flip(Body,3);
%Plan.OptROIVoxels_nominal = sparse(logical(double(temp(:)))); %Mask defining the voxels of the dose influence matrix to be loaded in the sparse matrix
Plan.OptROIVoxels_nominal = true(nvoxels,1);
%Plan.OptROIVoxels_robust = false(nvoxels,1); % initialize to zeros

Export_image(Plan.CTname,fullfile(Plan.output_path,'Outputs','ct_IMPT'),'dcm',handles);


%Compute the dose influence matrix. It will be used to compute the dose rate
%There is no hedgehog in the plan, so we can use the standard MIROPT function to compute the Piuj matrix
[warm_start_in,Plan] = ComputePijsMatrix(Plan, handles, Plan.CTname); %Compute the influence matrix
Plan.PlanExistentFile = Plan.output_path;


OptConfig.BeamletsMatrixPrecision = 'd';
Plan = ReadPijsMatrix(Plan, OptConfig, handles);

%Compute the dose rate and save the results to disk
path2beamResults = getOutputDir(Plan.output_path , 1);
copyfile (planFileName, fullfile(path2beamResults,'Plan.dcm')); %Copy the RS plan into the output folder to be used when creating the dose maps
handles = ComputeFinalDoseRate(Plan, handles, ROI);

%Compute the final dose with higher number of protons
DoseFileName = 'MCsquare_Dose.dcm';
handles = ComputeFinalDose(Plan, handles, DoseFileName, 'multi-energy');

end
