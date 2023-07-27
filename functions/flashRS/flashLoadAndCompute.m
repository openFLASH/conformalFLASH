%% flashLoadAndCompute
% Load a DICOM RT plan (ConformalFLASH plan) and:
%  * Optimize spot trajectory if  |BeamProp.FLAGOptimiseSpotOrder = true|
%  * Compute the dose rate in the |RTstruct.selected_ROIs|
%
%
%% Syntax
%
% |[handles, Plan, doseRatesCreated] = flashLoadAndCompute(  planFileName , CTname , rtstructFileName , output_path , BeamProp , ExternalROI , CEMprop , scanAlgoGW ,  spots)|
%
% |[handles, Plan, doseRatesCreated] = flashLoadAndCompute(  planFileName , CTname , []               , output_path , BeamProp , ExternalROI , CEMprop , scanAlgoGW ,  spots  , handles)|
%
%
%% Description
% |[handles, Plan] = flashLoadAndCompute(planFileName, CTname , rtstructFileName , output_path , BeamProp , RTstruct, DoseRate , CEMprop)| Description
%
%
%% Input arguments
%
% |planFileName| -_STRING_- Full path and file name to the RT DICOM plan
%
% |CTname| -_STRING_- If handles = [] or missing , file name and full path to the CT scan
%                     If handles is provided, name of the CT scan in |handles|
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
% |ExternalROI| -_STRING_- Name of the RT struct with the body contour
%
% |CEMprop| -_STRUCTURE_- [OPTIONAL. If absent, no CEM is present in the plan] Information for the computation of the CEM
%    * |makeSTL| -_BOOL_- If true, the STL file is saved in the output folder
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
% |handles| -_STRUCT_- [OPTIONAL] REGGUI image handles. If provided, the CT scan and the Rt struct are read from handles. Otherwise, they are loaded  from disk
%
%% Output arguments
%
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are stored.
%
% |handles| - _STRUCT_ - REGGUI data structure.
%
% |doseRatesCreated| -_CELL VECTOR of STRINGS_- List of names of dose rate maps saved to disk and in |handles|
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [handles, Plan, doseRatesCreated] = flashLoadAndCompute( planFileName , CTname , rtstructFileName , output_path , BeamProp , ExternalROI , CEMprop , scanAlgoGW ,  spots  , handles)

if nargin < 7
  CEMprop = struct;
end

if nargin < 8
  scanAlgoGW = [];
end

if nargin < 9
  spots = [];
end

warning('on') %Turn on warning messages


%Add MIROpt specific configuration
%----------------------------------------------------------------------
Plan = configMiropt_RS(BeamProp, CEMprop, output_path);

%Load the CT scan and update the data
%----------------------------------------------------------------------
if nargin < 10 || isempty(handles)
  %hanldes are not provided. Load CT from disk
  % CTname is a file name
  [handles, Plan] =  loadRSctScan(CTname , []      , Plan);
else
  %handles is provided. Do not load CT from disk
  %CTname is the image name in |handles|
  Plan.CTname = CTname;
  [handles, Plan] =  loadRSctScan([]    , handles , Plan);
end

%Load the plan
%----------------------------------------------------------------------
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

if CEMprop.makeSTL
  %Export the STL file of the CEM
  path2beamResults = getOutputDir(Plan.output_path , 1);
  filename = fullfile(path2beamResults,[matlab.lang.makeValidName(Plan.Beams.RangeModulator.AccessoryCode),'.stl']);
  exportCEM2STL(Plan.Beams.RangeModulator.CEMThicknessData  , Plan.Beams.RangeModulator.Modulator3DPixelSpacing , Plan.Beams.RangeModulator.ModulatorOrigin , Plan.Beams.RangeModulator.AccessoryCode , Plan.Beams.RangeModulator.ModulatorMountingPosition , filename)
end


%Load the RT structures
%----------------------
[handles , Plan ] = loadEmptyStructs( rtstructFileName, handles, Plan, ExternalROI );
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

%Compute dose map
%-----------------
% Compute the dose through the CEF using the high resolution CT scan
%Save the high resolution dose map of each beamlet in separate files
fprintf('Computing the dose map in high resolution CT scan \n')
for b = 1:numel(Plan.Beams)
  %TODO deal with plan contianing setup beams
  fprintf('Beam %d \n' , b)

  %Compute the dose of each beamlet
  path2beamResults = getOutputDir(Plan.output_path , b);
  Plan = computeDoseWithCEF(Plan , path2beamResults , handles , Plan.CTname , true);
  movefile (fullfile(Plan.output_path,'Outputs','Plan.dcm') , fullfile(path2beamResults,'Plan_CEM.dcm'));
end

%Compute the dose rate and save the results to disk
%--------------------------------------------------

%Compute dose rate in all structures
[handles, doseRatesCreated] = ComputeFinalDoseRate(Plan, handles, ROI);


end
