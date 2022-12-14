%% flashLoadAndCompute
% Load a DICOM RT plan (ConformalFLASH plan) and:
%  * Optimize spot trajectory if  |BeamProp.FLAGOptimiseSpotOrder = true|
%  * Compute the dose rate in the |RTstruct.selected_ROIs|
%
%
%% Syntax
% |[handles, Plan] = flashLoadAndCompute(planFileName, CTname , rtstructFileName , output_path , BeamProp , RTstruct, DoseRate , CEMprop)|
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
% |DoseRate| -_STRUCTURE_- [OPTIONAL: If empty, does not compute the dose rate histograms] Definition of parameters for the computation of dose rates
%   * |DoseRate.Dref| -_SCALAR_- Threshold dose (Gy / FRACTION) in OAR above which the dose rate condition must be respected. Voxels below the threshold dose are not included in DR condition
%   * |DoseRate.DMF| -_SCALAR_- Dose modifying factor of flash for this organ. Not used for dose constraint.
%   * |DoseRate.DR50| -_SCALAR_- Gy/s 50% dose rate of the sigmoid probability function. Not used for dose constraint.
%
% |CEMprop| -_STRUCTURE_- Information for the computation of the CEM
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

function [handles, Plan] = flashLoadAndCompute(planFileName, CTname , rtstructFileName , output_path , BeamProp , RTstruct, DoseRate , CEMprop , scanAlgoGW , spots)

if nargin < 9
  scanAlgoGW = [];
end

if nargin < 10
  spots = [];
end


%Add MIROpt specific configuration
%-------------
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

if CEMprop.makeSTL
  %Export the STL file of the CEM
  path2beamResults = getOutputDir(Plan.output_path , 1);
  filename = fullfile(path2beamResults,[matlab.lang.makeValidName(Plan.Beams.RangeModulator.AccessoryCode),'.stl']);
  exportCEM2STL(Plan.Beams.RangeModulator.CEM3Dmask  , Plan.Beams.RangeModulator.Modulator3DPixelSpacing , Plan.Beams.RangeModulator.ModulatorOrigin , Plan.Beams.RangeModulator.AccessoryCode , filename)
end


%Load the RT structures
%----------------------
[handles , Plan , ROI, usedROI] = loadRSrtStructs(rtstructFileName , handles, Plan, RTstruct, DoseRate);
Plan  = updatePlanCTparam(handles , Plan  );


% Add the aperture in the CT scan
%--------------------------------
fprintf('Adding aperture to high resolution CT\n')
[Plan , handles ] = setApertureInCT(handles , Plan , Plan.CTname); %Add an apertrue block in the CT scan
[handles , Plan , ROI] = updateROI(handles , Plan , ROI); %Update the ROI mask used to load the dose influence matrices

%Compute the optimum spot trajectory
%----------------------------------
if BeamProp.FLAGOptimiseSpotOrder
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
      Plan.SpotTrajectoryInfo.beam{b}.TimePerSpot = Plan.Beams(b).Layers.duration; %|TimePerSpot(s)| Duration (ms) of the s-th spot
      deltaT = diff(Plan.Beams(b).Layers.time); %time difference between spots
      deltaT = deltaT - Plan.Beams(b).Layers.duration(2:end);
      Plan.SpotTrajectoryInfo.beam{b}.dT = [0; deltaT]; %|dT(st)| Sweep (ms) to move from the s-1-th spot to the s-th spot
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

  %Make sure that the spot trajectory defined in the plan is the same as MIROPT
  if BeamProp.FLAGcheckSpotOrdering
        fprintf('Checking that trajectory of original plan is the same than optimised with MIROPT \n')
        SpotTrajectoryInfo = optimizeTrajectory(Plan  , ROI); %Optimsie the trajectory with simple model or scanAlgo

        if (~prod(diff(SpotTrajectoryInfo.beam{1}.sobpSequence)==1))
            %If the spot trajectory is not a monotonically increasing sequence, then the spot trajectory of the plan
            %was optimised with a different algorithm than MIROPT (which was supposed to be used in the RaysStation script)
            %This is an error condition
            Plan.Beams(1).Layers(1).SpotPositions
            Plan.SpotTrajectoryInfo.sobpPosition{1}
            SpotTrajectoryInfo.weight2spot
            error('Spot ordering has been changed from the order in the plan')
        end
    end

end

PlanMono = Plan; %This is a Monolayer plan already

%Compute dose map
%-----------------
% Compute the dose through the CEF using the high resolution CT scan
%Save the high resolution dose map of each beamlet in separate files
fprintf('Computing the dose map in high resolution CT scan \n')
for b = 1:numel(Plan.Beams)
  fprintf('Beam %d \n' , b)
  %Find the PTV in the structure set
  for i = 1:length(ROI)
      if (strcmp(ROI(i).name, Plan.TargetROI))
          idxPTV = i;
      end
  end

  %Compute the dose of each beamlet
  path2beamResults = getOutputDir(Plan.output_path , b);
  computeDoseWithCEF(Plan , path2beamResults , handles , ROI(idxPTV).mask3D.value , 'CTwithAperture' , true);
  movefile (fullfile(Plan.output_path,'Outputs','Plan.dcm') , fullfile(path2beamResults,'Plan_CEM.dcm'));
end


%Compute the dose influence matrix
%---------------------------------
[handles, Plan] = getRS_doses(Plan, handles);

%Compute the dose rate and save the results to disk
%--------------------------------------------------
%Include all structure in the dose rate computation
for  optFidx = 1:numel(Plan.optFunction)
    Plan.optFunction(optFidx).ID = 8; %Activate the dose rate computation
end

%Compute dose rate in all structures
handles = ComputeFinalDoseRate(Plan, handles, ROI);


end
