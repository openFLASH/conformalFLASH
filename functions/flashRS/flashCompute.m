%% flashCompute
% Receive CT, RTStructs and DICOM RT plan (ConformalFLASH plan) and:
%  * Optimize spot trajectory if  |BeamProp.FLAGOptimiseSpotOrder = true|
%  * Compute the dose rate in the |RTstruct.selected_ROIs|
%
%
%% Syntax
% |[handles, Plan] = flashCompute(handles, output_path, BeamProp, RTstruct, DoseRate, CEMprop)|
%
%
%% Description
% |[handles, Plan] = flashCompute(handles, output_path, BeamProp, RTstruct, DoseRate, CEMprop)| Description
%
%
%% Input arguments
%
% |handles| -_STRUCT_- structure containing flash plan, CT scan and RTstructs
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
% Authors : L. Hotoiu (open.reggui@gmail.com)

function [handles, Plan, doseRatesCreated] = flashCompute(handles, SimuParam, BeamProp, RTstruct, DoseRate, CEMprop, scanAlgoGW, spots)

spots

    if nargin < 7
      scanAlgoGW = [];
    end

    if nargin < 8
      spots = [];
    end

    warning('on') %Turn on warning messages


    %Add MIROpt specific configuration
    %----------------------------------------------------------------------
    Plan = configMiropt_RS(BeamProp, CEMprop, SimuParam.OutputDirectory);


    %Load the CT scan and update the data
    %----------------------------------------------------------------------
    [Plan, MCsqExecPath, ScannersPath] = getScannerCalib(handles, SimuParam.CT, Plan);
    %[handles, Plan] =  loadRSctScan(CTname , [], Plan);


    %Load the plan
    %----------------------------------------------------------------------
    if ~isempty(scanAlgoGW)
      %Define the parameters for the conextion to scnaAlgo gateway
      Plan.scanAlgoGW.scanalgoGateway_IP = scanAlgoGW.scanalgoGateway_IP;
      Plan.scanAlgoGW.room_id = scanAlgoGW.room_id;
      Plan.scanAlgoGW.snout_id = scanAlgoGW.snout_id; %TODO get this infor from plan
      Plan.scanAlgoGW.spot_id = scanAlgoGW.spot_id; %TODO get this infor from plan
    end

    [handles, Plan] = parseFLASHplan(SimuParam.RefPlanName, Plan, handles);

    %If records from logs are provided, then overwrite the spot info in plan
    % with the log records
    if ~isempty(spots)
      fprintf('Overwriting spot info with record from logs \n')
      Plan = overwriteWithLogs(Plan, spots);
    end

    if (~exist(fullfile(Plan.output_path,'Outputs'),'dir'))
      mkdir (fullfile(Plan.output_path,'Outputs'))
    end
    copyfile (SimuParam.RefPlanName, fullfile(Plan.output_path,'Outputs','Plan.dcm')); %Copy the RS plan into the output folder to be used when creating the dose maps

    if CEMprop.makeSTL
      %Export the STL file of the CEM
      path2beamResults = getOutputDir(Plan.output_path , 1);
      filename = fullfile(path2beamResults,[matlab.lang.makeValidName(Plan.Beams.RangeModulator.AccessoryCode),'.stl']);
      exportCEM2STL(Plan.Beams.RangeModulator.CEMThicknessData, Plan.Beams.RangeModulator.Modulator3DPixelSpacing, Plan.Beams.RangeModulator.ModulatorOrigin, Plan.Beams.RangeModulator.AccessoryCode, Plan.Beams.RangeModulator.ModulatorMountingPosition , filename)
    end


    %Load the RT structures
    %----------------------------------------------------------------------
    %[handles , Plan , ROI, usedROI] = loadRSrtStructs(rtstructFileName , handles, Plan, RTstruct, DoseRate);
    [handles, Plan, ROI, usedROI] = make_DR_opt_func(handles, Plan, RTstruct, DoseRate);
    Plan = updatePlanCTparam(handles, Plan);


    %Compute the optimum spot trajectory
    %----------------------------------------------------------------------
    if Plan.FLAGOptimiseSpotOrder
      %Optimise spot trajectory
      % If |Plan.scanAlgoGW| is defined, then use scanAlgo. In the output |SpotTrajectoryInfo.scanAlgoGW| will be defined
      % Otherwise, use simple model
      fprintf('Optimising spot trajectory \n')
      [SpotTrajectoryInfo , Plan] = optimizeTrajectory(Plan, ROI);
      Plan.SpotTrajectoryInfo = SpotTrajectoryInfo;

    else
      %Do not optimise spot trajectory. Use the one from the plan
      fprintf('Using spot trajectory from plan\n')

      [sobpPosition , weight2spot ] =  collectSpotsinBeamlets(Plan, [])
                  %Gather the spots into beamlets for form SOBP
                  %This is a mono layer plan, so we should get back the same spots
      Plan.SpotTrajectoryInfo.sobpPosition =  sobpPosition;
      Plan.SpotTrajectoryInfo.weight2spot = weight2spot;

      for b = 1:numel(Plan.Beams)
        Plan.SpotTrajectoryInfo.beam{b}.sobpSequence  = weight2spot(:,2);

Plan.Beams(b).Layers

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

Plan.SpotTrajectoryInfo
pause

        if Plan.showGraph
          figure(100+b)
          hold on
          plot(Plan.SpotTrajectoryInfo.sobpPosition{b}(:,1) , Plan.SpotTrajectoryInfo.sobpPosition{b}(:,2) , '-k')
        end
      end
    end


    %Compute dose map
    %----------------------------------------------------------------------
    % Compute the dose through the CEF using the high resolution CT scan
    %Save the high resolution dose map of each beamlet in separate files
    fprintf('Computing the dose map in high resolution CT scan \n')

    for b = 1:numel(Plan.Beams)
      %TODO deal with plan containing setup beams
      fprintf('Beam %d \n' , b)

      %Compute the dose of each beamlet
      path2beamResults = getOutputDir(Plan.output_path , b);
      Plan.Scenario4D(1).RandomScenario(Plan.rr_nominal).RangeScenario(Plan.rs_nominal).P = computeDoseWithCEF(Plan, path2beamResults , handles , Plan.CTname, true);
      movefile (fullfile(Plan.output_path,'Outputs','Plan.dcm') , fullfile(path2beamResults,'Plan_CEM.dcm'));
    end


    %Include all structure in the dose rate computation
    %----------------------------------------------------------------------
    for  optFidx = 1:numel(Plan.optFunction)
        Plan.optFunction(optFidx).ID = 8; %Activate the dose rate computation
    end

    %Compute dose rate in all structures and save the results to disk
    [handles, doseRatesCreated] = ComputeFinalDoseRate(Plan, handles, ROI);

end
%##########################################################################



function [Plan, MCsqExecPath, ScannersPath] = getScannerCalib(handles, CTimageName, Plan)

    [~,info] = Get_reggui_data(handles, CTimageName);

    %Read the scanner model from the plan in order to identify the correct scanner calibration folder for MCsquare
    if isfield(info.OriginalHeader , 'ManufacturerModelName')
      Plan.ScannerDirectory = remove_bad_chars(info.OriginalHeader.ManufacturerModelName);
    else
      warning('(008,1090) ManufacturerModelName undefined in CT scan')
      Plan.ScannerDirectory = 'default';
      fprintf('Using default folder for CT scanner calibration : %s \n', Plan.ScannerDirectory);
    end
    fprintf('CT scanner calibration file : %s \n', Plan.ScannerDirectory);

    %Check that the scanner directory exists.
    [pluginPath , MCsqExecPath , BDLpath , MaterialsPath , ScannersPath] = get_MCsquare_folders();
    if (~exist(Plan.ScannerDirectory,'dir'))
      %TODO This should be an error and the program should stop here in a well configured system
      warning('(008,1090) ManufacturerModelName references an unknown scanner')
      fprintf('Folder for CT scanner calibration : %s \n', Plan.ScannerDirectory);
      Plan.ScannerDirectory = 'default';
      fprintf('Using default folder for CT scanner calibration : %s \n', Plan.ScannerDirectory);
    end

    Plan.CTinfo = info;
    Plan.CTname = CTimageName;

    Plan.CTinfo.Spacing = info.Spacing;
    Plan.DoseGrid.resolution = Plan.CTinfo.Spacing;
    Plan.CTinfo.ImagePositionPatient = info.ImagePositionPatient;
end
%##########################################################################



function [handles, Plan, ROI, usedROI] = make_DR_opt_func(handles, Plan, RTstruct, DoseRate)

    strName = unique({RTstruct.selected_ROIs{:}, RTstruct.TargetROI, RTstruct.ExternalROI}); %Make sure there is no duplicate of name so that we load the structure only once

    % Create ROI specific properties for all structures
    Plan.ExternalROI = [remove_bad_chars(RTstruct.ExternalROI)];
    Plan.TargetROI = [remove_bad_chars(RTstruct.TargetROI)];
    Plan.TargetROI_ID = numel(strName) - 1;
    Plan.ExternalROI_ID = numel(strName);

    Plan.DoseGrid.size = handles.size';
    %nvoxels = prod(handles.size);
    %Plan.OptROIVoxels_nominal = false(nvoxels,1); % initialize to zeros
    %Plan.OptROIVoxels_robust = false(nvoxels,1); % initialize to zeros
    %ROItotalMask = zeros(Plan.DoseGrid.size); %Create a structure which is the sum of all ROI

    %Create the ROI structure
    for i = 1:numel(strName)
        % Set ROI name and create list of ROIs
        RTstruct.selected_ROIs{i} = [remove_bad_chars(strName{i})];
        ROImask = Get_reggui_data(handles, RTstruct.selected_ROIs{i});
        ROI(i) = createROI(RTstruct.selected_ROIs{i}, ROImask);

        if (sum(strcmp(strName{i} , RTstruct.DRCritical_ROIs)))
            %This is a dose rate critical structure to include in the spot trajectory
            Plan.optFunction(i) = createOptFunctionStruct(remove_bad_chars(strName{i}), ROI, DoseRate , true);
          else
            %This is not a dose rate critical structure. Ignore it for the trajectory computation
            Plan.optFunction(i) = createOptFunctionStruct(remove_bad_chars(strName{i}), ROI, DoseRate , false);
        end

        %if(i<numel(strName)) %Do not add EXTERNAL to the mask
        %    Plan.OptROIVoxels_nominal = Plan.OptROIVoxels_nominal | ROI(i).mask1D;
        %    ROItotalMask = ROItotalMask  | ROImask;
        %end
    end

    %Add the mask for the total volume in which the dose rate is to be computed
    %TotalMaskName = [remove_bad_chars('TotalMask')];
    %RTstruct.selected_ROIs{end+1} = TotalMaskName;
    %handles = Set_reggui_data(handles,TotalMaskName,ROItotalMask , [] , 'images');
    %ROImask = Get_reggui_data(handles, TotalMaskName);
    %ROI(end+1) = createROI(TotalMaskName, ROImask);
    %Plan.optFunction(end+1) = createOptFunctionStruct(remove_bad_chars('TotalMask') , ROI , DoseRate , false);

    Plan.TargetROI_ID = getROIByName(ROI, [remove_bad_chars(RTstruct.TargetROI)]);

    %Get the list of structure names
    usedROI = RTstruct.selected_ROIs;
end
%##########################################################################



% Create a structure defining the dose rate objective function
function optFunction = createOptFunctionStruct(ROIName , ROI , DoseRate , isInTrajectory)

  if ~isInTrajectory
    optFunction.ID = 1; %This is not an objective on dose rate. Ignore from trajectory computation
  else
    optFunction.ID = 8; %Select this ID because the optimisation function ID=8 requires the computation of dose rate
  end

  optFunction.ROIindex = getROIByName(ROI, ROIName);
  optFunction.ROIname = ROIName;

  if isfield(DoseRate,'Dref')
    optFunction.Dref  = DoseRate.Dref;
  end

  if isfield(DoseRate,'DMF')
    optFunction.DMF  = DoseRate.DMF;
  else
    optFunction.DMF = [];
  end
  if isfield(DoseRate,'DR50')
    optFunction.DR50  = DoseRate.DR50;
  else
    optFunction.DR50 = [];
  end
end
%##########################################################################



% Add the ROI to the list of ROIs
function ROI = createROI(ROIname, ROImask)

    ROI.name = ROIname;
    % Calculate mask for each ROI
    ROI.mask3D.value = ROImask;

    % Convert to 1D sparse mask with z inverse format (to be consistent
    % with beamlets format)
    temp = flip(ROI.mask3D.value,3);
    ROI.mask1D = sparse(logical(double(temp(:))));
end
%##########################################################################
