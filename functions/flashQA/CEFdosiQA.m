%% Contributors
% Authors : Lucian Hotoiu (open.reggui@gmail.com)
% -------------------------------------------------------------------------

tic
clearvars
close all


%Inputs / Outputs

% CEF images
scanCEF_path = 'C:\Users\lhotoiu\Downloads\20230628\D58_NO_bubbles\CT.1.3.12.2.1107.5.1.4.83552.30000023062822093699000003556.dcm';
%refCEF_path = 'C:\Users\lhotoiu\Downloads\ctCEF505030_DCM\CEF505030_90degrees\reggui_CEF_origin_rot3_0001.dcm';

% CEM image names in reggui handles
scan_CEF_imageName = 'scan_CEF';
%ref_CEF_imageName = 'ref_CEF';

% Patient CT, Plan, RTstructs
path_patientCT = 'D:\MATLAB\Data\Tests\Raysearch\Dshape\CT_air\CT_air_0001.dcm';
path_rtstructs = 'D:\MATLAB\Data\Tests\Raysearch\Dshape\CT_air\RS1.2.752.243.1.1.20230221164349509.2960.70223.dcm';
path_RS_plan = 'D:\MATLAB\Data\Tests\Raysearch\Dshape\D_shape_matrix_measurements_14092023\Plan\FP-D58.dcm';

% Output path
outputPath = 'D:\MATLAB\Data\Tests\CEFQA\OUTPUT_dosi_test';

% Get structs
RTstruct.selected_ROIs = {'WaterCube'}; %Name of the RT structs for which the gamma index is to be computed
RTstruct.ExternalROI = 'WaterCube'; %name for external ROI - the body contour
RTstruct.TargetROI = 'D58'; %name for target ROI
% -------------------------------------------------------------------------

%Config code flags
%saveFilesToDisk = false;
%b = 1; % which beam we run the QA for
% -------------------------------------------------------------------------


% Initialise reggui handles
handles = Initialize_reggui_handles();
handles.dataPath = outputPath;
% -------------------------------------------------------------------------


% Initialise reggui handles to store hi-res CEM scan CT 
cts_handles = Initialize_reggui_handles();
cts_handles.dataPath = outputPath;
% -------------------------------------------------------------------------



% Config Plan variable to run dose with cef in MIROpt
Plan = configPlan(outputPath);

% Get the pre-calculated patient plan, CT, RTstructs
% Import plan
%-------------------------------------
[~, Plan] = parseFLASHplan(path_RS_plan, Plan, handles);   

%Load the scan and ref image (if any)
%-------------------------------------
cts_handles = loadCTDataSets(cts_handles, scanCEF_path, scan_CEF_imageName);

%Load the CT scan and update the data
%-------------------------------------
[handles, Plan] =  loadRSctScan(path_patientCT, handles, Plan);

%Load the RT structures
%-------------------------------------
[handles, Plan ] = loadEmptyStructs(path_rtstructs, handles, Plan, RTstruct.ExternalROI ); % this is called to precondition computeDOseWithCEF
% -------------------------------------------------------------------------



% Initialise reggui handles
%trans_handles = Initialize_reggui_handles();
%trans_handles.dataPath = outputPath;
% -------------------------------------------------------------------------



% Load data from scan and reference sources
% trans_handles = loadDataSets(scanCEF_path, scan_CEF_imageName, refCEF_path, ref_CEF_imageName, ref_source);

%[scan_cef_data, scan_cef_info, ~] = Get_reggui_data(trans_handles, scan_CEF_imageName);
%[ref_cef_data, ref_cef_info, ~] = Get_reggui_data(trans_handles, ref_CEF_imageName);

% % Upsample reference image with nearest-neighbor interpolation to match the
% % scan resolution. We want to use the reference image as the fixed image
% % for registration to have aligned voxels to the IEC axes.
% % Resample reference image to same spacing as scan - nearest interpolation to avoid smoothing
% 
% % Resample opt image to same spacing as scan - nearest interpolation to avoid smoothing
% cur_sz = size(ref_cef_data);
% target_sz = round(cur_sz .* (ref_cef_info.Spacing./scan_cef_info.Spacing)');
% [X, Y, Z] = meshgrid(linspace(1,cur_sz(2),target_sz(2)),linspace(1,cur_sz(1),target_sz(1)),linspace(1,cur_sz(3),target_sz(3)));
% X = single(X);
% Y = single(Y);
% Z = single(Z);
% ref_cef_data_resample = interp3(ref_cef_data,X,Y,Z,'nearest',-1024);
% ref_cef_info = scan_cef_info;
% ref_cef_resampled_imageName = 'ref_cef_resampled';
% trans_handles = Set_reggui_data(trans_handles, ref_cef_resampled_imageName, ref_cef_data_resample, ref_cef_info, 'images', 1); %Add the interpolated CT scan to the list of CT scan
% 
% 
% 
% % Rigid registration - reference fixed - moving scan cef
% fprintf('Computing rigid registration between reference (fixed) and scan (moving) images \n');
% scan_rigid_def_imageName = 'scan_rigid_def';
% scan_rigid_transName = 'scan_rigid_trans';
% 
% 
% % Ensure same coordinate system inside handles.origin as specified in the
% % mono plan, to force the registration and autothreshold to work with that origin system
% trans_handles.origin = Plan.Beams(b).CEMorigin;
% trans_handles = Registration_ITK_rigid_multimodal(ref_cef_resampled_imageName, scan_CEF_imageName, scan_rigid_def_imageName, scan_rigid_transName, trans_handles);
% [scan_rigid_def_data,~,~] = Get_reggui_data(trans_handles, scan_rigid_def_imageName,'images');
% % -------------------------------------------------------------------------
% 
% 
% % Get scan mask
% fprintf('Computing scan and reference binary mask images \n');
% 
% scan_cef_mask = 'scan_cef_mask';
% trans_handles = AutoThreshold(scan_rigid_def_imageName, [128], scan_cef_mask, trans_handles);
% [scan_cef_mask_data, scan_cef_mask_info, ~] = Get_reggui_data(trans_handles, scan_cef_mask);
% 
% ref_cef_mask = 'ref_cef_mask';
% trans_handles = AutoThreshold(ref_cef_resampled_imageName, [128], ref_cef_mask, trans_handles);
% [ref_cef_mask_data, ref_cef_mask_info, ~] = Get_reggui_data(trans_handles, ref_cef_mask);
% % -------------------------------------------------------------------------
% 
% 
% % Save registration files to disk
% if saveFilesToDisk
%     trans_handles = save2Disk(trans_handles, scan_rigid_def_data, size(scan_rigid_def_data), scan_cef_info, scan_rigid_def_imageName, fullfile(trans_handles.dataPath));
%     trans_handles = save2Disk(trans_handles, ref_cef_data_resample, size(ref_cef_data_resample), ref_cef_info, ref_cef_resampled_imageName, fullfile(trans_handles.dataPath));
%     trans_handles = save2Disk(trans_handles, scan_cef_mask_data, size(scan_cef_mask_data), scan_cef_mask_info, scan_cef_mask, fullfile(trans_handles.dataPath));
%     trans_handles = save2Disk(trans_handles, ref_cef_mask_data, size(ref_cef_mask_data), ref_cef_mask_info, ref_cef_mask, fullfile(trans_handles.dataPath));
% end
% -------------------------------------------------------------------------


% Copy CEF images as data into main handles containing the original patient CT
% This will include the CEF CT scan and the patient CT scan at different
% resolutions and origins into the same handles that are passed to the functions of
% "computeDoseWithCEF"
% -------------------------------------------------------------------------
%handles.mydata.name = trans_handles.images.name;
%handles.mydata.data = trans_handles.images.data;
%handles.mydata.info = trans_handles.images.info;



%Compute the dose through the reference CEM
%-------------------------------------
% Compute the dose through the CEM using the high resolution CT scan
fprintf('Computing dose map in high resolution CT using reference CEM \n')
for b = 1:numel(Plan.Beams)
  fprintf('Beam %d \n' , b)
   
  %Compute the dose
  path2beamResults_ref = getOutputDir(fullfile(Plan.output_path, 'Ref'), b);
  Plan = computeDoseWithCEF(Plan, path2beamResults_ref, handles, Plan.CTname, true);
end
% -------------------------------------------------------------------------


% Compute the dose through the scan CEM
%-------------------------------------
% Compute the dose through the CEM using the high resolution CT scan
fprintf('Computing dose map in high resolution CT using 3D-printed CEM scan \n')
for b = 1:numel(Plan.Beams)
  fprintf('Beam %d \n' , b)
  
  % Replace reference CEM with mask of 3D_printed scan  
  Plan = addCEMScanToPlan(cts_handles, Plan, scan_CEF_imageName, b);

  %Compute the dose
  path2beamResults_scan = getOutputDir(fullfile(Plan.output_path, 'Scan'), b);
  Plan = computeDoseWithCEF(Plan, path2beamResults_scan, handles, Plan.CTname, true);
end
% -------------------------------------------------------------------------


% Reimport dose maps
scan_CEM_doseName = 'Scan_CEM_dose';
scan_CEM_dose_onDisk_Name = 'Dose_withCEF_beam1.dcm';
handles = Import_image(path2beamResults_scan, scan_CEM_dose_onDisk_Name, 1, scan_CEM_doseName, handles);
[scan_CEM_doseData, scan_CEM_doseInfo,~] = Get_reggui_data(handles, scan_CEM_doseName,'images');

ref_CEM_doseName = 'Ref_CEM_dose';
ref_CEM_dose_onDisk_Name = 'Dose_withCEF_beam1.dcm';
handles = Import_image(path2beamResults_ref, ref_CEM_dose_onDisk_Name, 1, ref_CEM_doseName, handles);
[ref_CEM_doseData, ref_CEM_doseInfo,~] = Get_reggui_data(handles, ref_CEM_doseName,'images');


% Import contours to compute gamma index
strName = unique({RTstruct.selected_ROIs{:}, RTstruct.TargetROI, RTstruct.ExternalROI}); %Make sure there is no duplicate of name so that we load the structure only once
handles = Import_contour(path_rtstructs, strName, Plan.CTname, 1, handles);

% Compare scan dose to reference dose via gamma index
options.dd = 3; % (mm)
options.DD = 3; % (%)
options.FI = 5; % number of points for internal interpolation
options.global_ref = 1; % local (0) or global (1 - default) WET difference will be used
options.threshold = 5; % (% unit) specify a WET threshold (in % of the reference WET) under which the gamma is not computed
%options.search_distance = 5; % specify a reference WET to be used for global WET computation and thresholding instead of the max WET
fprintf('Computing dose gamma index (ref - scan) \n');

gamma_name = 'gamma_ref_scan';
ROI_name = 'ct_D58';
disp(['Computing gamma index for dose through scanned vs. reference CEM for ', ROI_name]);
[handles, gamma, passing_rate, average_gamma] = Gamma_index(ref_CEM_doseName, scan_CEM_doseName, ROI_name, gamma_name, handles, options);
handles = save2Disk(handles, gamma, size(gamma), ref_CEM_doseInfo, gamma_name, fullfile(Plan.output_path));
% -------------------------------------------------------------------------

toc



function Plan = addCEMScanToPlan(handles, Plan, scanImageName, b)
    
    scan_cef_mask = 'scan_cef_mask';
    handles = AutoThreshold(scanImageName, [128], scan_cef_mask, handles);
    [scan_cef_mask_data, scan_cef_mask_info, ~] = Get_reggui_data(handles, scan_cef_mask);

    % Rotate scan CEM to align with reference
    permOrder1 = [1,3,2]; %Permutation 1 of the dimensions of CT scan to align it with plan
    permOrder2 = [2,1,3]; %Permutation 2 of the dimensions of CT scan to align it with plan
    permOrder = [permOrder1; permOrder2];
    flipAxis = 2; %Which axis index should be flipped (after permute)

    % Permute the dimension of CEM image and flip some dimension
    % in order to get the smae orientation for the reference CT and the meausrmenet CT
    if flipAxis
        scan_cef_mask_data = flip(scan_cef_mask_data, flipAxis);
    end
    if ~isempty(permOrder)
        for i = 1:size(permOrder, 1)
            scan_cef_mask_data = permute(scan_cef_mask_data , permOrder(i,:));
            mask_spacing = scan_cef_mask_info.Spacing(permOrder(i,:));
        end
    end
    origin = (- round(size(handles.images.data{2}) ./2) .* mask_spacing')'; %Move origin back to middle of CEM   

    Plan.Beams(b).RangeModulator.CEM3Dmask = scan_cef_mask_data;
    Plan.Beams(b).RangeModulator.Modulator3DPixelSpacing = mask_spacing;
    Plan.Beams(b).RangeModulator.ModulatorOrigin = origin;
end
%--------------------------------------------------------------------------


% function Plan = addCEFtoPlan(Plan, b, CEFbinaryMask, CEF_info, RTstruct)
%     pxlSizeZ = CEF_info.Spacing(3); % Different in Z direction
%     pxlSizeXY = [CEF_info.Spacing(1), CEF_info.Spacing(2)]; % Same in X and Y directions
%     CEFmaxThickness = max(max(sum(CEFbinaryMask, 3, 'omitnan').* pxlSizeZ));
%     CompensatorPixelSpacing = round([pxlSizeXY, pxlSizeZ], 2); % 2 - keep only two decimals from pixel spacing value
%     ModulatorOrigin = CEF_info.ImagePositionPatient';
% 
% 
%     Plan.Beams(b).Iso2Skin = getIsocentreToSkinDistance(Plan.Beams(b), RTstruct.ExternalROI, 5 , Plan.DoseGrid.resolution , Plan.CTinfo.ImagePositionPatient) + 2; %Isocentre (mm) To Skin Distance along the proton beam axis
%     Plan.Beams(b).CompensatorThicknessData = CEFmaxThickness; %Compute the maximum height of the CEF. This will be used to estimate the isocentre to Range modulator distance
%     Plan.Beams(b).IsocenterToRangeModulatorDistance = getIsocenterToRangeModulatorDistance(Plan.Beams(b), Plan.BDL); %First guess at CEf thickness data
%     iso2Nozzle = getNozzle2IsoDistanceFromBDL(Plan.BDL);
%     if(Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance > iso2Nozzle)
%         fprintf('Isocenter to Modulator Tray Distance must be %f mm \n',Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance)
%         fprintf('Maximum distance between isocentre and nozzle frame (in BDL) is %f mm \n',iso2Nozzle)
%         fprintf('Select a BDL (= a treatment room) with more space between isocentre and nozzle frame \n')
%         fprintf('Alternatively, clear some space by moving the isocentre using the parameter isocenter in YAML \n')
%         error('CEM does not fit on the nozzle')
%     end
% 
%     %Store CEF data into the structure that will be exported in DICOM plan
%     Plan.Beams(b).CompensatorCEFmask = CEFbinaryMask;
%     Plan.Beams(b).CompensatorPixelSpacing = CompensatorPixelSpacing;
%     Plan.Beams(b).CompensatorOrigin = CompensatorOrigin;
%     Plan.Beams(b).IsocenterToRangeModulatorDistance = getIsocenterToRangeModulatorDistance(Plan.Beams(b), Plan.BDL); %Compute the distance to base of CEF now that we know the exact thickness of CEF
% end
%--------------------------------------------------------------------------


function Plan = configPlan(output_path)

    BeamProp.protonsHighResDose = 1e2; % Number of protons in the dose in high resolution CT    
    BeamProp.NbScarves = 1; % Number of scarves to paint on the BEV
    BeamProp.CEFDoseGrid = {3, 3, 3}; % Size (mm) of final dose scoring grid. Compute the final dose through CEF on a different grid than the high-res
    BeamProp.FLAGOptimiseSpotOrder = false;
    BeamProp.FLAGcheckSpotOrdering = false;
    BeamProp.BDL = 'D:\MATLAB\flash\openMCsquare\lib\BDL\BDL_default_UN1_G0_Al_RangeShifter_tilted.txt';
    BeamProp.ScannerDirectory = 'D:\MATLAB\REGGUI\plugins\openMCsquare\lib\Scanners\default';
    BeamProp.MCsqExecPath = 'D:\MATLAB\REGGUI\plugins\openMCsquare\lib';

    CEMprop.makeSTL = false;

    Plan = configMiropt_RS(BeamProp, CEMprop, output_path);

    Plan.showGraph = true;  % true = display the graphs during computation
    Plan.SaveDoseBeamlets = 'dcm';  % Save the dose of each beamlet in the reference frame of the CT with aperture: dcm (DICOM format) , sparse (sparse matrix) , false (not saved)
    Plan.SaveHighResCT = true;         % Do not save the high resolution CT for each beamlet in the reference frame of the beamlet
    Plan.SaveHighResDoseMap = false;   % Do not                     save the dose map at CEFDoseGrid resolution in the reference frame of the beamlet

end
%--------------------------------------------------------------------------


function [handles] = loadCTDataSets(handles, scan_CEM_path, scan_CEF_imageName, ref_CEM_path, ref_CEF_imageName)
    
    if nargin > 3
        % Import ref CEF
        [ref_CEF_directory, ref_CEF_fileName, EXT] = fileparts(ref_CEM_path);
        handles = Import_data(ref_CEF_directory, [ref_CEF_fileName EXT], 1, ref_CEF_imageName, handles);  
    end

    % Import scanned CEF
    [scan_CEF_directory, scan_CEF_fileName, EXT] = fileparts(scan_CEM_path);
    handles = Import_image(scan_CEF_directory, [scan_CEF_fileName EXT], 1, scan_CEF_imageName, handles);  
end
