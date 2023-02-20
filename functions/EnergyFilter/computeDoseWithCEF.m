%% computeDoseWithCEF
% Compute the dose distribution using MCsquare for a plan with a single beam and a single energy layer
% and a conformal energy modulator.
% The planning CT scan is re-interpolated on a pixel grid with the pixel size of the CEM in order to properly model the CEM.
% The dose is tallied on a grid with resolution |Plan.CEFDoseGrid|.
% The dose map is saved with the resolution |Plan.CEFDoseGrid|
%
% The intermediate high resolution dose map and high resolution CT scan are optionally saved.
% The final dose map is saved in a DICOM file in the folder |outputPath| and a DICOM plan is saved in the folder to which the dose map are referenced
%
%% Syntax
% |computeDoseWithCEF(Plan, outputPath, handles, PTV, CTName , FLAGdosePerSpot)|
%
%
%% Description
% |computeDoseWithCEF(Plan, outputPath, handles, PTV, CTName , FLAGdosePerSpot)| Description
%
%
%% Input arguments
%
% |Plan| -_STRUCTURE_- Structure defining a multi energy layer PBS plan
%
% |outputPath| -_STRING_- Path to the output folder where the DICOM file of the dose map and the mono-layer plan will be saved
%
% |handles| -_STRUCTURE_- REggui data handle. The CT scan is stored in |handles| in the image with name |Plan.CTname|.
%
% |PTV| - _SCALAR MATRIX_ - Mask defining the position of the PTV |PTV(x,y,z)=1| if the voxel is inside the PTV
%
% |CTname| -_STRING_- Name of the CT image in handles.images
%
% |FLAGdosePerSpot| -_BOOLEAN_- TRUE : compute high res dose per beamlet. FALSE : compute high res dose in whole volume
%
%
%% Output arguments
%
% None
%
%
%% Contributors
% Authors : R. Labarbe, Lucian Hotoiu (open.reggui@gmail.com)

function computeDoseWithCEF(Plan, outputPath, handles, PTV, CTName , FLAGdosePerSpot)


    global g_HUair;
    global g_HUcem; %Define HU as a global variable that will be visible inside the function getHighResDose
                  %getMaterialSPR reads the disk. Only make the reading once and not at every iteration of the loop for each beamlet in order to make computation faster
    g_HUair =  getMaterialSPR('Schneider_Air' , Plan.ScannerDirectory) + 1; %Hounsfield unit associated to air in the material file
    g_HUcem =  getMaterialSPR(Plan.Spike.MaterialID , Plan.ScannerDirectory) +1 ; %Hounsfield unit associated to CEM in the material file
        %The g_HUair and g_HUcem are now stored in a global varialbe that is accessible by the sub-function getHighResDose

    if numel(Plan.Beams) > 1
      error('computeDoseWithCEF requires plans with single beam')
    end
    if (numel(Plan.Beams.Layers) ~= 1)
      error('computeDoseWithCEF requires plans with single energy layer')
    end

    if ~isfield(Plan , 'SaveDoseBeamlets')
      %If the field is not defined, this means we do not want to save the dose map of the beamlets
      % in the reference frame of the CT scan
      Plan.SaveDoseBeamlets = false;
    end

    %Save the monolayer plan to the output folder
    %This will be used to reference to dose map to it
    if (~exist(fullfile(outputPath,'CEF_beam'),'dir'))
      %The folder to save the CT does not exist. Create it
      mkdir (fullfile(outputPath,'CEF_beam'))
    end
    Plan.FileName = 'Plan';
    regguiPath = fileparts(which('reggui'));
    dictionary = fullfile(regguiPath ,'plugins','openMIROpt','functions','io','dicom-dict.txt');
    createDICOMPlan(Plan,Plan.CTinfo,outputPath,dictionary)

    %Insert the range shifter in the high resolution CT
    %The insertion is done in the low resolution CT to avoid wasting time re-inserting the range shifter in all high reoslution beamlets
    fprintf('Adding Range shifter to low resolution CT \n')
    [Plan , handles ] = setRangeShifterinCT(handles , Plan , CTName);

     if ~FLAGdosePerSpot
        %The interpolated pixels are large. Compute the dose in the whole volume in one go
        [Plan , ~ , ~ , ~ , Zdistal ,minField , maxField, sDoseIECg] = getHRdoseGridInfo(Plan , handles , PTV);
        [DoseIECg , DoseFileName , handlesDoseIECg] = getHighResDose(Plan, outputPath , handles , PTV , CTName , 0 , minField , maxField , Zdistal , true);
        [iCTgntX , iCTgntY , iCTgntZ ] = getCTaxes(handlesDoseIECg.origin , handlesDoseIECg.spacing , handlesDoseIECg.size , [0,0,0]); %Cooridnates of dose map aligned with IEC gantry
        [DoseOrig , handlesDose]  = doseIECg2DICOMcs(DoseIECg , outputPath , handlesDoseIECg , handles , Plan, iCTgntX, iCTgntY , iCTgntZ); %Dose map aligned to orignal Ct with resolution of |Plan.CEFDoseGrid|

     else
        %The voxels are small. Compute the dose in one beamlet at a time
        fprintf('Computing high resolution dose map one spot at a time \n')
        [DoseOrig, DoseFileName , handlesDose] = getFullHighResDosemap(Plan , handles , PTV , outputPath, CTName);

    end

    %Save dose map in original grid
    planFullPath = fullfile(outputPath,Plan.FileName);
    handlesDose = save2Disk(handlesDose, DoseOrig , handlesDose.size , Plan.CTinfo , 'Dose_withCEF' , fullfile(outputPath,'CEF_beam') , planFullPath , 'RTDOSE');

    %Copy the output file to the 'Output' folder
    movefile (fullfile(outputPath,'CEF_beam', 'Dose_withCEF.dcm' ) , fullfile(outputPath,'Dose_withCEF.dcm'));

  end % This is the end of the computeDoseWithCEF


%-----------------------------------------------
%Define the parameters of the high resolution dose map
%-----------------------------------------------
function [Plan , iCTgntX , iCTgntY , iCTgntZ, Zdistal , minField , maxField, sDoseIECg] = getHRdoseGridInfo(Plan , handles , PTV)

  % Set coarse grid for dose calculation.
  Plan.Independent_scoring_grid = 1; % Enable a different scoring grid than the base CT for the dose calculation
  Plan.resampleScoringGrid = false; %Do not resample the dose map. Leave it on the coarse scoring grid
  Plan.Scoring_voxel_spacing = cell2mat(Plan.CEFDoseGrid); % In [mm]. Set dose calcuation scorinng grid to 1mm spacing and overwrite the CT grid. This would reduce computation time if the CT resolution is very high.
  format = ['Dose scoring grid: ' repmat(' %1.0f',1,numel(Plan.Scoring_voxel_spacing)) '\n'];
  fprintf(format,Plan.Scoring_voxel_spacing);

  %Create a dose map at the resolution |Plan.Scoring_voxel_spacing|
  [minField , maxField] = getMaxBEVsize(Plan.Beams);
  [iCTgntX , iCTgntY , iCTgntZ, Zdistal] =  getDoseMapCoordInIECg(Plan.Beams , PTV , handles.spacing , handles.origin , minField , maxField ,  Plan.Scoring_voxel_spacing);
  sDoseIECg = [numel(iCTgntX) , numel(iCTgntY) , numel(iCTgntZ)];

end


%==============================================
% Compute high resolution dose map in the IEC gantry CS, one beamlet at a time
% Do not intrepolate for each beamlet, just insert at the proper place in the map
% then only rotate at the end to place the dose map in the original CT orientation
%while keeping the dose map resolution
%
% The function return a dose map with a aspatial resolution of |Plan.Scoring_voxel_spacing|
% and aligned wit the coordinate axes of the original CT scan
%================================================
function [DoseOrigCT, DoseFileName , handlesDose] = getFullHighResDosemap(Plan , handles , PTV, outputPath , CTName)

  FieldSize = 40; % mm This radius must be larger than the spot radius at the distal surface of PTV. Otherwise the lattice structure will be visible in the dose map

  NbBeamlets = numel(Plan.Beams.Layers(1).SpotWeights); %Number of PBS spots in plan with energy monolayer

  % Set coarse grid for dose calculation.
  [Plan , iCTgntX , iCTgntY , iCTgntZ, Zdistal ] = getHRdoseGridInfo(Plan , handles , PTV);
  DoseIECg = zeros(numel(iCTgntX) , numel(iCTgntY) , numel(iCTgntZ) ); %Create an empty dose map in which the dose from each bemalet will be saved

  for spt = 1:NbBeamlets

      SptPos = Plan.Beams.Layers(1).SpotPositions(spt,:);

      %Define a field size that match the pixel position of the CEM
      %This will avoid the alising problems when accumulating the dose
      minField = alignFieldWithCEM(SptPos - FieldSize , Plan.Beams.RangeModulator);
      maxField = SptPos + FieldSize;
      Zdistal = getClosePixel(iCTgntZ , Zdistal);

      fprintf('Computing dose for PBS spots %d of %d at [%f mm, %f mm] \n',spt,NbBeamlets,SptPos(1),SptPos(2))
      [doseBeamlet, DoseFileName , handlesDoseIECg] = getHighResDose(Plan, outputPath, handles, PTV, CTName, spt, minField , maxField , Zdistal , false); %Add the dose of each spot to the global high resolution dose map

      %Align the parital dose map on the pixel of the full dose map
      [~ , IndexStart(1)] = getClosePixel(iCTgntX , handlesDoseIECg.origin(1));
      [~ , IndexStart(2)] = getClosePixel(iCTgntY , handlesDoseIECg.origin(2));
      [~ , IndexStart(3)] = getClosePixel(iCTgntZ , handlesDoseIECg.origin(3));

      DoseIECg(IndexStart(1):IndexStart(1)+handlesDoseIECg.size(1)-1 , IndexStart(2):IndexStart(2)+handlesDoseIECg.size(2)-1 , IndexStart(3):IndexStart(3)+handlesDoseIECg.size(3)-1) = ...
                         DoseIECg(IndexStart(1):IndexStart(1)+handlesDoseIECg.size(1)-1 , IndexStart(2):IndexStart(2)+handlesDoseIECg.size(2)-1 , IndexStart(3):IndexStart(3)+handlesDoseIECg.size(3)-1) + doseBeamlet;


      %If required, save the dose map of each beamlet in the reference frame of the orignial CT scan
      DosePath = fullfile(outputPath,'CEF_beam','Outputs','SpotDoseInCT');
      switch Plan.SaveDoseBeamlets
        case 'dcm'
          if (~exist(DosePath,'dir'))
            %The folder to save the CT does not exist. Create it
            mkdir (DosePath)
          end
          planFullPath = fullfile(outputPath,Plan.FileName);
          save2Disk(handlesDoseIECg, doseBeamlet , handlesDoseIECg.size , Plan.CTinfo , DoseFileName , DosePath , planFullPath, 'RTDOSE'); %Save the beamlet

        case 'sparse'
            if (~exist(DosePath ,'dir'))
              %The folder to save the CT does not exist. Create it
              mkdir (DosePath )
            end
            save(fullfile(DosePath , [DoseFileName '.mat']) , 'doseBeamlet')

         otherwise
            %Do not save anything
      end %switch Plan.SaveDoseBeamlets

  end %end for spt

  [DoseOrigCT, handlesDose] = doseIECg2DICOMcs(DoseIECg , outputPath , handlesDoseIECg , handles , Plan, iCTgntX, iCTgntY , iCTgntZ);


end

%------------------------------------------
%REsample the dose map with pixels aligned on the IEC gantry CS
% to a dose map aligned on the DICOM CS
%------------------------------------------
function [DoseOrigCT, handlesDose] = doseIECg2DICOMcs(DoseIECg , outputPath , handlesDoseIECg , handles , Plan, iCTgntX, iCTgntY , iCTgntZ)

  %Create a fake handle for the dose map in the original CT CS
  handlesDose = struct;
  handlesDose.path = outputPath;
  handlesDose = Initialize_reggui_handles(handlesDose); % Initialize the handles structure
  handlesDose.spacing = handlesDoseIECg.spacing;
  handlesDose.origin = handles.origin;
  handlesDose.size = ceil(handles.size .* handles.spacing ./ handlesDoseIECg.spacing); %Compute the number of pixels at dose map resolution to fill the volume of original CT

  %Compute the coordinate of all pixels of the full dose map
  %The DoseOrig dose map is aligned with the axes of the original DICOM CT
  [oCTdcmX, oCTdcmY , oCTdcmZ , X , Y , Z] = getCTaxes(handlesDose.origin , handlesDose.spacing , handlesDose.size , Plan.Beams.isocenter); %Create coordinate system of original CT, centered on original isocenter
  A = [X(:), Y(:), Z(:)];
  A = [A , ones(numel(X),1)]; %Coordinate of all voxels of full dose map in DICOM CS
  Beam2 = Plan.Beams;
  oCTgntr = DICOM2gantry(A' , Beam2); %Coordinate of all voxels of the original CT scan, expressed in IEC gantry. Rotate arouind isocenter at [0,0,0], using gantry angle =0, PPS=0

  %Interpolate the dose map into the original CT scan
  DoseOrigCT = interpolateCT(DoseIECg, iCTgntX, iCTgntY , iCTgntZ , oCTgntr(1,:) , oCTgntr(2,:) , oCTgntr(3,:) , 0 , handlesDose.size  , 'linear'); %The dose in the orignal grid

end


  %===========================================================
  % GEt the dose map at the resolution specified in |CEFDoseGrid| within the field defined by |minField| & |maxField|  and for the specified spot
  % The CT scan is interpolated at the same resolution as the CEM mask
  % The axes of the dose map are aligned with the IEC gantry CS
  %
  % INPUT
  %
  % |Plan| -_STRUCTURE_- Structure defining a multi energy layer PBS plan
  %
  % |outputPath| -_STRING_- Path to the folder where the DICOM file with the mono-layer plan will be saved
  %
  % |handles| -_STRUCTURE_- REggui data handle. The CT scan is stored in |handles| in the image with name |Plan.CTname|.
  %
  % |PTV| - _SCALAR MATRIX_ - Mask defining the position of the PTV |PTV(x,y,z)=1| if the voxel is inside the PTV
  %
  % |CTname| -_STRING_- Name of the CT image in handles.images
  %
  % |spt| -_SCALAR_- [OPTIONAL. Default: compute dose on whole CT] Index of the spot to process.
  %
  % |minField| -_SCALAR VECTOR_- [x,y,z] coordinate (mm) in the IEC gantry CS of the first corner of the field where the dose is to be computed. Must be aligned with a pixel of the CEM
  %
  % |maxField| -_SCALAR VECTOR_- [x,y,z] coordinate (mm) in the IEC gantry CS of the last corner of the field where the dose is to be computed.
  %
  % |Zdistal| -_SCLAR_-  Z Coordinate (mm) in the IEC gantry CS of the deepest plane in which the dose is to be computed
  %
  % |WholeField| -_BOOLEAN_- |true| if the dose is to be computed for the all bemalets. |false| if the dose is to be computed for one isngle beamlet
  %
  % OUTPUT
  %
  % |Dose| -_SCALAR MATRIX_- Dose(x,y,z) Dose (Gy) at the pixel [x,y,z]. Matrix axes are aligned with IEC gantry CS
  %
  % |DoseFileName| -_STRING_- File name of the dose map computed by MCsquare at the resolution |scoring_voxel_spacing| and in the IEC gantry CS
  %
  % |DoseHR| -_STRUCTURE_- Handles with the image properties of the dose map at the resolution |scoring_voxel_spacing| and in the IEC gantry CS
  %===========================================================

  function [Dose, DoseFileName , DoseHR] = getHighResDose(Plan, outputPath , handles , PTV , CTName , spt , minField , maxField , Zdistal , WholeField)


    global g_HUair;
    global g_HUcem; %Define HU as a global variable that will be visible inside the function getHighResDose

    %Define the name of the output dose file
    if WholeField
        DoseFileName = 'dose_HighRes';
    elseif (isfield(Plan, 'CEFDoseGrid') && all(cell2mat(Plan.CEFDoseGrid)))
        SptPos = Plan.Beams.Layers(1).SpotPositions(spt,:);
        DoseFileName = ['dose_UserGrid_spt' , num2str(spt) , '_spt_X_',num2str(round(SptPos(1))),'_Y_',num2str(round(SptPos(2)))];
    else
        SptPos = Plan.Beams.Layers(1).SpotPositions(spt,:);
        DoseFileName = ['dose_HighRes_spt' , num2str(spt) , '_spt_X_',num2str(round(SptPos(1))),'_Y_',num2str(round(SptPos(2)))];
    end
    DosePath = fullfile(outputPath,'CEF_beam','Outputs');
    DoseFullFileName = fullfile(DosePath,[DoseFileName,'.dcm']);

    %Insert the CEM into the high resolution CT scan
    % Define min and max field so that we do not expand the CT scan in the Xg and Yg direction to fit the CEM
    hrCTName = 'highResCTbev';

    %check that the field position is aligned with the pixel resolution and origin of CEM
    a = (minField - Plan.Beams.RangeModulator.ModulatorOrigin(1:2)) ./ Plan.Beams.RangeModulator.Modulator3DPixelSpacing(1:2);
    b = round( (minField - Plan.Beams.RangeModulator.ModulatorOrigin(1:2)) ./ Plan.Beams.RangeModulator.Modulator3DPixelSpacing(1:2) );
    if (max(abs(a - b)) > 1e-4)
      a
      b
      a-b
      error('Small field not aligned with CEM grid')
    end

    %Generate the high resolution CT scan
    %The IEC gantry is aligned with the Y axis of the CT scan
    %the spatial resolution of the Ct is the same as CEM to avoid aliasing problems when inserting CEM in high res CT
    fprintf('Interpolating high resolution CT \n')
    [handlesHR , BeamHR ] = createHighResCT(handles , CTName , hrCTName , PTV, Plan.Beams , Plan.Beams.RangeModulator.Modulator3DPixelSpacing , g_HUair , minField , maxField , Zdistal , Plan.CTinfo);
    PlanHR.Beams = BeamHR;
    PlanHR = copyFields(PlanHR , Plan); %Overwrite the YAML parameter in the reloaded |Plan|
    PlanHR.CTname =  hrCTName;
    PlanHR  = updatePlanCTparam(handlesHR , PlanHR );

    %Add the CEM into the high resolution CT
    fprintf('Adding CEM to high resolution CT\n')
    [PlanHR , handlesHR ] = setCEMinCT(handlesHR , PlanHR , hrCTName , minField , maxField, g_HUcem , g_HUair);

    %Save the interpolated CT on disk
    %The IEC gantry axis of the beamlet is aligned with the Y axis of this high resolution Ct scan
    if (~exist(fullfile(outputPath,'CEF_beam'),'dir') && Plan.SaveHighResCT)
      %The folder to save the CT does not exist. Create it
      mkdir (fullfile(outputPath,'CEF_beam'))
    end
    if (~exist(DoseFullFileName) && Plan.SaveHighResCT)
      %The dose file does not already exist. Save the high resolution CT if we were asked to save it
      CTintrp = Get_reggui_data(handlesHR, hrCTName ,'images'); %Get the high res CT
      infoCTinterp = PlanHR.CTinfo;
      infoCTinterp.ImagePositionPatient = handlesHR.origin;
      infoCTinterp.Spacing = handlesHR.spacing;
      save_Image(CTintrp , infoCTinterp , fullfile(outputPath,'CEF_beam','CT_with_CEF'),'dcm');
      CTintrp = []; %Clean memory
    else
      %The output folder exists. Skip saving the CT
    end

    %Save the monolayer plan to disk
    PlanMono = saveMonoLayerPlan(Plan , spt ,  handlesHR.isocenter , handlesHR.origin , handlesHR.spacing , WholeField  , outputPath);

    if ~WholeField
      %We are working with one spot at a time.
      %Save the folder with the partial high resolution CT scan
      folderName = fullfile(outputPath,'CEF_beam',['CT_with_CEF_spt_X_',num2str(round(SptPos(1),1)),'_Y_',num2str(round(SptPos(2),1))]);
      if(~exist(DoseFullFileName)  && Plan.SaveHighResCT)
        %The dose file does not already exist. Make a backup of the high resolution CT scan into a folder
        % If the dose file already exist, then the backup of the CT scan already exist and we do not want to overwrite it
        copyfile (fullfile(outputPath,'CEF_beam','CT_with_CEF') , folderName);
      end
    end

    %Prepare data for monolayer plan for MCsquare computation
    Plan2 = struct;
    Plan2.CTname = hrCTName;
    Plan2.output_path = fullfile(outputPath,'CEF_beam');
    Plan2.ScannerDirectory = Plan.ScannerDirectory;
    Plan2.BDL = Plan.BDL;
    Plan2.protonsFullDose = Plan.protonsHighResDose;
    Plan2.CTinfo = PlanMono.CTinfo;

    %Compute the dose using MCsquare on the high resolution CT scan
    if (isfield(Plan, 'CEFDoseGrid') && all(cell2mat(Plan.CEFDoseGrid)))
        %Define the MCsquare parameter to do the dose scoring on a different pixel size than the high resolution CT scan
        fprintf('----Dose scoring grid---- \n')
        Plan2.Independent_scoring_grid = Plan.Independent_scoring_grid; % Enable a different scoring grid than the base CT for the dose calculation
        Plan2.resampleScoringGrid = Plan.resampleScoringGrid;
        Plan2.Scoring_voxel_spacing = Plan.Scoring_voxel_spacing; % In [mm]. Set dose calcuation scorinng grid to 1mm spacing and overwrite the CT grid. This would reduce computation time if the CT resolution is very high.
    end

    % Compute the dose map if no dose file already exists
    if (~exist(DoseFullFileName))
      %The output folder does not exist. Compute the high resolution dose
      if Plan.SaveHighResDoseMap
        %The user wants the high resolution dose in IEC gantry at dcm format
        %ComputeFinalDose uses the plan info form the fille Plan.dcm that is saved in the folder Plan2.output_path
        handlesHR = ComputeFinalDose(Plan2, handlesHR , DoseFileName);
      else
        %The user does not want to save the high resolution dose map in IEC gantry
        handlesHR = ComputeFinalDose(Plan2, handlesHR , []);
      end
      [Dose , DoseInfo ] = Get_reggui_data(handlesHR,'dose_final_miropt');
                      %Get the image from handlesHR.myData. This is the idx-th element with |handlesHR.myData.name{idx} = 'dose_final_miropt'|
                      %This is an image at the resolution of the scoring grid defined in |CEFDoseGrid|
                      %It has not been resampled at the resolution of the HR CT scan.
    else
      %The output folder exists. Skip dose computation. Simply relaod the file
      fprintf('File %s .dcm already exists. Skipping dose computation \n' , DoseFileName)
      [Dose , DoseInfo ] = load_Image(DosePath,DoseFileName,'dcm',0);
    end

    %Rotate the dose map to be aligned with the IEC gantry again
    Dose = flip(Dose,2); %NB: flip (C code) is faster than flipdim (.m).
    Dose = permute(Dose,[1,3,2]);

    %Define the image properties for the dose map in the IEC gantry CS
    DoseHR = struct;
    DoseHR.path = outputPath;
    DoseHR = Initialize_reggui_handles(DoseHR);
    DoseHR.size = size(Dose);
    origZ = -DoseInfo.ImagePositionPatient(2) - DoseInfo.Spacing(2) .* (DoseHR.size(3) - 1);


    DoseHR.origin = [DoseInfo.ImagePositionPatient(1) , DoseInfo.ImagePositionPatient(3) , origZ ]'; %Second index is minus (because flipped) Zg
    DoseHR.isocenter = [handlesHR.isocenter(1) , handlesHR.isocenter(3) , handlesHR.isocenter(2)]; %permuted the CT scan
    DoseHR.spacing = [ DoseInfo.Spacing(1) , DoseInfo.Spacing(3) , DoseInfo.Spacing(2)]'; %permutation has no effect here: all elements are identical

  end


%-----------------------------------
% Compute the maximum BEV field size (IEC gantry) in which the dose should be computed
% The field size is chosen to be aligned with the pixels of the CEM mask
% in order to avoid aliasing problems when inserting CEM into high reoslution CT
%------------------------------------
function [minField , maxField] = getMaxBEVsize(Beam)

      %There is no aperture block. The field size is defined by CEF size
      [~ , ~ , ~ , BlockBorder] = findApertureBlockSize({}); %Retrieve the border size from this function
      BlockBorder = rounding([BlockBorder,BlockBorder] , Beam.RangeModulator.Modulator3DPixelSpacing(1:2)); %Round the border to a multiple of the number of voxels

      NbPxlCEF = size(Beam.RangeModulator.CEM3Dmask); %[Nx,Ny,Nz] number of pixels along X and Y IEC gantry
      SizeCEF = Beam.RangeModulator.Modulator3DPixelSpacing .* (NbPxlCEF - 1); %[Sx,Sy] (mm) dimension of the CEF block along X and Y IEC gantry
      minCEF = Beam.RangeModulator.ModulatorOrigin; %One point of the main diagonal of the CEF
      maxCEF = Beam.RangeModulator.ModulatorOrigin + SizeCEF; %Other point of the main diagonal of the CEF

      minField = minCEF(1:2) - BlockBorder; %Add the border which is a multiple of the number of voxels
      maxField = maxCEF(1:2) + BlockBorder;

end

%========================
% Convert coordinates of points from the DICOM CS to the IEC gantry CS
%========================
function Adcm = DICOM2gantry(A , Beam)
    M = matDICOM2IECgantry(Beam.GantryAngle,Beam.PatientSupportAngle);
    Adcm = M * A ;
    %The rotation is done orond the origin of the CS. The isocenter is at origin. No need to worry about the + isocenter. Make computation faster
    %Coordinate of voxels in IEC gantry
end


%===========================
% Save the mono layer plan to use with the high resolution CT
%
% The monolayer plan is linked to the high resolution Ct which assumes gantry angle 0°
% The modifications made to |Plan| are local to this function and will not affect the calling function
%=========================
function PlanMono = saveMonoLayerPlan(Plan , spt , isocenterINTER , ImagePositionPatientINTER , spacingINTER , WholeField , outputPath)

  fprintf('Saving DICOM plan with monolayer \n')

  %The interpolated CT was rotated so that the CEM axis is paralell to the Ydicom axis
  %We will use a virtualm gantry at 0° and no couch rotation to shout through the CEM in the interpolated CT
  %TODO This asusmes that the BDL does not depend on the gantry angle !!!!!!!!!!!!!
  Plan.Beams.GantryAngle = 0;
  Plan.Beams.PatientSupportAngle = 0;
  Plan.Beams.isocenter = isocenterINTER; %Redefine the position of isocentre in the resampled CT

  %Remove the range shifter definition. It is now included inside the CT scan
  Plan.Beams.RSinfo.RangeShifterSetting ='OUT';
  Plan.Beams = rmfield(Plan.Beams,'RSinfo');
  Plan.Beams.NumberOfRangeShifters = 0;

  PlanMono = Plan;
  PlanMono.CTinfo.ImagePositionPatient = ImagePositionPatientINTER;
  PlanMono.CTinfo.Spacing = spacingINTER;

  if ~WholeField
    %Compute the dose for one single spot only
    PlanMono.Beams.Layers(1).SpotPositions = PlanMono.Beams.Layers(1).SpotPositions(spt,:);
    PlanMono.Beams.Layers(1).SpotWeights = PlanMono.Beams.Layers(1).SpotWeights(spt);
  end

  regguiPath = fileparts(which('reggui'));
  dictionary = fullfile(regguiPath ,'plugins','openMIROpt','functions','io','dicom-dict.txt');
  createDICOMPlan(PlanMono,PlanMono.CTinfo,fullfile(outputPath,'CEF_beam'),dictionary)

end

%================================================
% Get the pixel coordinates of the dose map in IEC CS
% The first pixel of the dose map is algined with the pixels of the CEM
%================================================
function [iCTgntX , iCTgntY , iCTgntZ , Zdistal] =  getDoseMapCoordInIECg(Beam , PTV , Spacing , ImagePositionPatient , minField , maxField ,  PixelSizeIECg)

  DepthExtension = 50; % mm Extend the computation to this depth beyond PTV distal surface

  %Get the maximum Zg extension of the CEM
  %This will defined one of the maximum extension of the interpolated CT scan
  switch Beam.RangeModulator.ModulatorMountingPosition
    case 'SOURCE_SIDE'
      % The CEM is pointing towards the source.
      %Add the CEM height to the postion of the base
      NbPxlCEF = size(Beam.RangeModulator.CEM3Dmask); %[Nx,Ny,Nz] number of pixels along X and Y IEC gantry
      SizeCEF = Beam.RangeModulator.Modulator3DPixelSpacing .* NbPxlCEF; %[Sx,Sy,Sz] (mm) dimension of the CEM
      maxCEF = Beam.RangeModulator.IsocenterToRangeModulatorDistance + SizeCEF(3);

    case  'PATIENT_SIDE'
      %The base of the CEM is at the maximum Zg
      maxCEF = Beam.RangeModulator.IsocenterToRangeModulatorDistance;
  end

  Bdcm = getDICOMcoord(PTV , Spacing , ImagePositionPatient , Beam.isocenter);
  M = matDICOM2IECgantry(Beam.GantryAngle,Beam.PatientSupportAngle);
  Bbev = M * Bdcm; %Coordinate of the voxels of the body in IC gantry
  Zdistal = min(Bbev(3,:)) - DepthExtension; %Position on beam axis of the distal surface of PTV + additional depth

  %Define the coordinate of the voxels in the IEC gantry CS for the interpolated Ct scan
  %Make the dose map a bit larger to be sire to fit all bemalets
  iCTgntX = (minField(1) - 15 .* PixelSizeIECg(1) )      : PixelSizeIECg(1) : (maxField(1) + 15 .* PixelSizeIECg(1));
  iCTgntY = (minField(2) - 15 .* PixelSizeIECg(2) )      : PixelSizeIECg(2) : (maxField(2) + 15 .* PixelSizeIECg(2));
  iCTgntZ = double( (Zdistal - 15 .* PixelSizeIECg(3) )  : PixelSizeIECg(3) : maxCEF + 30 .* PixelSizeIECg(3));

end


%-------------------------------------------
% find the closest pixel coordinate to a given number
%-------------------------------------------
function [pxlVal , pixIdx] = getClosePixel(iCTgntX , value)
  [~ , pixIdx] = min((iCTgntX - value).^2);
  pxlVal = iCTgntX(pixIdx);
end

%----------------------------------------------
% Align the field with the origin of the CEM and the CEM pixel resolution
%----------------------------------------------
function minField = alignFieldWithCEM(minField , RangeModulator)
  NbPxl = round( (minField - RangeModulator.ModulatorOrigin(1:2)) ./ RangeModulator.Modulator3DPixelSpacing(1:2)); %Number of pixels between origin of CEM and corner of the field
  minField = RangeModulator.ModulatorOrigin(1:2) + NbPxl .* RangeModulator.Modulator3DPixelSpacing(1:2); %Shift the origiin of the field to be aligned with CEM grid
end
