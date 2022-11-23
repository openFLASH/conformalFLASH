%% computeDoseWithCEF
% Compute the dose distribution using MCsquare using a plan with a single beam and a single energy layer
% and a conformal energy modulator.
% The planning CT scan is re-interpolated on a pixel grid with small pixel size.
% The dose is tallied on a grid with resolution |Plan.CEFDoseGrid|.
% The dose map is then re-interpolated on the grid of the planning CT scan
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


    %Create a fake handle to save the dose maps
    handles2 = struct;
    handles2.path = outputPath;
    handles2 = Initialize_reggui_handles(handles2); % Initialize the handles structure
    handles2.spacing = handles.spacing;
    handles2.origin = handles.origin;

     if ~FLAGdosePerSpot
        %The interpolated pixels are large. Compute the dose in the whole volume in one go
        DoseOrig = getHighResDose(Plan, outputPath, handles , PTV, CTName);
     else
        %The voxels are small. Compute the dose in one beamlet at a time
        fprintf('Computing high resolution dose map one spot at a time \n')
        PlanMono = Plan;
        NbBeamlets = numel(PlanMono.Beams.Layers(1).SpotWeights); %Number of PBS spots in plan with energy monolayer

        sCT = handles.size;
        DoseOrig = zeros(sCT(1),sCT(2),sCT(3)); %Create an empty dose map in which the dose from each bemalet will be saved

        for spt = 1:NbBeamlets
            SptPos = PlanMono.Beams.Layers(1).SpotPositions(spt,:);
            fprintf('Computing dose for PBS spots %d of %d at [%f mm, %f mm] \n',spt,NbBeamlets,SptPos(1),SptPos(2))

            if (isfield(Plan, 'CEFDoseGrid') && all(cell2mat(Plan.CEFDoseGrid)))
                fprintf('Computing final dose with a coarser grid on high-res CT of PBS spot \n');
                % Set coarse grid for dose calculation.
                Plan.Independent_scoring_grid = 1; % Enable a different scoring grid than the base CT for the dose calculation
                Plan.Scoring_voxel_spacing = cell2mat(Plan.CEFDoseGrid); % In [mm]. Set dose calcuation scorinng grid to 1mm spacing and overwrite the CT grid. This would reduce computation time if the CT resolution is very high.
                format = ['Dose scoring grid: ' repmat(' %1.0f',1,numel(Plan.Scoring_voxel_spacing)) '\n'];
                fprintf(format,Plan.Scoring_voxel_spacing);
            end

            [doseBeamlet, DoseFileName]= getHighResDose(Plan, outputPath, handles, PTV, CTName, spt, SptPos); %Add the dose of each spot to the global high resolution dose map
            DoseOrig = DoseOrig + doseBeamlet;

            %If required, save the dose map of each beamlet in the reference frame of the orignial CT scan
            DosePath = fullfile(outputPath,'CEF_beam','Outputs','SpotDoseInCT');
            switch Plan.SaveDoseBeamlets
              case 'dcm'
                if (~exist(DosePath,'dir'))
                  %The folder to save the CT does not exist. Create it
                  mkdir (DosePath)
                end
                planFullPath = fullfile(outputPath,Plan.FileName);
                handles2 = save2Disk(handles2, doseBeamlet , Plan.DoseGrid.size , Plan.CTinfo , DoseFileName , DosePath , planFullPath, 'RTDOSE'); %Save the beamlet

              case 'sparse'
                  if (~exist(DosePath ,'dir'))
                    %The folder to save the CT does not exist. Create it
                    mkdir (DosePath )
                  end
                  save(fullfile(DosePath , [DoseFileName '.mat']) , 'doseBeamlet')

               otherwise
                  %Do not save anything
            end %switch Plan.SaveDoseBeamlets

        end
    end
    %Save dose map in original grid
    planFullPath = fullfile(outputPath,Plan.FileName);
    handles2 = save2Disk(handles2, DoseOrig , Plan.DoseGrid.size , Plan.CTinfo , 'Dose_withCEF' , fullfile(outputPath,'CEF_beam') , planFullPath , 'RTDOSE');

    %Copy the output file to the 'Output' folder
    movefile (fullfile(outputPath,'CEF_beam', 'Dose_withCEF.dcm' ) , fullfile(outputPath,'Dose_withCEF.dcm'));
end

%===========================================================
% GEt the high resolution dose within the requested field size and for the specified spot
%
% INPUT
%
% |Plan| -_STRUCTURE_- Structure defining a multi energy layer PBS plan
%
% |CTname| -_STRING_- Name of the CT image in handles.images
%
% |outputPath| -_STRING_- Path to the folder where the DICOM file with the mono-layer plan will be saved
%
% |handles| -_STRUCTURE_- REggui data handle. The CT scan is stored in |handles| in the image with name |Plan.CTname|.
%
% |PTV| - _SCALAR MATRIX_ - Mask defining the position of the PTV |PTV(x,y,z)=1| if the voxel is inside the PTV
%
% |SptPos| -_SCALAR VECTOR_- [x,y] coordinate (mm) of the single PBS spot for which the dose map is to be computed
%
% OUTPUT
%
% |DoseOrig| -_SCALAR MATRIX_- DoseOrig(x,y,z) Dose map in the same matrix size as the input CT scan
%
% |DoseFileName| -_STRING_- File name of the high resolution dose map saved on disk
%===========================================================

function [DoseOrig, DoseFileName] = getHighResDose(Plan, outputPath , handles , PTV , CTName, spt , SptPos)

  if nargin < 7
    WholeField = true;
    spt = 0;
  else
    WholeField = false;
  end

  %Define the name of the output dose file
  if WholeField
      DoseFileName = 'dose_HighRes';
  elseif (isfield(Plan, 'CEFDoseGrid') && all(cell2mat(Plan.CEFDoseGrid)))
      DoseFileName = ['dose_UserGrid_spt' , num2str(spt) , '_spt_X_',num2str(round(SptPos(1))),'_Y_',num2str(round(SptPos(2)))];
  else
      DoseFileName = ['dose_HighRes_spt' , num2str(spt) , '_spt_X_',num2str(round(SptPos(1))),'_Y_',num2str(round(SptPos(2)))];
  end
  DosePath = fullfile(outputPath,'CEF_beam','Outputs');
  DoseFullFileName = fullfile(DosePath,[DoseFileName,'.dcm']);

  %Add the CEF to the CT scan
  if WholeField
    %Interpolate the whole CT scan in one go
    [minField , maxField] = getMaxBEVsize(Plan.Beams);
  else
    %Interpolate only around the chosen beamlet
    FieldSize = 30; % mm This radius must be larger than the spot radius at the distal surface of PTV. Otherwise the lattice structure will be visible in the dose map
    minField = SptPos - [FieldSize , FieldSize];
    maxField = SptPos + [FieldSize , FieldSize];
  end

  %Insert the CEM into the high resolution CT scan
  % Define min and max field so that we do not expand the CT scan in the Xg and Yg direction to fit the CEM
  HUair =  getMaterialSPR('Schneider_Air' , Plan.ScannerDirectory) + 1; %Hounsfield unit associated to air in the material file
  hrCTName = 'highResCTbev';

  %Generate the high resolution CT scan
  %The IEC gantry is aligned with the Y axis of the CT scan
  fprintf('Interpolating high resolution CT \n')
  [handlesHR , BeamHR , iCTgntY , iCTgntX , iCTgntZ] = createHighResCT(handles , CTName , hrCTName , PTV, Plan.Beams , Plan.Spike.intrpCTpxlSize , HUair , minField , maxField , Plan.CTinfo);

  PlanHR.Beams = BeamHR;
  PlanHR = copyFields(PlanHR , Plan); %Overwrite the YAML parameter in the reloaded |Plan|
  PlanHR.CTname =  hrCTName;
  PlanHR  = updatePlanCTparam(handlesHR , PlanHR );

  %Add the CEM into the high resolution CT
  fprintf('Adding CEM to high resolution CT\n')
  [PlanHR , handlesHR ] = setCEMinCT(handlesHR , PlanHR , hrCTName , minField , maxField);

  %Insert the range shifter in the high resolution CT
  fprintf('Adding Range shifter to high resolution CT \n')
  [PlanHR , handlesHR] = setRangeShifterinCT(handlesHR , PlanHR , hrCTName, minField , maxField);

  %Get the high res CT
  CTintrp = Get_reggui_data(handlesHR, hrCTName ,'images');
  ImagePositionPatientINTER = handlesHR.origin;
  spacingINTER = handlesHR.spacing;
  isocenterINTER =  handlesHR.isocenter;

  infoCTinterp = PlanHR.CTinfo;
  infoCTinterp.ImagePositionPatient = ImagePositionPatientINTER;
  infoCTinterp.Spacing = handlesHR.spacing;

  %Save the interpolated CT on disk
  %The IEC gantry axis of the beamlet is aligned with the Y axis of this high resolution Ct scan
  if (~exist(fullfile(outputPath,'CEF_beam'),'dir'))
    %The folder to save the CT does not exist. Create it
    mkdir (fullfile(outputPath,'CEF_beam'))
  end
  if (~exist(DoseFullFileName))
    %The dose file does not already exist. Save the high resolution CT
    save_Image(CTintrp , infoCTinterp , fullfile(outputPath,'CEF_beam','CT_with_CEF'),'dcm');
  else
    %The output folder exists. Skip saving the CT
  end

  %Clean memory
  CTintrp = [];

  %Save the monolayer plan to disk
  PlanMono = saveMonoLayerPlan(Plan , spt ,  isocenterINTER , ImagePositionPatientINTER , spacingINTER , WholeField  , outputPath);


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
  if WholeField
      %Do nothing
  elseif (isfield(Plan, 'CEFDoseGrid') && all(cell2mat(Plan.CEFDoseGrid)))
      %Define the MCsquare parameter to do the dose scoring on a different pixel size than the high resolution CT scan
      fprintf('Dose scoring grid \n')
      Plan2.Independent_scoring_grid = Plan.Independent_scoring_grid; % Enable a different scoring grid than the base CT for the dose calculation
      Plan2.Scoring_voxel_spacing = Plan.Scoring_voxel_spacing; % In [mm]. Set dose calcuation scorinng grid to 1mm spacing and overwrite the CT grid. This would reduce computation time if the CT resolution is very high.
  end

  % Compute the dose map if no dose file alreadyt exists
  if (~exist(DoseFullFileName))
    %The output folder does not exist. Compute the high resolution dose
    handlesHR = ComputeFinalDose(Plan2, handlesHR , DoseFileName);
    [Dose , DoseInfo ] = Get_reggui_data(handlesHR,'dose_final_miropt');
                    %Get the image from handlesHR.myData. This is the idx-th element with |handlesHR.myData.name{idx} = 'dose_final_miropt'|
                    %This is an image at the resolution of the scoring grid defined in |CEFDoseGrid|
                    %It has not been resample at the resolution of the HR CT scan.
  else
    %The output folder exists. Skip dose computation. Simply relaod the file
    fprintf('File %s .dcm already exists. Skipping dose computation \n' , DoseFileName)
    [Dose , DoseInfo ] = load_Image(DosePath,DoseFileName,'dcm',0);
  end

  %Cleanup temporary file if required
  if ~Plan.SaveHighResDoseMap
    %If the user does not want to keep the high res dose map in the refernece frame of the beamlet
    %then delete it
    fileName = fullfile(DosePath,[DoseFileName,'.dcm']);
    fprintf('Deleting %s \n',fileName)
    delete (fileName)
  end

  %Rotate the dose map to be aligned with the IEC gantry again
  Dose = flipdim(Dose,2);
  Dose = permute(Dose,[1,3,2]);

  DoseHR.size = size(Dose);
  origZ = -DoseInfo.ImagePositionPatient(2) - DoseInfo.Spacing(2) .* (DoseHR.size(3) - 1);
  DoseHR.origin = [DoseInfo.ImagePositionPatient(1) , DoseInfo.ImagePositionPatient(3) , origZ ]'; %Second index is minus (because flipped) Zg
  DoseHR.isocenter = [handlesHR.isocenter(1) , handlesHR.isocenter(3) , handlesHR.isocenter(2)]; %permuted the CT scan
  DoseHR.spacing = [ DoseInfo.Spacing(1) , DoseInfo.Spacing(3) , DoseInfo.Spacing(2)]'; %permutation has no effect here: all elements are identical

  %Compute the coordinate of all pixels of the original CT
  %The DoseOrig dose map is aligned with the axes of the original DICOM CT
  [oCTdcmX, oCTdcmY , oCTdcmZ , X , Y , Z] = getCTaxes(handles.origin , handles.spacing , handles.size , Plan.Beams.isocenter); %Create coordinate system of original CT, centered on original isocenter

  A = [X(:), Y(:), Z(:)];
  A = [A , ones(numel(X),1)]; %Coordinate of all voxels of original CT expressed in DICOM CS
  Beam2 = Plan.Beams;
  Beam2.isocenter = [0,0,0];
  oCTgntr = DICOM2gantry(A' , Beam2); %Coordinate of all voxels of the original CT scan, expressed in IEC gantry. Rotate arouind isocenter at [0,0,0], using gantry angle =0, PPS=0

  %The pixels of the high resolution dose map are aligned on the IEC gantry coordinates
  %Let's recompute the coordinates of original CT scan from DICOM to IEC gantry
  [iCTgntX, iCTgntY , iCTgntZ ] = getCTaxes(DoseHR.origin , DoseHR.spacing , DoseHR.size , DoseHR.isocenter);

  %Interpolate the dose map into the original CT scan
  DoseOrig = interpolateCT(Dose, iCTgntX, iCTgntY , iCTgntZ , oCTgntr(1,:) , oCTgntr(2,:) , oCTgntr(3,:) , 0 , handles.size , handles.spacing , 'linear'); %The dose in the orignal grid

end

%-----------------------------------
% Compute the maximum BEV field size (IEC gantry)
% in which the dose should be computed
% If there is an aperture, use the max aperture size
% If there is no aperture, then use the CEF size
%------------------------------------
function [minField , maxField] = getMaxBEVsize(Beam)

  if (isfield(Beam,'BlockData'))
      %There is an aperture block. The field size is defined by the aperture block
      [minField , maxField ] = findApertureBlockSize(Beam.BlockData);

  else
      %There is no aperture block. The field size is defined by CEF size
      [~ , ~ , ~ , BlockBorder] = findApertureBlockSize({}); %Retrieve the border size from this function

      NbPxlCEF = size(Beam.RangeModulator.CEMmask); %[Nx,Ny,Nz] number of pixels along X and Y IEC gantry
      SizeCEF = Beam.RangeModulator.Modulator3DPixelSpacing(1:2) .* (NbPxlCEF-1); %[Sx,Sy] (mm) dimension of the CEF block along X and Y IEC gantry
      minCEF = Beam.RangeModulator.ModulatorOrigin; %One point of the main diagonal of the CEF
      maxCEF = Beam.RangeModulator.ModulatorOrigin + SizeCEF; %Other point of the main diagonal of the CEF

      minField = minCEF(1:2) - BlockBorder;
      maxField = maxCEF(1:2) + BlockBorder;
  end
end

%========================
% Convert coordinates of points from the DICOM CS to the IEC gantry CS
%========================
function Adcm = DICOM2gantry(A , Beam)
    M = matDICOM2IECgantry(Beam.GantryAngle,Beam.PatientSupportAngle);
    Adcm = M * A ;

    %Coordinate of voxels in IEC gantry
    Adcm(1,:) = Adcm(1,:)  + Beam.isocenter(1);
    Adcm(2,:) = Adcm(2,:)  + Beam.isocenter(2);
    Adcm(3,:) = Adcm(3,:)  + Beam.isocenter(3);

end


%===========================
% Save the mono layer plan to use with the high resolution CT
%
% The monolayer plan is linked to the high resolution Ct which assumes gantry angle 0°
% The modifications made to |Plan| are local to this function and will not affect the calling function
%=========================
function PlanMono = saveMonoLayerPlan(Plan , spt , isocenterINTER , ImagePositionPatientINTER , spacingINTER , WholeField , outputPath)

  fprintf('Saving DICOM plan with monolayer \n')
  Plan.Beams = Plan.Beams; %Only keep the chosen beam
  %The interpolated CT was rotated so that the CEM axis is paralell to the Ydicom axis
  %We will use a virtualm gantry at 0° and no couch rotation to shout through the CEM in the interpolated CT
  %TODO This asusmes that the BDL does not depend on the gnatry angle !!!!!!!!!!!!!
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
