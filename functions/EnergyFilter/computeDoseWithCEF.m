%% computeDoseWithCEF
% Compute the dose distribution using MCsquare for a plan with a single beam and a single energy layer
% and a conformal energy modulator.
% The planning CT scan is re-interpolated on a pixel grid with the pixel size of the CEM in order to properly model the CEM and with the axes oriented along the IEC gantry CS.
% In order to save memory and avoid out-of-memory errors, the interpolation of the high resolution CT is done around each PBS bemlet sequentially.
% The dose delivered by each beamlet is tallied on a grid with resolution |Plan.CEFDoseGrid| and with ortientation ofthe IEC gantry CS.
% The dose of each beamlet is then added to the global dose map (a grid with resolution |Plan.CEFDoseGrid|) and with axes oriented along the IEC gantry CS.
% Note that for this intermediate dose computation, the 3 components of |Plan.CEFDoseGrid| represent pixels sizer in IEC gantry CS.
%
% The global dose map is then rotated (re-interpolated) to have pixels aligned with the DICOM axes of the original CT scan and with pixel resolution |Plan.CEFDoseGrid|.
% Note that following this re-interpolation the 3 components of |Plan.CEFDoseGrid| represent pixels sizer along the axes of the DICOM CS of the original Ct scan.
% Note therefore that |Plan.CEFDoseGrid| define the spatial resolution of two different dose map, aligned along different coordinate systems.
%
% The intermediate PBS beamlet dose map (aligned with IEC gantry) is optionally saved depending on the flag  |Plan.SaveHighResDoseMap|
% The intermediate  high resolution CT scan (aligned with IEC gantry) is optionally saved depending on the flag |Plan.SaveHighResCT|.
% If is possible to save intermediate PBS beamlet resolution dose map (aligned with axes of the original CT scan) with the flag |Plan.SaveDoseBeamlets = 'dcm'|. NOgte however that this increases the computation time.
% Note that saving the intermediate dose maps and CT scan increase the execution time of the function.
%
%% Syntax
% |[Plan , MinDose , MaxDose , DoseOrig] = computeDoseWithCEF(Plan, outputPath, handles, CTName , FLAGdosePerSpot)|
%
%
%% Description
% |[Plan , MinDose , MaxDose , DoseOrig] = computeDoseWithCEF(Plan, outputPath, handles, CTName , FLAGdosePerSpot)| Compute dose maps and save them to disk
%
%
%% Input arguments
%
% |Plan| -_STRUCTURE_- Structure defining a mono energy layer PBS plan
%
% |outputPath| -_STRING_- Path to the output folder where the DICOM file of the dose map and the mono-layer plan will be saved
%
% |handles| -_STRUCTURE_- REggui data handle. The CT scan is stored in |handles| in the image with name |Plan.CTname|.
%
% |CTname| -_STRING_- Name of the CT image in handles.images
%
% |FLAGdosePerSpot| -_BOOLEAN_- TRUE : compute high res dose per beamlet. FALSE : compute high res dose in whole volume
%
%
%% Output arguments
%
% |Plan| -_SCALAR MATRIX_- Update the value of |Plan.Scenario4D(1).RandomScenario(Plan.rr_nominal).RangeScenario(Plan.rs_nominal).P| with the dose influence matrix: |Pij(vox,spot)| The dose contribution to voxel |vox| of the spot number |spot|
%
% |MinDose| -_SCALAR_- Dose (Gy) of the 0.1% lower percentile
%
% |MaxDose| -_SCALAR_- Dose (Gy) of the 0.1% higher percentile
%
% |DoseOrig| -_SCALAR MATRIX_- Dose map (Gy) with the spatial resolution of |Plan.CEFDoseGrid| and with pixels aligned with the DICOM CS of the original CT scan
%
%% Contributors
% Authors : R. Labarbe, Lucian Hotoiu (open.reggui@gmail.com)

function [Plan , MinDose , MaxDose , DoseOrig] = computeDoseWithCEF(Plan, outputPath, handles, CTName , FLAGdosePerSpot)

    global g_HUair;
    global g_HUbrass;
    global g_HUcem; %Define HU as a global variable that will be visible inside the function getHighResDose
    global g_HUrangeshifter;
    %getMaterialPropCT reads the disk. Only make the reading once and not at every iteration of the loop for each beamlet in order to make computation faster
    g_HUair =  getMaterialPropCT('Schneider_Air' , Plan.ScannerDirectory ) + 1; %Hounsfield unit associated to air in the material file
    g_HUcem =  getMaterialPropCT(Plan.Spike.MaterialID , Plan.ScannerDirectory) +1 ; %Hounsfield unit associated to CEM in the material file
    g_HUbrass = getMaterialPropCT('Brass' , Plan.ScannerDirectory) + 1 ; %Hounsfield unit associated to brass in the material file
    g_HUrangeshifter =  getMaterialPropCT(Plan.Beams(1).RSinfo.RangeShifterMaterial , Plan.ScannerDirectory) + 1 ; %HU and relative stopping power of the range shifter
        %The g_HUair and g_HUcem and g_HUrangeshifter are now stored in a global varialbe that is accessible by the sub-function getHighResDose
        %We assume that all beams use the same material for the range shifter

    %Make some sanity check on the treatment plan
    if ~isfield(Plan, 'CEFDoseGrid')
      error('Plan.CEFDoseGrid missing : dose map resolution is undefined')
    end

    if nargout > 0
      %We need to compute and return the dose influence matrix
      PijFlag = true;
      FLAGdosePerSpot = true; %If we want the Pij matrix, we need to compute the dose per spot
    else
      %Do not waste time computing the dose influence matrix
      PijFlag = false ;
    end

    if ~isfield(Plan , 'SaveDoseBeamlets')
      %If the field is not defined, this means we do not want to save the dose map of the beamlets
      % in the reference frame of the CT scan
      Plan.SaveDoseBeamlets = false;
    end

    %Save the monolayer plan to the output folder
    %This will be used to reference to dose map to it
    % if (~exist(fullfile(outputPath,'CEF_beam'),'dir'))
    %   %The folder to save the CT does not exist. Create it
    %   mkdir (fullfile(outputPath,'CEF_beam'))
    % end
    Plan.FileName = 'Plan';
    regguiPath = fileparts(which('reggui'));
    dictionary = fullfile(regguiPath ,'plugins','openMIROpt','functions','io','dicom-dict.txt');
    createDICOMPlan(Plan,Plan.CTinfo,outputPath,dictionary)


    %REcord the parameters of the original CT scan
    hCT = createREGGUIhandles(handles.size , handles.origin , handles.spacing); %Make first a copy of the size of the original CT scan. We want the dose map to match this size

    if ~FLAGdosePerSpot
        % Compute the dose in the whole volume in one go
        error('Not implemented')
        % [Plan , ~ , ~ , ~ , minField , maxField] = getHRdoseGridInfo(Plan , handles, Zdistal , b);
        % [DoseIECg , DoseFileName , handlesDoseIECg] = getHighResDose(Plan, outputPath , handles , CTName , 0 , minField , maxField , Zdistal , true);
        % [iDoseGntX , iDoseGntY , iDoseGntZ ] = getCTaxes(handlesDoseIECg.origin , handlesDoseIECg.spacing , handlesDoseIECg.size , [0,0,0]); %Cooridnates of dose map aligned with IEC gantry
        %
        % hD = createREGGUIhandles(ceil(hCT.size .* hCT.spacing ./ handlesDoseIECg.spacing) , hCT.origin , handlesDoseIECg.spacing);
        %             %The final dose map will have the spatial resolution of |Plan.CEFDoseGrid|
        %             %Compute the number of pixels at dose map resolution to fill the volume of original CT
        % [DoseOrig , handlesDose]  = doseIECg2DICOMcs(DoseIECg  , handlesDoseIECg , hD , Plan , outputPath); %Dose map aligned to orignal Ct with resolution of |Plan.CEFDoseGrid|
     else
        % Compute the dose in one beamlet at a time
        fprintf('Computing high resolution dose map one spot at a time \n')
        [DoseOrig, DoseFileName , handlesDose, Pij] = getFullHighResDosemap(Plan , handles , outputPath, CTName , PijFlag , hCT);
    end

    %Find the intensity of the 0 99.9% percentile dose
    [MinDose , MaxDose] = get_image_scale({DoseOrig},0.1);

    % Save to plan the Pij matrix with beamlets through CEM
    Plan.Scenario4D(1).RandomScenario(Plan.rr_nominal).RangeScenario(Plan.rs_nominal).P = Pij;


    %Save dose map in original grid
    planFullPath = fullfile(outputPath,Plan.FileName);
    handlesDose = save2Disk(handlesDose, DoseOrig , handlesDose.size , Plan.CTinfo , 'Dose_withCEF' , outputPath , planFullPath , 'RTDOSE');
    %handlesDose = save2Disk(handlesDose, DoseOrig , handlesDose.size , Plan.CTinfo , 'Dose_withCEF' , fullfile(outputPath,'CEF_beam') , planFullPath , 'RTDOSE');

    %Copy the output file to the 'Output' folder
    %movefile (fullfile(outputPath,'CEF_beam', 'Dose_withCEF.dcm' ) , fullfile(outputPath,'Dose_withCEF.dcm'));

  end % This is the end of the computeDoseWithCEF


%-----------------------------------------------
%Define the parameters of the high resolution dose map
%
% OUTPUT
% |iDoseGntX| -_SCALAR VECTOR_- |iDoseGntX(i)| X coordinate (mm) in IEC gantry CS of the i-th voxel of the DOSE MAP
% |iDoseGntY| -_SCALAR VECTOR_- |iDoseGntY(i)| Y coordinate (mm) in IEC gantry CS of the i-th voxel of the DOSE MAP
% |iDoseGntZ| -_SCALAR VECTOR_- |iDoseGntZ(i)| Z coordinate (mm) in IEC gantry CS of the i-th voxel of the DOSE MAP
% |minField| -_SCALAR VECTOR_- [x,y,z] coordinate (mm) in the IEC gantry CS of the first corner of the field where the dose is to be computed. Must be aligned with a pixel of the CEM
% |maxField| -_SCALAR VECTOR_- [x,y,z] coordinate (mm) in the IEC gantry CS of the last corner of the field where the dose is to be computed.
% |hD| -_STRUCTURE_- Dimensions of the high resolution dose map
%     * |hD.size| -_SCALAR VECTOR_- [Nx,Ny,Nz] Number of voxels of the |DoseDCMcs|
%     * |hD.origin| -_SCALAR VECTOR_- [x,y,z] Coordinate (mm) of the first pixel of the |DoseDCMcs| in the original CT CS
%     * |hD.spacing|-_SCALAR VECTOR_- [dx,dy,dz] Physical size (mm) of the voxels of |DoseDCMcs|
% |Plan| -_STRUCTURE_- Structure defining a mono energy layer PBS plan
%-----------------------------------------------
function [Plan , iDoseGntX , iDoseGntY , iDoseGntZ, minField , maxField, hD , Zdistal] = getHRdoseGridInfo(Plan , handles, b)

  % Set coarse grid for dose calculation.
  Plan.Independent_scoring_grid = 1; % Enable a different scoring grid than the base CT for the dose calculation
  Plan.resampleScoringGrid = false; %Do not resample the dose map in MC_compute. Leave it on the coarse scoring grid
  Plan.Scoring_voxel_spacing = cell2mat(Plan.CEFDoseGrid); % In [mm]. Set dose calcuation scorinng grid to 1mm spacing and overwrite the CT grid. This would reduce computation time if the CT resolution is very high.
  fprintf('Dose scoring grid: [%f , %f , %f ] mm \n',Plan.Scoring_voxel_spacing(1),Plan.Scoring_voxel_spacing(2),Plan.Scoring_voxel_spacing(3));

  %GEt the deepest Zg at which the dose should be computed
  Body = Get_reggui_data(handles , Plan.ExternalROI);
  Zdistal = getZdistal(Body , handles.spacing , handles.origin , Plan.Beams(b)); % |Zdistal| -_SCLAR_-  Z Coordinate (mm) in the IEC gantry CS of the deepest plane in which the dose is to be computed

  %Create a dose map at the resolution |Plan.Scoring_voxel_spacing|
  [minField , maxField] = getMaxBEVsize(Plan.Beams(b));
  [iDoseGntX , iDoseGntY , iDoseGntZ] =  getDoseMapCoordInIECg(Plan.Beams(b) , Zdistal , handles.spacing , handles.origin , minField , maxField ,  Plan.Scoring_voxel_spacing);

  hD = createREGGUIhandles([numel(iDoseGntX) , numel(iDoseGntY) , numel(iDoseGntZ)] , [iDoseGntX(1) , iDoseGntY(1) , iDoseGntZ(1)] , Plan.Scoring_voxel_spacing);

end


%==============================================
% Compute high resolution dose map in the IEC gantry CS, one beamlet at a time
% Do not intrepolate for each beamlet, just insert at the proper place in the map
% then only rotate at the end to place the dose map in the original CT orientation
%while keeping the dose map resolution
%
% The function return a dose map with a aspatial resolution of |Plan.Scoring_voxel_spacing|
% and aligned wit the coordinate axes of the original CT scan
%
% INPUT
% |Plan| -_STRUCTURE_- Structure defining a multi energy layer PBS plan
%
% |handles| -_STRUCTURE_- REggui data handle. The CT scan is stored in |handles| in the image with name |Plan.CTname|.
%
% |Zdistal| -_SCLAR_-  Z Coordinate (mm) in the IEC gantry CS of the deepest plane in which the dose is to be computed
%
% |outputPath| -_STRING_- Path to the folder where the DICOM file with the mono-layer plan will be saved
%
% |CTname| -_STRING_- Name of the CT image in handles.images
%
% |PijFlag| -_BOOL_- true = compute and return the dose influence matrix
%
% |hCT| -_STRUCTURE_- dimension of the original CT scan
%     * |hCT.size| -_SCALAR VECTOR_- [Nx,Ny,Nz] Number of voxels
%     * |hCT.origin| -_SCALAR VECTOR_- [x,y,z] Coordinate (mm) of the first pixel
%     * |hCT.spacing|-_SCALAR VECTOR_- [dx,dy,dz] Physical size (mm) of the voxels
%================================================
function [DoseOrigCT, DoseFileName , handlesDose , Pij] = getFullHighResDosemap(Plan , handles , outputPath , CTName , PijFlag, hCT)

  FieldSize = 60; % mm This radius must be larger than the spot radius at the distal surface of PTV. Otherwise the lattice structure will be visible in the dose map


  %Loop for each beam in the plan
  %The setup beams are removed when loading the plan. Only treatment beams are present here.
  DoseOrigCT = 0;

  if PijFlag
    %If we need to compute the dose influence matrix, prepare an empty sparse matrix
    [Pij , w] = preparePij(Plan , hCT);
    jIndex = 0;
  else
    %No one will use Pij. It can be empty
    Pij = [];
  end

  for b = 1:numel(Plan.Beams)

      %Path where to save the dose maps for the individual beams
      path2beamResults = getOutputDir(outputPath , b);

      %Make some sanity check on the treatment plan
      if (numel(Plan.Beams(b).Layers) ~= 1)
        error('computeDoseWithCEF requires plans with single energy layer')
      end

      NbBeamlets = numel(Plan.Beams(b).Layers(1).SpotWeights); %Number of PBS spots in plan with energy monolayer

      % Set coarse grid for dose calculation.
      [Plan , iDoseGntX , iDoseGntY , iDoseGntZ , ~ , ~ , hDoseIECg , Zdistal] = getHRdoseGridInfo(Plan , handles, b);
      DoseIECg = zeros(numel(iDoseGntX) , numel(iDoseGntY) , numel(iDoseGntZ) ); %Create an empty dose map in which the dose from each bemalet will be saved

      %Check whether we need to save beamlet dose map in the original CT CS
      switch Plan.SaveDoseBeamlets
          case {'dcm' , 'sparse'}
            saveBeamlets = true;
          otherwise
            %We do not want to save the beamlets. Do not waste time reinterpoalting dose map in oriiginal CT scan CS
            saveBeamlets = false;
      end

      for spt = 1:NbBeamlets

          SptPos = Plan.Beams(b).Layers(1).SpotPositions(spt,:);
          jIndex = jIndex +1;

          %Define a field size that match the pixel position of the CEM
          %The CEM can be copied pixel wise from the mask into the high resolution CT
          %This will avoid the alising problems when accumulating the dose
          minField = alignFieldWithCEM(SptPos - FieldSize , Plan.Beams(b).RangeModulator);
          maxField = SptPos + FieldSize; %The maxfield will be automatically aligned with the CEM pixel when calling minField:step:maxfield
          Zdistal = getClosePixel(iDoseGntZ , Zdistal);  %Zdistal is aligned on the global dose map grid

          fprintf('Computing dose for beam %d , PBS spots %d of %d at [%f mm, %f mm] \n',b , spt,NbBeamlets,SptPos(1),SptPos(2))
          [doseSpt, DoseFileName , hDoseSpt] = getHighResDose(Plan, b , path2beamResults, handles, CTName, spt, minField , maxField , Zdistal , false , iDoseGntX , iDoseGntY , iDoseGntZ); %Add the dose of each spot to the global high resolution dose map

          %Align the partial dose map on the pixel of the full dose map
          % The center of the voxels of the  partial dose map are already aligned with the voxels of the global dose map.
          %We only need to find the index of the voxel where to add the partial dose map
          [~ , IndexStart(1)] = getClosePixel(iDoseGntX , hDoseSpt.origin(1));
          [~ , IndexStart(2)] = getClosePixel(iDoseGntY , hDoseSpt.origin(2));
          [~ , IndexStart(3)] = getClosePixel(iDoseGntZ , hDoseSpt.origin(3));

          DoseIECg(IndexStart(1):IndexStart(1)+hDoseSpt.size(1)-1 , IndexStart(2):IndexStart(2)+hDoseSpt.size(2)-1 , IndexStart(3):IndexStart(3)+hDoseSpt.size(3)-1) = ...
                             DoseIECg(IndexStart(1):IndexStart(1)+hDoseSpt.size(1)-1 , IndexStart(2):IndexStart(2)+hDoseSpt.size(2)-1 , IndexStart(3):IndexStart(3)+hDoseSpt.size(3)-1) + doseSpt;

          if PijFlag | saveBeamlets
            %Reinterpolate the dose map (in IECg) into the CS of the original CT scan
            % The final dose map will have the resolution of the original CT scan
            [DoseSptDCMcs, hDoseDCM] = doseIECg2DICOMcs(doseSpt , hDoseSpt , hCT , Plan , b , path2beamResults);
          end

          if PijFlag
            %Add this beamlet to the dose influence matrix
            Pij = addBeamlet2Pij(Pij , spt , jIndex , DoseSptDCMcs , Plan , w);
          end

          %If required, save the dose map of each beamlet in the reference frame of the original CT scan
          DosePath = fullfile(path2beamResults,'CEF_beam','Outputs','SpotDoseInCT');
          switch Plan.SaveDoseBeamlets
            case 'dcm'
              if (~exist(DosePath,'dir'))
                %The folder to save the CT does not exist. Create it
                mkdir (DosePath)
              end
              planFullPath = fullfile(path2beamResults,'CEF_beam',Plan.FileName);
              save2Disk(hDoseDCM, DoseSptDCMcs , hDoseDCM.size , Plan.CTinfo , DoseFileName , DosePath , planFullPath, 'RTDOSE'); %Save the beamlet

            case 'sparse'
                if (~exist(DosePath ,'dir'))
                  %The folder to save the CT does not exist. Create it
                  mkdir (DosePath )
                end
                save(fullfile(DosePath , [DoseFileName '.mat']) , 'DoseSptDCMcs')

             otherwise
                %Do not save anything
          end %switch Plan.SaveDoseBeamlets

      end %end for spt

      %Rotate the full dose map to make axes paralell to original CT scan
      hD = createREGGUIhandles(ceil(hCT.size .* hCT.spacing ./ hDoseIECg.spacing) , hCT.origin , hDoseIECg.spacing);
                %The final dose map will have the spatial resolution of |Plan.CEFDoseGrid|
                %Compute the number of pixels at dose map resolution to fill the volume of original CT

      [DoseOrigCTloc, handlesDose] = doseIECg2DICOMcs(DoseIECg  , hDoseIECg , hD , Plan , b , path2beamResults);
      DoseOrigCT = DoseOrigCT + DoseOrigCTloc;
  end %for b (beams)



end

%------------------------------------------
%REsample the dose map with pixels aligned on the IEC gantry CS
% to a dose map aligned on the DICOM CS
%
% INPUT
% |DoseIECg| -_SCALAR MATRIX_- DoseIECg(x,y,z) Dose (Gy) at the pixel [x,y,z]. Matrix axes are aligned with IEC gantry CS
% |handlesDoseIECg| -_STRUCTURE_- REGGUI handle defining the spacing, origin and size of the |DoseIECg| matrix
% |hD| -_STRUCTURE_- dimensions of the output dose |DoseDCMcs|
%     * |hD.size| -_SCALAR VECTOR_- [Nx,Ny,Nz] Number of voxels of the |DoseDCMcs|
%     * |hD.origin| -_SCALAR VECTOR_- [x,y,z] Coordinate (mm) of the first pixel of the |DoseDCMcs| in the original CT CS
%     * |hD.spacing|-_SCALAR VECTOR_- [dx,dy,dz] Physical size (mm) of the voxels of |DoseDCMcs|
% |Plan| -_STRUCTURE_- Structure defining a mono energy layer PBS plan
% |outputPath| -_STRING_- Path to the output folder where the DICOM file of the dose map and the mono-layer plan will be saved
%
% OUTPUT
% |DoseDCMcs| -_SCALAR MATRIX_- DoseDCMcs(x,y,z) Dose (Gy) at the pixel [x,y,z]. Matrix axes are aligned with the DICOM CS of the CT scan
%------------------------------------------
function [DoseDCMcs, handlesDose] = doseIECg2DICOMcs(DoseIECg , handlesDoseIECg , hD , Plan , b , outputPath)

  %Create a fake handle for the dose map in the original CT CS
  handlesDose = struct;
  handlesDose.path = outputPath;
  handlesDose = Initialize_reggui_handles(handlesDose); % Initialize the handles structure
  handlesDose.spacing = hD.spacing;
  handlesDose.origin = hD.origin;
  handlesDose.size = hD.size ;
  handlesDose.spatialpropsettled = true;

  %Compute the coordinate of all pixels of the full dose map
  %The DoseOrig dose map is aligned with the axes of the original DICOM CT
  [oCTdcmX, oCTdcmY , oCTdcmZ , X , Y , Z] = getCTaxes(handlesDose.origin , handlesDose.spacing , handlesDose.size , Plan.Beams(b).isocenter); %Create coordinate system of original CT, centered on original isocenter
  A = [X(:), Y(:), Z(:)];
  A = [A , ones(numel(X),1)]; %Coordinate of all voxels of full dose map in DICOM CS
  Beam2 = Plan.Beams(b);
  Beam2.isocenter = [0,0,0]; %The oCTdcm CS has its orgin at isocenter. We need to reset the isocenter to the origin of the CS.
  oCTgntr = DICOM2gantry(A' , Beam2); %Coordinate of all voxels of the original CT scan, expressed in IEC gantry. Rotate arouind isocenter at [0,0,0], using gantry angle =0, PPS=0

  %Cooridnates of dose map aligned with IEC gantry
  % The |DoseIECg| has its origin at  the isocenter. So isocenter = [0,0,0]
  [iDoseGntX , iDoseGntY , iDoseGntZ ] = getCTaxes(handlesDoseIECg.origin , handlesDoseIECg.spacing , handlesDoseIECg.size , [0,0,0]);

  %Interpolate the dose map into the original CT scan
  DoseDCMcs = interpolateCT(DoseIECg, iDoseGntX, iDoseGntY , iDoseGntZ , oCTgntr(1,:) , oCTgntr(2,:) , oCTgntr(3,:) , 0 , handlesDose.size  , 'linear'); %The dose in the orignal grid

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
  % |iDoseGntX| -_SCALAR VECTOR_- |iDoseGntX(i)| X coordinate (mm) in IEC gantry CS of the i-th voxel of the DOSE MAP
  % |iDoseGntY| -_SCALAR VECTOR_- |iDoseGntY(i)| Y coordinate (mm) in IEC gantry CS of the i-th voxel of the DOSE MAP
  % |iDoseGntZ| -_SCALAR VECTOR_- |iDoseGntZ(i)| Z coordinate (mm) in IEC gantry CS of the i-th voxel of the DOSE MAP
  %
  % OUTPUT
  %
  % |Dose| -_SCALAR MATRIX_- Dose(x,y,z) Dose (Gy) at the pixel [x,y,z]. Matrix axes are aligned with IEC gantry CS
  %
  % |DoseFileName| -_STRING_- File name of the dose map computed by MCsquare at the resolution |scoring_voxel_spacing| and in the IEC gantry CS
  %
  % |DoseHR| -_STRUCTURE_- Handles with the image properties of the dose map at the resolution |scoring_voxel_spacing| and in the IEC gantry CS
  %===========================================================

  function [Dose, DoseFileName , DoseHR] = getHighResDose(Plan, b , outputPath , handles , CTName , spt , minField , maxField , Zdistal , WholeField , iDoseGntX , iDoseGntY , iDoseGntZ)

    global g_HUair;
    global g_HUcem; %Define HU as a global variable that will be visible inside the function getHighResDose

    %Define the name of the output dose file
    if WholeField
        DoseFileName = 'dose_HighRes';
    else
        SptPos = Plan.Beams(b).Layers(1).SpotPositions(spt,:);
        DoseFileName = ['dose_UserGrid_spt' , num2str(spt) , '_spt_X_',num2str(round(SptPos(1))),'_Y_',num2str(round(SptPos(2)))];
    end
    DosePath = fullfile(outputPath,'CEF_beam','Outputs');
    DoseFullFileName = fullfile(DosePath,[DoseFileName,'.dcm']);

    %Insert the CEM into the high resolution CT scan
    % Define min and max field so that we do not expand the CT scan in the Xg and Yg direction to fit the CEM
    hrCTName = 'highResCTbev';

    %Check that the field position is aligned with the pixel resolution and origin of CEM
    a = (minField - Plan.Beams(b).RangeModulator.ModulatorOrigin(1:2)) ./ Plan.Beams(b).RangeModulator.Modulator3DPixelSpacing(1:2);
    b1 = round( (minField - Plan.Beams(b).RangeModulator.ModulatorOrigin(1:2)) ./ Plan.Beams(b).RangeModulator.Modulator3DPixelSpacing(1:2) );
    if (max(abs(a - b1)) > 1e-4)
      a
      b1
      a-b1
      error('Small field not aligned with CEM grid')
    end

    %Generate the high resolution CT scan
    %The proton beam is aligned with the **Y axis** of the CT scan : Xg = Xct; Yct = -Zg ; Zct = Yg
    % The X-Z spatial resolution of the Ct is the same as CEM to avoid aliasing problems when inserting CEM in high res CT
    % The Y resolution is 0.5mm. The CEM height is a multiple of 1mm and the rnage shifter are also multiple of 1mm
    % The dose map and original CT scan resolution are coarser tan 0.5mm.
    % Therefore with Z resolution 0.5mm, we can describes the fine structures of the CEM and range shifter.
    CTresolution = Plan.Beams(b).RangeModulator.Modulator3DPixelSpacing;
    CTresolution(3) = 0.5; %mm

    fprintf('Interpolating CT with pixels [%f , %f , %f ] mm \n', CTresolution(1) , CTresolution(2) , CTresolution(3))

    [handlesHR , BeamHR ] = createHighResCT(handles , CTName , hrCTName , Plan.Beams(b) , CTresolution , g_HUair , minField , maxField , Zdistal , Plan.CTinfo);

    PlanHR = Plan; %Copy all plan info
    BeamHR = copyFields(BeamHR , Plan.Beams(b)); %Copy all other beam info from original plan
    PlanHR.Beams = BeamHR; %Make sure there is only one beam
    PlanHR.CTname =  hrCTName;
    PlanHR  = updatePlanCTparam(handlesHR , PlanHR );

    % Add the aperture in the CT scan
    fprintf('Adding aperture to high resolution CT\n')
    handlesHR = setApertureinHRCT(handlesHR , PlanHR , hrCTName); %Add an apertrue block in the CT scan


    %Add the CEM into the high resolution CT
    fprintf('Adding CEM to high resolution CT\n')
    handlesHR = setCEMinhrCT(handlesHR , PlanHR , hrCTName , g_HUcem , g_HUair);

    %Add range shifter in the high resolution CT
    %This must be done in the high resolution CT to avoid the RS thickness to be aliased by the Z pixel resolution of the CT
    %As the CT is aligned with the IEC gantry CS, adding the RS can be done quickly by using the Matlab indexing.
    fprintf('Adding Range shifter to high resolution CT \n')
    handlesHR = setRangeShifterinHRCT(handlesHR , PlanHR , hrCTName);
    PlanHR.Beams.NumberOfRangeShifters = 0;  %Remove the range shifter from the MCsquare beam model. The range shifter is now inserted in the CT scan

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
    PlanMono = saveMonoLayerPlan(Plan , b , spt ,  handlesHR.isocenter , handlesHR.origin , handlesHR.spacing , WholeField  , outputPath);

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
    Plan2.MCsqExecPath = Plan.MCsqExecPath;
    Plan2.CTname = hrCTName;
    Plan2.output_path = fullfile(outputPath, 'CEF_beam');
    Plan2.ScannerDirectory = Plan.ScannerDirectory;
    Plan2.BDL = Plan.BDL;
    Plan2.protonsFullDose = Plan.protonsHighResDose;
    Plan2.CTinfo = PlanMono.CTinfo;

    %Compute the dose using MCsquare on the high resolution CT scan
    %Define the MCsquare parameter to do the dose scoring on a different pixel size than the high resolution CT scan
    Plan2.Independent_scoring_grid = Plan.Independent_scoring_grid; % Enable a different scoring grid than the base CT for the dose calculation
    Plan2.resampleScoringGrid = Plan.resampleScoringGrid; %Do not resample the dose map in MC_compute. Leave it on the coarse scoring grid. This save time in MC_compute
    Plan2.Scoring_voxel_spacing = [Plan.Scoring_voxel_spacing(1) , Plan.Scoring_voxel_spacing(3) , Plan.Scoring_voxel_spacing(2)];
                % In [mm]. Set dose calculation scorinng grid to 1mm spacing and overwrite the CT grid. This would reduce computation time if the CT resolution is very high.
                                %Plan.Scoring_voxel_spacing dimension are along the high res CT scan CS. The Zg is the Y axis of the Ct CS

    %Define the origin of the dose map scoring grid so that it is aligned with the voxels of the global dose map
    %This will avoid aliasing problem when copyting the spot dose map into the global dose map
    DoseMapOrig(1) = getClosePixel(iDoseGntX , handlesHR.origin(1));
    DoseMapOrig(2) = getClosePixel(-iDoseGntZ , handlesHR.origin(2));
    DoseMapOrig(3) = getClosePixel(iDoseGntY , handlesHR.origin(3));
    Plan2.Scoring_origin = [DoseMapOrig(1) , DoseMapOrig(2) , DoseMapOrig(3)]; %The first pixel of the spot dose map is algined with the pixels of the global dose map

    %Make the dose scoring grid as large as possible to fit into the high resolution CT scan
    Plan2.Scoring_grid_size =  round(((handlesHR.origin'  + handlesHR.spacing' .* (handlesHR.size - 1) ) - DoseMapOrig ) ./ Plan2.Scoring_voxel_spacing);

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
    DoseHR.spatialpropsettled = true;

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
function PlanMono = saveMonoLayerPlan(Plan , b , spt , isocenterINTER , ImagePositionPatientINTER , spacingINTER , WholeField , outputPath)

  fprintf('Saving DICOM plan with monolayer \n')

  %The interpolated CT was rotated so that the CEM axis is paralell to the Ydicom axis
  %We will use a virtualm gantry at 0° and no couch rotation to shout through the CEM in the interpolated CT
  %TODO This asusmes that the BDL does not depend on the gantry angle !!!!!!!!!!!!!
  Plan.Beams = Plan.Beams(b); %Choose the correct beam
  Plan.Beams.GantryAngle = 0;
  Plan.Beams.PatientSupportAngle = 0;
  Plan.Beams.isocenter = isocenterINTER; %Redefine the position of isocentre in the resampled CT

  %Remove the range shifter definition. It is now included inside the CT scan
  Plan.Beams.RSinfo.RangeShifterSetting ='OUT';
  Plan.Beams = rmfield(Plan.Beams,'RSinfo');
  Plan.Beams.NumberOfRangeShifters = 0;

  PlanMono = Plan;
  PlanMono.fractions = Plan.fractions(b);
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
%
% INPUT
% |Zdistal| -_SCLAR_-  Z Coordinate (mm) in the IEC gantry CS of the deepest plane in which the dose is to be computed
%================================================
function [iDoseGntX , iDoseGntY , iDoseGntZ ] =  getDoseMapCoordInIECg(Beam , Zdistal , Spacing , ImagePositionPatient , minField , maxField ,  PixelSizeIECg)

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

  %Define the coordinate of the voxels in the IEC gantry CS for the interpolated Ct scan
  %Make the dose map a bit larger to be sire to fit all bemalets
  iDoseGntX = (minField(1) - 15 .* PixelSizeIECg(1) )      : PixelSizeIECg(1) : (maxField(1) + 15 .* PixelSizeIECg(1));
  iDoseGntY = (minField(2) - 15 .* PixelSizeIECg(2) )      : PixelSizeIECg(2) : (maxField(2) + 15 .* PixelSizeIECg(2));
  iDoseGntZ = double( (Zdistal - 15 .* PixelSizeIECg(3) )  : PixelSizeIECg(3) : maxCEF       + 30 .* PixelSizeIECg(3));

end


%-------------------------------------------
% find the closest pixel coordinate to a given number
%-------------------------------------------
function [pxlVal , pixIdx] = getClosePixel(iDoseGntX , value)
  [~ , pixIdx] = min((iDoseGntX - value).^2);
  pxlVal = iDoseGntX(pixIdx);
end

%----------------------------------------------
% Align the field with the origin of the CEM and the CEM pixel resolution
%----------------------------------------------
function minField = alignFieldWithCEM(minField , RangeModulator)
  NbPxl = round( (minField - RangeModulator.ModulatorOrigin(1:2)) ./ RangeModulator.Modulator3DPixelSpacing(1:2)); %Number of pixels between origin of CEM and corner of the field
  minField = RangeModulator.ModulatorOrigin(1:2) + NbPxl .* RangeModulator.Modulator3DPixelSpacing(1:2); %Shift the origiin of the field to be aligned with CEM grid
end


%----------------------------------------------
% Create an empty sparse matrix
%This will store the dose influence matrix
% INPUT
%
% |Plan| -_STRUCTURE_- Structure defining a multi energy layer PBS plan
%
% |handles| -_STRUCTURE_- REggui data handle. The CT scan is stored in |handles| in the image with name |Plan.CTname|.
%
% OUTPUT
% |Pij| -_SCALAR MATRIX_- dose influence matrix: |Pij(vox,spot)| The dose contribution to voxel |vox| of the spot number |spot|
% |w| -_SCALAR VECTOR_- |w(spot)| weight per fraction of the j-th spot
%----------------------------------------------
function [Pij , w]  = preparePij(Plan , handles)

  NbBeamlets = 0;
  for b = 1:numel(Plan.Beams)
    NbBeamlets = NbBeamlets + numel(Plan.Beams(b).Layers(1).SpotWeights); %Number of PBS spots in plan with energy monolayer
  end
  nvoxels = prod(handles.size);

  Pij = sparse(nvoxels , NbBeamlets);  %|Pij(vox,spot)|. Create  zeros sparse matrix
          %Create empty sparse matrix with proper size
          %This will make reading faster because it avoids resizing of matrices

  [~ , w] = plan2weight(Plan); %The beamlet weight in the monolayer plan. This weight is the sum over all fractions

end

%------------------------------------
% Add one bemalet to the Pij dose influence matrix
%
% INPUT
% |Pij| -_SCALAR MATRIX_- dose influence matrix: |Pij(vox,spot)| The dose contribution to voxel |vox| of the spot number |spot|
% |spt| -_INTEGER_- Index of the beamlet to be added to Pij
% |DoseDCMcs| -_SCALAR MATRIX_- DoseDCMcs(x,y,z) Dose (Gy) at the pixel [x,y,z]. Matrix axes are aligned with the DICOM CS of the CT scan
% |Plan| -_STRUCTURE_- Structure defining a mono energy layer PBS plan
% |w| -_SCALAR VECTOR_- |w(spot)| weight per fraction of the j-th spot
%
% OUTPUT
% |Pij| -_SCALAR MATRIX_- [OPTIONAL] dose influence matrix: |Pij(vox,spot)| The dose contribution to voxel |vox| of the spot number |spot|
%------------------------------------
function Pij = addBeamlet2Pij(Pij , spt , jIndex , DoseSptDCMcs , Plan , w)

    DoseSptDCMcs = DoseSptDCMcs ./ w(spt) ; % w is Normalise the dose for w=1 and for 1 fraction
            %In SpotWeightsOptimization at line 85, the weight are divided by the number of fractions.

    temp = flip(DoseSptDCMcs,3);
    dose1D = temp(:); %spare matrix with the dose influence of the spot spt

    flagMasked = false;
    if isfield(Plan, 'OptROIVoxels_nominal')
        if (~isempty(Plan.OptROIVoxels_nominal) || ~isempty(find(Plan.OptROIVoxels_nominal,1)))
            dose1D = sparse(double(full(Plan.OptROIVoxels_nominal) .* dose1D)); %Apply the mask to force to zero the voxels outside of the RT struct of interrest. This saves memory
            %the .* product should use full matrices. Product .* with sparse matrices seems buggy in Matlab
            flagMasked = true;
        end
    end

    %If the matrix was not clipped with the mask, clip it with a dose threshold
    %This will save memory
    if max(dose1D,[],'all')
        %There is some dose in this bemalet. Add it to the full dose map
        if ~flagMasked
          thresh = max(dose1D,[],'all') ./ 1000; %Threshold at 0.1% of the maximum dose of the spot
          dose1D = sparse( (full(dose1D)>thresh) .* full(dose1D) );
                    %The product .* is applied to the full matrices
                    %when the product is appleid on the sparse matrices, the size of |dose1D| doubles instead of decreasinf (as there are more zeros).
                    %So the product is applied to double matrices. And then sparse is applied afterwards to remove all the zeros
        end
        Pij(:,jIndex) =  dose1D; %|Pij(vox,spot)|
    end
end


%-------------------------
%compute the deepest Z in IECg CS
% to which the dose should be computed
%
% OUTPUT
% |Zdistal| -_SCLAR_-  Z Coordinate (mm) in the IEC gantry CS of the deepest plane in which the dose is to be computed
%-------------------------
function Zdistal = getZdistal(PTV , Spacing , ImagePositionPatient , Beam)

  Bdcm = getDICOMcoord(PTV , Spacing , ImagePositionPatient , Beam.isocenter);
  M = matDICOM2IECgantry(Beam.GantryAngle,Beam.PatientSupportAngle);
  Bbev = M * Bdcm; %Coordinate of the voxels of the body in IC gantry
  Zdistal = min(Bbev(3,:)); %Position on beam axis of the distal surface of body
end


%------------------------------------------------
% Add the range shifter in the high resolution Ct scan
%This must be done in the high resolution CT to avoid the RS thickness to be aliased by the Z pixel resolution of the CT
%As the CT is aligned with the IEC gantry CS, adding the RS can be done quickly by using the Matlab indexing.
%------------------------------------------------
function handlesHR = setRangeShifterinHRCT(handlesHR , PlanHR , hrCTName)

  global g_HUrangeshifter;

  CT = Get_reggui_data(handlesHR,hrCTName,'images'); %Update the CT scan with the aperture block in handles

  if isfield( PlanHR.Beams , 'RSinfo')

      for slab =  PlanHR.Beams.RSinfo.NbSlabs : -1 : 1
          fprintf('Slab # %d \n',slab)

          %%Compute Zg of upstream side of range shifter
          ZgUp =  PlanHR.Beams.RSinfo.IsocenterToRangeShifterDistance +  ... %Donwstream side of slab close to isocenter
                 PlanHR.Beams.RSinfo.SlabOffset(slab) - ...  %distance from downstream side of 1st slab to upstream side of |slab|
                 handlesHR.spacing(3) ./2 ; %Remove half a pixel to getthe coordinate of the center of last pixel

          %+Zg is aligned with -Y CT
          [ ~, ~ , Zup ]= DICOM2PXLindex([] , handlesHR.spacing , handlesHR.origin , true, 0 , -ZgUp , 0 );
          ZgDwn = ZgUp -  PlanHR.Beams.RSinfo.RSslabThickness(slab) + handlesHR.spacing(3) ./2; %Add half a pixel to get coordinate of center of first pixel
          [ ~, ~ , Zdn ]= DICOM2PXLindex([] , handlesHR.spacing , handlesHR.origin , true, 0 , -ZgDwn , 0 );

          stp = 1 .* sign(Zup-Zdn) ;
          for slice = Zdn:stp:Zup
            CT(:,slice,:) = g_HUrangeshifter;
          end % for slice
      end %for slab
  end % if isfield( PlanHR.Beams , 'RSinfo')

  handlesHR = Set_reggui_data(handlesHR , hrCTName , CT ,  PlanHR.CTinfo , 'images',1); %Update the CT scan with the aperture block in handles


end


%------------------------------------------------
% Add the aperture in the high resolution Ct scan
%------------------------------------------------
function handlesHR = setApertureinHRCT(handlesHR , PlanHR , hrCTName)

  global g_HUbrass;

  CT = Get_reggui_data(handlesHR,hrCTName,'images'); %Update the CT scan with the aperture block in handles
        %The proton beam is aligned with the **Y axis** of the CT scan : Xg = Xct; Yct = -Zg ; Zct = Yg
  CTsize = size(CT);
  ImagePositionPatient = handlesHR.origin;

  if PlanHR.Beams.ApertureBlock

        Abrass  = apertureContour2Block(PlanHR.Beams ,  handlesHR.spacing  , PlanHR.BDL); %Coordinate of the brass voxels in IEC gantry

        %Xg = Xct; Yct = -Zg ; Zct = Yg
        MinC = ImagePositionPatient;
        MaxC = ImagePositionPatient + (CTsize'-1) .* handlesHR.spacing;

        %                      Xg = Xct                    Yg = Zct                 Zg  = -Yct                 Xg = Xct                    Yg = Zct                 Zg  = -Yct (inverse min and max because of - sign)
        GoodIdx = find((Abrass(:,1) > MinC(1)) .* (Abrass(:,2) > MinC(3)) .* (Abrass(:,3) > -MaxC(2)) .* (Abrass(:,1) <= MaxC(1)) .* (Abrass(:,2) <= MaxC(3)) .* (Abrass(:,3) <= -MinC(2)) ); %This is the list of indices that fit in the CT
        [~, X , Y , Z] = DICOM2PXLindex([] , handlesHR.spacing , ImagePositionPatient , 1, Abrass(GoodIdx,1) , -Abrass(GoodIdx,3) , Abrass(GoodIdx,2));
        Aidx = sub2ind(size(CT) , X , Y , Z); %Only get the indices of the aperture that fits in the CT scan
        CT(Aidx) = g_HUbrass; %Put Hu of the device in the voxels of the device

  end

  %Update the handles and plan
  handlesHR = Set_reggui_data(handlesHR,hrCTName,CT,PlanHR.CTinfo,'images',1); %Create a new image with only the aperture

end


%------------------------------------
% Create a new REGGUI handles with the apstial properties of the image
%
% INPUT
% |size| -_SCALAR VECTOR_- [Nx , Ny , Nz] Number of voxel along each axis
% |origin|  -_SCALAR VECTOR_- [x,y,z] coordiante of the first voxel of the images
% |spacing| -_SCALAR VECTOR_- [dx , dy , dz] Dimension (mm) of each axis of the voxels
%------------------------------------
function handles = createREGGUIhandles(size , origin , spacing)
  handles = struct;
  handles.size = size; %Make first a copy of the size of the original CT scan. We want the dose map to match this size
  handles.origin = origin;
  handles.spacing = spacing;
  handles.spatialpropsettled = true;
end
