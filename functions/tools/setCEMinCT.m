%% setCEMinCT
% Add the Conformal Energy Modulator (CEM) in the CT scan for all beams with field |CEM3Dmask|.
% The cross section width of the CEM is defined by the size of the matrix |CEM3Dmask|.
% If |minField| and |maxField| are provided, only the cross section defined by these parameter is inserted in the CT.
% The CT scan is the image with name |CTname| stored in |handles|
% Expand the CT if it is too small.
% Save the updated CT in the |handles| with name |CTname|
% Update the |Plan| structure with the updated CT size
% Update all the images in |handles| with the updated CT size
%
%% Syntax
% |[Plan , handles ] = setCEMinCT(handles , Plan , CTname )|
%
% |[Plan , handles ] = setCEMinCT(handles , Plan , CTname , minField , maxField )|
%
% |[Plan , handles ] = setCEMinCT(handles , Plan , CTname , minField , maxField , HUcem )|
%
% |[Plan , handles ] = setCEMinCT(handles , Plan , CTname , minField , maxField , HUcem , HUair)|
%
%% Description
% |[Plan , handles ] = setCEMinCT(handles , Plan , CTname , minField , maxField , HUcem , HUair)| Description
%
%
%% Input arguments
%
% |handles| -_STRUCTURE_- REggui data handle. The CT scan is stored in |handles| in the image with name |Plan.CTname|.
%
% |Plan| -_STRUCTURE_- Structure defining a multi energy layer PBS plan
%     * |Plan.Spike.MaterialID| - _STRING_ - Name of the CEM material, as defined in the file "plugins\openMCsquare\lib\Materials\list.dat"
%     * |Plan.ScannerDirectory| -_STRING_- Ct scanner HU to SPR conversion parameters. Name of the folder in REGGUI\plugins\openMCsquare\lib\Scanners
%     * |Plan.DoseGrid| - _struct_ - Structure containing the information about the dose grid. The following data must be present:
%     * |Plan.DoseGrid.nvoxels| - _scalar_ - Total number of voxels of the dose grid, i.e., number of voxels in the CT image
%     * |Plan.CTinfo| -_STRUCTURE_- DICOM header of the CT scan
%     * |Plan.Beams| -_STRUCTURES_- Information about the beam
%        * |Plan.Beams.GantryAngle| -_SCALAR_- Gantry Angle (deg)
%        * |Plan.Beams.PatientSupportAngle| -_SCALAR_- Couch Angle (deg)
%        * |Plan.Beams.isocenter| -_SCALAR VECTOR_- [x,y,z] Coordiantes (mm) of the isocentre in the planning CT scan
%        * |Plan.Beams(b).addPrinterErrors| -_STRING_- Type of error to introduce in the CEM: false , 'dilate' , 'erode' , 'TallSlim','rough'
%        * |Plan.Beams(b).RangeModulator.CEM3Dmask| -_SCALAR MATRIX_- 3D mask of the CEM. |CEM3Dmask(x,y,z)=1| if the voxel at location (x,y,z)  in the plane of the CEM for beam b belongs to the CEM.
%                                               Z=0 at the base of CEM. Z increase in the same way as Zg if the spike point toward the proton source
%        * |Plan.Beams(b).RangeModulator.Modulator3DPixelSpacing| -_SCALAR VECTOR_- |CompensatorPixelSpacing = [x,y,z]| Pixel size (mm) in the plane of the CEM
%        * |Plan.Beams(b).RangeModulator.ModulatorMountingPosition| -_STRING_- Direction in which the CEM is pointing from the tray. 'PATIENT_SIDE' or 'SOURCE_SIDE'
%        * |Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance| -_SCALAR_- Distance (mm) from isocentre to the base of the CEM.
%
% |CTname| -_STRING_- Name of the CT image in handles.images into which the device is to be inserted
%
% |minField| -_SCALAR VECTOR_- [OPTIONAL, only needed if part of the CEM is to be used] [X, Y] Coordinate (mm, in IEC gantry) of [-x,-y] the corner of the field
%
% |maxField| -_SCALAR VECTOR_- [OPTIONAL, only needed if part of the CEM is to be used] [X, Y] Coordinate (mm, in IEC gantry) of [+x,+y] the corner of the field
%
% |HUcem| -_SCALAR_- [OPTIONAL. If absent, read from disk] Hounsfield unitof the  CEM to insert in the CT scan
%
% |HUair| -_SCALAR_- [OPTIONAL. If absent, read from disk] Hounsfield unit of air air to insert in CT scan
%
%% Output arguments
%
% |handles| -_STRUCTURE_- Updated REggui data handle. The CT scan is stored in |handles| in the image with name |Plan.CTname|.
%
% |Plan| -_STRUCTURE_- Updated structure defining a multi energy layer PBS plan
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [Plan , handles ] = setCEMinCT(handles , Plan , CTname , minField , maxField , HUcem , HUair)

  if nargin < 4
    minField = [];
    maxField = [];
  end

  if nargin < 6
    HUcem = getMaterialSPR(Plan.Spike.MaterialID , Plan.ScannerDirectory) +1 ; %Hounsfield unit associated to CEM in the material file
  end
  if nargin < 7
    HUair =  getMaterialSPR('Schneider_Air' , Plan.ScannerDirectory) +1 ; %Hounsfield unit associated to air in the material file
  end

  CT = Get_reggui_data(handles, CTname ,'images'); %Update the CT scan with the aperture block in handles
  ImagePositionPatient = handles.origin;


  for b = 1: size(Plan.Beams,2) %Loop for each beam

    if isfield(Plan.Beams(b), 'addPrinterErrors')
      % ROBUSTESS TEST
      % Dilate the CEM shape to estimate the effect of the inacuracies of the printer that dilate the spikes of the CEM
      CEM = addPrinterErrors(Plan.Beams(b).RangeModulator.CEM3Dmask , Plan.Beams(b).addPrinterErrors); % Type of error to introduce: false , 'dilate' , 'erode' , 'TallSlim'
    else
      %No definition of the printer error in the YAML file. Do not add any error
      CEM = Plan.Beams(b).RangeModulator.CEM3Dmask;
    end

    %Add the CEM in the CT scan
    %------------------------------------
    switch Plan.Beams(b).RangeModulator.ModulatorMountingPosition
      case 'SOURCE_SIDE'
        %The Z_cem is pointing in the same direction as the Zg
        Origin = Plan.Beams(b).RangeModulator.ModulatorOrigin + [0,0, Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance];

      case  'PATIENT_SIDE'
        Origin = [Plan.Beams(b).RangeModulator.ModulatorOrigin(1:2) , 0];

        %The Z_cem is pointing opposite to the Zg
        CEM = flip(CEM,3); %flip the matrix to have it ordered in increasing Zg
        Origin(3) = Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance - ( Plan.Beams(b).RangeModulator.ModulatorOrigin(3) + (size(CEM,3) - 1) .* Plan.Beams(b).RangeModulator.Modulator3DPixelSpacing(3) );
            %the origin of the 3D image is now at the tip of the CEM
    end


    %Re-interpolate the CEM at the resolution of the CT scan
    M = matDICOM2IECgantry(Plan.Beams(b).GantryAngle , Plan.Beams(b).PatientSupportAngle , Plan.Beams(b).isocenter); %Rotate around isocentre
    ZgVec = M * [0,0,0,1 ; 0,0,Plan.Beams(b).RangeModulator.Modulator3DPixelSpacing(3),1]'; %Define the orientation and length the Zg vector in the CT scan CS
    ZgVec = ZgVec(:,2) -  ZgVec(:,1);  %Define the orientation and length the Zg vector in the CT scan CS
    ZgVec = round(ZgVec(1:3) ./ handles.spacing); %Number of CT scan pixels for each Z step in CEM
    ResFac = max(ZgVec);

    if ResFac > 1
      %Make CEM pixels smaller in Z to match CT scan reoslution
      %Otherwise there will be missing layers in CEM
      fprintf('resizing \n')
      CEMi = imresize3(CEM , size(CEM) .* [1,1,ResFac] , 'nearest');
      CEMres = Plan.Beams(b).RangeModulator.Modulator3DPixelSpacing;
      CEMres(3) = CEMres(3) ./ ResFac;
    else
      %Keep current CEM resolution
      CEMi = CEM;
      CEMres = Plan.Beams(b).RangeModulator.Modulator3DPixelSpacing;
    end

    %Get the coordinates of the voxels of the CEM in the IEC gantry CS
    Acem = getDICOMcoord(CEMi, CEMres , Origin , [0,0,0]);
    Acem = Acem';
    Acem = Acem(:,1:3);

    if ~isempty(minField) & ~isempty(maxField)
          %Remove the points outside of the field of interest
          Acem( Acem(:,1) < minField(1) ,:) = [];
          Acem( Acem(:,2) < minField(2) ,:) = [];
          Acem( Acem(:,1) > maxField(1) ,:) = [];
          Acem( Acem(:,2) > maxField(2) ,:) = [];
    end

    Acem = [Acem , ones(size(Acem,1),1)];

    %Insert the range compensator in the CT
    [CT , ImagePositionPatient] = insertDeviceInCT(CT , Acem, HUcem , Plan.Beams(b) , handles.spacing , ImagePositionPatient, HUair);

  end %for b

  %Update the handles and plan
  [handles , Plan] = updateAllImages(handles , Plan , CT , ImagePositionPatient , HUair  , CTname);

end



%=================================================
% Add errors to the 3D CEM shape
% to simulate the effect of the inacuracies of the printer that dilate the spikes of the
%
% INPUT
% |hedgehog| -_SCALAR MATRIX_- |hedgehog(x,y,z) = 1| if the pixel at (x,y,z) belongs to the hedgehog
% |typeEr| -_STRING_- Type of error to introduce: 'dilate' , 'erode' , 'TallSlim'
%=================================================
function hedgehog = addPrinterErrors(hedgehog , typeEr)

  switch typeEr
    case 'dilate'
        %Uniform dilate in all direction
        NbPxl = 1; %Number of pixels of the dilation. This gives the magnitude of the printer error
        fprintf('Adding printing error. Type : %s with %d pixels \n', typeEr,NbPxl)
        se = strel('sphere',NbPxl);
        hedgehog = imdilate(hedgehog,se);

    case 'erode'
        NbPxl = 1; %Number of pixels of the dilation. This gives the magnitude of the printer error
        fprintf('Adding printing error. Type : %s with %d pixels \n', typeEr,NbPxl)
        se = strel('sphere',NbPxl);
        hedgehog = imerode(hedgehog,se);

    case 'TallSlim'
        %The printer makes the spikes slimer and taller than requested
        NbPxl = 2; %Number of pixels of the dilation. This gives the magnitude of the printer error
        fprintf('Adding printing error. Type : %s with %d pixels \n', typeEr,NbPxl)
        se = strel('sphere',NbPxl);
        hedgehog = imerode(hedgehog,se); %Erode in all direction

        k = zeros(3,3,3); %Kernel to expand in the Z direction
        k(2,2,2)=1;
        k(2,2,1)=1; %The isocentre is located towards z =0
        k(2,2,3)=1;
        se = strel(k);
        %Cancel the erosion in the Z direction
        for i=1:NbPxl
          hedgehog = imdilate(hedgehog,se); %make one pixel taller towards isocentre
        end

        k = zeros(3,3,3); %Kernel to expand in the Z direction
        k(2,2,2)=1;
        k(2,2,1)=1; %The isocentre is located towards z =1
        %k(2,2,3)=1; %TODO when the spike go towards the source, we should expend towards Z=3
        %Dilate the height of the spike towards the isocentre by 3 pixels
        se = strel(k);
        for i=1:3
          hedgehog = imdilate(hedgehog,se); %make one pixel taller towards isocentre
        end

    case 'rough'
        %The printer add chunk of material at different places on the spikes
        BlisterSize = 0.6; % Fraction of potential 'blister' voxels that will become actual rough surface of the spike
        NbPxl = 2; %Number of pixels of the dilation. This gives the magnitude of the printer error
        fprintf('Adding printing error. Type : %s with %d pixels and fraction %f \n', typeEr,NbPxl,BlisterSize)
        se = strel('sphere',NbPxl);
        blister = imdilate(hedgehog,se) - hedgehog; %These are all the voxels where a blister could be located
        blister(blister < 0) = 0; %Make sure there is no negative number. There shouldn't be any.
        Pidx = find(blister); %Find the location of all the potential blister voxels
        Pidx = Pidx .* (rand(size(Pidx)) <= BlisterSize); %Randomly select the voxels that are not part of the blister
        Pidx(Pidx==0)=[]; %REmove the zeros from the vecotr
        hedgehog(Pidx) = 1; %Add the rough surface to the hedgehog

     otherwise
      fprintf('No printer error added \n')


  end

end
