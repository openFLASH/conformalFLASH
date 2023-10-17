%% setCEMinhrCT
% Add the Conformal Energy Modulator (CEM) in the high resolution CT scan for the single beams with field |CEM3Dmask|.
% The high resolution CT has axes aligned with the IEC gantry CT scan: Xg = Xct; Yct = -Zg ; Zct = Yg
% The high reoslution CT must have the same voxel size in Xct and Zct as the CEM. It must be sufficiently long in Yct to inser the CEM.
%
% The cross section width of the CEM is defined by the size of the  CT scan or matrix |CEM3Dmask|, whatever is smaller.
%
% The CT scan is the image with name |CTname| stored in |handles|
% Save the updated CT in the |handles| with name |CTname|
%
%% Syntax
% |handles = setCEMinhrCT(handles , Plan , CTname , HUcem , HUair)|
%
%
%% Description
% |handles = setCEMinhrCT(handles , Plan , CTname , HUcem , HUair)| Description
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
%        * |Plan.Beams.addPrinterErrors| -_STRING_- Type of error to introduce in the CEM: false , 'dilate' , 'erode' , 'TallSlim','rough'
%        * |Plan.Beams.RangeModulator.CEM3Dmask| -_SCALAR MATRIX_- 3D mask of the CEM. |CEM3Dmask(x,y,z)=1| if the voxel at location (x,y,z)  in the plane of the CEM for beam b belongs to the CEM.
%                                               Z=0 at the base of CEM. Z increase in the same way as Zg if the spike point toward the proton source
%        * |Plan.Beams.RangeModulator.Modulator3DPixelSpacing| -_SCALAR VECTOR_- |CompensatorPixelSpacing = [x,y,z]| Pixel size (mm) in the plane of the CEM
%        * |Plan.Beams.RangeModulator.ModulatorMountingPosition| -_STRING_- Direction in which the CEM is pointing from the tray. 'PATIENT_SIDE' or 'SOURCE_SIDE'
%        * |Plan.Beams.RangeModulator.IsocenterToRangeModulatorDistance| -_SCALAR_- Distance (mm) from isocentre to the base of the CEM.
%
% |CTname| -_STRING_- Name of the CT image in handles.images into which the device is to be inserted
%
% |HUcem| -_SCALAR_- [OPTIONAL. If absent, read from disk] Hounsfield unitof the  CEM to insert in the CT scan
%
% |HUair| -_SCALAR_- [OPTIONAL. If absent, read from disk] Hounsfield unit of air air to insert in CT scan
%
%% Output arguments
%
% |handles| -_STRUCTURE_- Updated REggui data handle. The CT scan is stored in |handles| in the image with name |Plan.CTname|.
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function handles = setCEMinhrCT(handles , Plan , CTname , HUcem , HUair)


  CT = Get_reggui_data(handles, CTname ,'images'); %Update the CT scan with the aperture block in handles
  ImagePositionPatient = handles.origin;

  if sum(Plan.Beams.RangeModulator.Modulator3DPixelSpacing(1:2) - [handles.spacing(1),handles.spacing(3)])
    error('High resolution CT does not have same pixel resolution as CEM')
  end


  if isfield(Plan.Beams, 'addPrinterErrors')
    % ROBUSTESS TEST
    % Dilate the CEM shape to estimate the effect of the inacuracies of the printer that dilate the spikes of the CEM
    CEM = addPrinterErrors(Plan.Beams.RangeModulator.CEM3Dmask , Plan.Beams.addPrinterErrors); % Type of error to introduce: false , 'dilate' , 'erode' , 'TallSlim'
  else
    %No definition of the printer error in the YAML file. Do not add any error
    CEM = Plan.Beams.RangeModulator.CEM3Dmask;
  end

  %Add the CEM in the CT scan
  %------------------------------------
  switch Plan.Beams.RangeModulator.ModulatorMountingPosition
    case 'SOURCE_SIDE'
      %The Z_cem is pointing in the same direction as the Zg
      Origin = Plan.Beams.RangeModulator.ModulatorOrigin + [0,0, Plan.Beams.RangeModulator.IsocenterToRangeModulatorDistance];

    case  'PATIENT_SIDE'
      Origin = [Plan.Beams.RangeModulator.ModulatorOrigin(1:2) , 0];

      %The Z_cem is pointing opposite to the Zg
      CEM = flip(CEM,3); %flip the matrix to have it ordered in increasing Zg
      Origin(3) = Plan.Beams.RangeModulator.IsocenterToRangeModulatorDistance - ( Plan.Beams.RangeModulator.ModulatorOrigin(3) + size(CEM,3) .* Plan.Beams.RangeModulator.Modulator3DPixelSpacing(3) );
          %the origin of the 3D image is now at the tip of the CEM
  end

  %Re-interpolate the CEM at the resolution of the CT scan
  M = matDICOM2IECgantry(Plan.Beams.GantryAngle , Plan.Beams.PatientSupportAngle , Plan.Beams.isocenter); %Rotate around isocentre
  ZgVec = M * [0,0,0,1 ; 0,0,Plan.Beams.RangeModulator.Modulator3DPixelSpacing(3),1]'; %Define the orientation and length the Zg vector in the CT scan CS
  ZgVec = ZgVec(:,2) -  ZgVec(:,1);  %Define the orientation and length the Zg vector in the CT scan CS
  ZgVec = round(ZgVec(1:3) ./ handles.spacing); %Number of CT scan pixels for each Z step in CEM
  ResFac = max(ZgVec);

  if ResFac > 1
    %Make CEM pixels smaller in Z to match CT scan reoslution
    %Otherwise there will be missing layers in CEM
    fprintf('resizing \n')
    CEM = imresize3(CEM , size(CEM) .* [1,1,ResFac] , 'nearest');
    CEMres = Plan.Beams.RangeModulator.Modulator3DPixelSpacing';
    CEMres(3) = CEMres(3) ./ ResFac;
  else
    %Keep current CEM resolution
    CEM = CEM;
    CEMres = Plan.Beams.RangeModulator.Modulator3DPixelSpacing';
  end

  % Rotate the CEM axes to be aligned with CT axes
  % The high resolution CT has axes aligned with the IEC gantry CT scan: Xg = Xct; Yct = -Zg ; Zct = Yg
  CEM = permute(CEM, [1,3,2]);
  CEM = flipdim(CEM,2);
  CEMres = CEMres([1 3 2 ]);
  Origin = Origin([1 3 2 ]);
  Origin(2) = -(Origin(2) + size(CEM,2) .* CEMres(2));


  %Get the Yct position in the CT where to insert the CEM
  %Coordinate along the beam axis
  [ ~, ~ , Yct1 ]= DICOM2PXLindex([] , handles.spacing , handles.origin , true, 0 , Origin(2) , 0 );
  [ ~, ~ , Yct0 ]= DICOM2PXLindex([] , handles.spacing , handles.origin , true, 0 , Origin(2) , 0 );
  Yct1 = Yct0 + size(CEM,2) - 1;


  %Clip the CEM at the edge of the CT
  %Coordiante in a plane normal to beam axis
  [ ~, Xg0 , ~ , Yg0 ]= DICOM2PXLindex([] , CEMres , Origin , true, handles.origin(1) , handles.origin(2) , handles.origin(3) );
  Xg1 = Xg0 + size(CT,1) - 1;  %If the first one is #1 and there are 10 elements, the last one is #10 i.e. 1 + 10 -1
  Yg1 = Yg0 + size(CT,3) - 1;

  %Lateral pixel coordinates of the CT where to insert the CEM
  Xct0 = 1;
  Zct0 = 1;
  [Xct1 , ~ , Zct1] = size(CT);


  %Do we have enough pixels in the CEM to fit inside the CT ?
  if Xg0 < 1
    %Part of the Ct scan is outside of the CEM
    offSet = 2 - Xg0; %Number of voxels of CT scan out of CEM
    Xg0 = 1;
    Xct0 = offSet; %The CEM is not inserted at edge of CT scan
  end

  if Xg1 > size(CEM,1)
    %Part of the Ct scan is outside of the CEM
    offSet = Xg1 - size(CEM,1);   %Number of voxels of CT scan out of CEM
    Xg1 = size(CEM,1);
    Xct1 = Xct1 - offSet;
  end

  if Yg0 < 1
    %Part of the Ct scan is outside of the CEM
    offSet = 2 - Yg0; %Number of voxels of CT scan out of CEM
    Yg0 = 1;
    Zct0 = offSet; %The CEM is not inserted at edge of CT scan
  end

  if Yg1 > size(CEM,3)
    %Part of the Ct scan is outside of the CEM
    offSet = Yg1 - size(CEM,3);  %Number of voxels of CT scan out of CEM
    Yg1 = size(CEM,3);
    Zct1 = Zct1 - offSet;
  end

  %Insert the range compensator in the CT
  % The mask is 0 or 1.
  % 0 is replaced by HUair and 1 is replaced by HUcem
  CT(Xct0:Xct1 , Yct0:Yct1 , Zct0:Zct1) = CEM(Xg0:Xg1 , : , Yg0:Yg1) .* (HUcem - HUair) + HUair;

  %Update the handles with new CT
  handles = Set_reggui_data(handles , CTname , CT , Plan.CTinfo , 'images',1); %Update the CT scan with the aperture block in handles


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
