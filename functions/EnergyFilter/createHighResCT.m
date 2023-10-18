%% createHighResCT
% Create a high resolution CT with the **Y axis** aligned on the proton beam axis
% from a low resolution CT : Xg = Xct; Yct = -Zg ; Zct = Yg
% |minField| and |maxField| define the beam)s eye view dimension of the high resolution CT.
% The position of the CEM block and isocentre define the length of the high resolution CT scan
% The high res CT is rotated the interpolated CT so that its the Zg axis points towards the Ydicom axis
%
%% Syntax
% |[handlesHR , BeamHR , iCTgntY , iCTgntX , iCTgntZ] = createHighResCT(handles , lrCTname , hrCTname , Beam , PixelSize , HUair , minField , maxField , Zdistal , CTinfo)|
%
%
%% Description
% |[handlesHR , BeamHR , iCTgntY , iCTgntX , iCTgntZ] = createHighResCT(handles , lrCTname , hrCTname , Beam , PixelSize , HUair , minField , maxField , Zdistal , CTinfo)| Description
%
%
%% Input arguments
%
% |handles| -_STRUCTURE_- REggui data handle. The CT scan is stored in |handles|
%
% |lrCTname| -_STRING_- Name of the low resolution CT image in handles.images
%
% |hrCTname| -_STRING_- Name of the high resolution CT image to be saved in handles.images
%
% |Beam| -_STRUCTURES_- Information about the beam
%     * |Beam.RangeModulator.CEM3Dmask| -_SCALAR MATRIX_- 3D mask of the CEM. |CEM3Dmask(x,y,z)=1| if the voxel at location (x,y,z)  in the plane of the CEM for beam b belongs to the CEM.
%                                 Z=0 at the base of CEM. Z increase in the smae way as Zg if the spike point toward the proton source
%     * |Beam.RangeModulator.ModulatorMountingPosition| -_STRING_- Direction in which the CEM is pointing from the tray. 'PATIENT_SIDE' or 'SOURCE_SIDE'
%     * |Beam.RangeModulator.IsocenterToRangeModulatorDistance| -_SCALAR_- Distance (mm) from isocentre to the base of the CEM.
%     * |Beam.RangeModulator.Modulator3DPixelSpacing| -_SCALAR VECTOR_- |CompensatorPixelSpacing = [x,y,z]| Pixel size (mm) in the plane of the CEF for the |CompensatorThicknessData| matrix
%
% |PixelSize| -_SCALAR VECTOR_- The [x,y,z] pixel size (mm) of the re-interpolated CT scan
%
% |HUair| - _SCALAR VECTOR_ - If the CT must be expanded, the new voxels are filled with |HUair|

% |minField| -_SCALAR VECTOR_- [OPTIONAL, only needed if part of the CEM is to be used] [X, Y] Coordinate (mm, in IEC gantry) of [-x,-y] the corner of the field
%
% |maxField| -_SCALAR VECTOR_- [OPTIONAL, only needed if part of the CEM is to be used] [X, Y] Coordinate (mm, in IEC gantry) of [+x,+y] the corner of the field
%
% |Zdistal| -_SCLAR_-  Z Coordinate (mm) in the IEC gantry CS of the deepest plane in which the dose is to be computed
%
% |CTinfo| -_STRUCTURE_- DICOM header of the CT scan, to be added to the image in |handles|
%
%
%% Output arguments
%
% |handles| -_STRUCTURE_- Updated REggui data handle. The CT scan is stored in |handles|
%
% |BeamHR| -_STRUCTURE_- Updated beam structure with information about the high resolution CT scan
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [handlesHR , BeamHR , iCTgntY , iCTgntX , iCTgntZ] = createHighResCT(handles , lrCTname , hrCTname , Beam , PixelSize , HUair , minField , maxField , Zdistal , CTinfo)

  CT = Get_reggui_data(handles, lrCTname); %Update the CT scan with the aperture block in handles
  Spacing = handles.spacing;
  ImagePositionPatient = handles.origin;

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

  % Round the CT scan size to a multiple of the pixel size
  maxCEF = rounding(maxCEF , PixelSize(3));
  Zdistal = rounding(Zdistal , PixelSize(3));

  %Define the coordinate of the voxels in the IEC gantry CS for the interpolated Ct scan
  iCTgntX = minField(1) : PixelSize(1) : maxField(1) ; %Coordinate axes of the interpolated CT expressed in IEC gantry CS
  iCTgntY = minField(2) : PixelSize(2) : maxField(2) ;
  iCTgntZ = double(Zdistal : PixelSize(3) : maxCEF + 10 .* PixelSize(3));

  %Compute the coordinate of the voxels of the interpolated CT scan in the DICOM cs
  %Origin of the interpolated DICOM CS is at the isocentre
  [Y , X , Z] = meshgrid( iCTgntY , iCTgntX , iCTgntZ);     %meshgrid uses index order Xct(j,i,k) where CT(i,j,k). We must swap the X and Y axis to get the index ordering correctly linked with CT
  A = [X(:), Y(:), Z(:) , ones(numel(X),1)]; %Coordinate of all voxels of the interpolated CT scan, expressed in IEC gantry

  BeamINTER = Beam;
  BeamINTER.isocenter = [0,0,0]'; %The isocenter is at the origin of the coordinate system of the high res CT scan and the Z axis is aligned with the Z IEC  gantry CS
  iCTdcm = gantry2DICOM(A' , BeamINTER); %Coordinate of all voxels of interpolated CT expressed in DICOM CS

  %Compute the coordinate of the original CT in the DICOM CS
  [oCTdcmX, oCTdcmY , oCTdcmZ , X , Y , Z] = getCTaxes(handles.origin , handles.spacing , handles.size , Beam.isocenter);

  %Re-interpolate the original CT scan on a grid with smaller pixel size
  %and which is sufficiently extended so as to include aperture and ridge filter
  %The origin of the 2CT scans is at isocentre
  %The pixels axes are paralell to the IEC gantry CS
  CTintrp = interpolateCT(CT, oCTdcmX, oCTdcmY, oCTdcmZ , iCTdcm(1,:), iCTdcm(2,:), iCTdcm(3,:) , HUair , [numel(iCTgntX) , numel(iCTgntY) , numel(iCTgntZ)] );

  %Rotate the interpolated CT so that its the Zg axis points towards the Ydicom axis
  %In this way, we will shoot with a virtual gantry at 0° into the interpolated CT
  %Xg = Xct; Yct = -Zg ; Zct = Yg
  CTintrp = permute(CTintrp,[1,3,2]);
  CTintrp = flip(CTintrp,2); %NB: flip (C code) is faster than flipdim (.m).

  %Recompute the orign and spacing of the re-interpolated CT scan
  handlesHR = struct;
  handlesHR.path = handles.path;
  handlesHR.dataPath = handles.dataPath;
  handlesHR = Initialize_reggui_handles(handlesHR); % Initialize the handles structure
  handlesHR.origin = [iCTgntX(1) , -max(iCTgntZ) , iCTgntY(1) ]'; %Second index is minus (because flipped) Zg
  handlesHR.spacing = [PixelSize(1) , PixelSize(3) , PixelSize(2)]' ;
  handlesHR.size = size(CTintrp);
  handlesHR.isocenter = BeamINTER.isocenter; %The beam isocenter is at the origin of the CT scan
  handlesHR.isocenter = [handlesHR.isocenter(1) , handlesHR.isocenter(3) , handlesHR.isocenter(2)]; %permuted the CT scan
  handlesHR.spatialpropsettled = true; %This flag is set to force the re-interpolation of the dose map when loading it

  CTinfo.ImagePositionPatient = handlesHR.origin; %Update origin of the CT scan
  CTinfo.Spacing = handlesHR.spacing;

  handlesHR = Set_reggui_data(handlesHR , hrCTname , CTintrp , CTinfo , 'images');

  BeamHR.isocenter = handlesHR.isocenter;
  BeamHR.GantryAngle = 0;
  BeamHR.PatientSupportAngle = 0;

end


%========================
% Convert point coordonate from IEC gantry to DICOM CS
%========================
function Adcm = gantry2DICOM(A , Beam)

    M = matDICOM2IECgantry(Beam.GantryAngle,Beam.PatientSupportAngle);
    Adcm = inv(M) * A ; %Coordinate of voxels in DICOM CS (original CT)
    %The rotations took place around the isocentre which is at origin of the coordinate system of CT scan

end
