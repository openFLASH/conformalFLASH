%% insertDeviceInCT
% Virtually insert a beam line device into a CT scan by changing the Hounsfield units of the voxels defining the device
% If the CT scan is too small to fit the device at the requested position, the CT scan is expanded so as to fit the device
%
%
%% Syntax
% |[CT , ImagePositionPatient] = insertDeviceInCT(CT , Adevice, HUdevice , Beam , Spacing , ImagePositionPatient , HUair)|
%
%
%% Description
% |[CT , ImagePositionPatient] = insertDeviceInCT(CT , Adevice, HUdevice , Beam , Spacing , ImagePositionPatient , HUair)| Description
%
%
%% Input arguments
% |CT| - _SCALAR MATRIX_ - CT scan in which the aperture is to be added |CT(x,y,z)=HU| gives the intensity (Hounsfield unit) of the voxel at location (x,y,z)
%
% |Adevice| -_SCALAR MATRIX_- |A(i,:) = [x,y,z,1]| Matrix of 4 vector [x,y,z,1] defining the position of the i-th voxel ( mm, in IC gantry CS) of the voxels belonging to the device
%
% |HUdevice| -_SCALAR_- Hounsfield unit of the device to be added to the CT scan
%
% |Beams| -_STRUCTURES_- Information about the beam
%     * |Beams.GantryAngle| -_SCALAR_- Gantry Angle (deg)
%     * |Beams.PatientSupportAngle| -_SCALAR_- Couch Angle (deg)
%     * |Beams.isocenter| -_SCALAR VECTOR_- [x,y,z] Coordiantes (mm) of the isocentre in the planning CT scan
%
% |Spacing| - _SCALAR VECTOR_ - Pixel size (|mm|) of the displayed images in GUI
%
% |ImagePositionPatient| - _SCALAR VECTOR_ - Coordinate (in |mm|) of the first pixel of the image in the coordinate system of the image
%
% |HUair| - _SCALAR VECTOR_ - If the CT must be expanded, the new voxels are filled with |HUair|
%
%% Output arguments
%
% |CT| - _SCALAR MATRIX_ - Updated CT scan
%
% |ImagePositionPatient| - _SCALAR VECTOR_ - Updated origin of the Ct scan if the Ct was expanded
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [CT , ImagePositionPatient] = insertDeviceInCT(CT , Adevice, HUdevice , Beam , Spacing , ImagePositionPatient , HUair)

  Aidx = [];
  while isempty(Aidx)
    [Aidx , Axyz]= IECgantry2CTindex(Adevice, Beam , Spacing , ImagePositionPatient , size(CT) , 2); %Find pixels indices in the CT scan coordinate system

    if isempty(Aidx)
      %The CT is too small to fit the aperture. Make the CT larger
      fprintf('CT too small to fit object. Expanding CT \n')
      minA = min(min(Axyz,[],1) , [1,1,1] ); %The smallest index
      maxA = max(max(Axyz,[],1) , size(CT)); %The largest index
      sCT2 = maxA - minA + 1; % The minimum size of the CT
      ImagePositionPatientOLD =  ImagePositionPatient;
      ImagePositionPatient = ImagePositionPatient - abs(minA -1)' .* Spacing;
      CT2 = enlargeCT(CT , ImagePositionPatientOLD , ImagePositionPatient , Spacing , sCT2 , HUair);
      CT = CT2; %This is the expanded CT

    end
  end

  CT(Aidx) = HUdevice; %Put Hu of the device in the voxels of the device

end
