%% IECgantry2CTindex
% Convert IEC gantry CS coordinate into pixel index of the CT scan and
% linear index of the CT scan
% If some of the requested voxels are outside of the CT scan, the function returns empty outputs
%
%% Syntax
% |[Aidx , Axyz] = IECgantry2CTindex(A, Beam , Spacing , ImagePositionPatient , sizeCT)|
%
%
%% Description
% |[Aidx , Axyz] = IECgantry2CTindex(A, Beam , Spacing , ImagePositionPatient , sizeCT)| Description
%
%
%% Input arguments
% |A| -_SCALAR MATRIX_- |A(i,:) = [x,y,z,1]| Matrix of 4 vector [x,y,z,1] defining the position of the i-th voxel in IC gantry CS
%
% |Beam| -_STRUCTURES_- Information about the beam
%     * |Beam.GantryAngle| -_SCALAR_- Gantry Angle (deg)
%     * |Beam.PatientSupportAngle| -_SCALAR_- Couch Angle (deg)
%     * |Beam.isocenter| -_SCALAR VECTOR_- [x,y,z] Coordiantes (mm) of the isocentre in the planning CT scan
%
% |Spacing| - _SCALAR VECTOR_ - Pixel size (|mm|) of the displayed images in GUI
%
% |ImagePositionPatient| - _SCALAR VECTOR_ - Coordinate (in |mm|) of the first pixel of the image in the coordinate system of the image
%
% |sizeCT| -_SCALAR VECTOR_- Nulber of voxel in the [x,y,z] dimension for the Ct scan
%
% |warn| -_INTEGER_- What to do if some voxels in |A| are outside the CT scan.
%           * |warn = 0| : return empty |Axyz| and |Aidx|
%           * |warn = 1| : return empty |Axyz| and |Aidx| and issue a warning
%           * |warn = 2| : return empty |Aidx|. |Axyz| contains all indices, even those out of bounds
%
%% Output arguments
%
% |Aidx| -_SCLAR VECTOR_- |Aidx(i)| Linear index of the i-th voxel in the CT scan
%
% |Axyz| -_SCALAR MATRIX_- |Axyz(i,:) = [xi , yi , zi]| The pixel indices along each axis of the CT scan for the i-th pixel
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [Aidx , Axyz] = IECgantry2CTindex(A, Beam , Spacing , ImagePositionPatient , sizeCT , warn)

  if nargin < 6
    warn = 0;
  end

  M = matDICOM2IECgantry(Beam.GantryAngle,Beam.PatientSupportAngle,Beam.isocenter); %Rotate around isocentre
  Adcm = inv(M) * A' ; %Convert coordinates from IEC gantry into DICOM CS. Origin is at isocentre
  Axyz = DICOM2PXLindex(Adcm' , Spacing , ImagePositionPatient); %Shift the origin of the CS from isocentre to |ImagePositionPatient| and then convert to pixel index

  %Check that the requested indices are within the provided CT scan
  if (sum(min(Axyz,[],1) < 1) | sum((sizeCT - max(Axyz,[],1)) < 0) )
    switch warn
      case 0
        %return empty output
        Axyz = [];
        Aidx = [];
      case 1
        %return empty output
        Axyz = [];
        Aidx = [];

        min(Axyz,[],1)
        max(Axyz,[],1)
        sizeCT
        warning('Some of the coordinates are outside of the CT scan')
      case 2
        %return empty Aidx
        Aidx = [];
      end
      
    return
  end

  Aidx = sub2ind(sizeCT,Axyz(:,1),Axyz(:,2),Axyz(:,3));

end
