%% enlargeCT
% Create a large CT scan |CT2| and insert a smaller CT scan |CT| into it
%
%% Syntax
% |CT2 = enlargeCT(CT , ImagePositionPatientOLD , ImagePositionPatientNEW , Spacing , sCT2 , HUvoid)|
%
%
%% Description
% |CT2 = enlargeCT(CT , ImagePositionPatientOLD , ImagePositionPatientNEW , Spacing , sCT2 , HUvoid)| Description
%
%
%% Input arguments
% |CT| - _SCALAR MATRIX_ - Small CT scan |CT(x,y,z)=HU| gives the intensity (Hounsfield unit) of the voxel at location (x,y,z)
%
% |ImagePositionPatientOLD| - _SCALAR VECTOR_ - Coordinate (in |mm|) of the first pixel of the image in the coordinate system of the small CT scan |CT|
%
% |ImagePositionPatientNEW| - _SCALAR VECTOR_ - Coordinate (in |mm|) of the first pixel of the image in the coordinate system of the large CT scan |CT2|
%
% |Spacing| - _SCALAR VECTOR_ - Pixel size (|mm|) of the displayed images in GUI
%
% |sCT2| -_SCALAR VECTOR_- [x,y,z] number of pixel along each dimension of the new larger CT scan
%
% |HUvoid| -_SCALAR_- Hounsfield unit to use for unoccupied voxels in larger CT scan
%
%% Output arguments
%
% |CT2| - _SCALAR MATRIX_ - Large CT scan |CT2(x,y,z)=HU| gives the intensity (Hounsfield unit) of the voxel at location (x,y,z)
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function CT2 = enlargeCT(CT , ImagePositionPatientOLD , ImagePositionPatientNEW , Spacing , sCT2 , HUvoid)

    CT2 = ones(sCT2) .* HUvoid;
    sCT0 = size(CT);
    ori0 = DICOM2PXLindex(ImagePositionPatientOLD' , Spacing , ImagePositionPatientNEW , true); %Get the index coordinate of the old CT into the new CT

    CT2(ori0(1):ori0(1)+sCT0(1)-1 , ori0(2):ori0(2)+sCT0(2)-1 , ori0(3):ori0(3)+sCT0(3)-1) = CT; %insert original CT in expanded CT
    CT = CT2; %This is the expanded CT

end
