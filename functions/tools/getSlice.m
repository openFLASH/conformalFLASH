%% getSlice
% Interpolate a 2D slice in the 3D dose map (in DICOM CS) from a Beam's eye view at a specified gantry angle and couch angle
% If some of the requested points are outside of the CT scan, the function returns an empty dose map
%
%% Syntax
% |DoseMap = getSlice(PtsG , GantryAngle , table_angle , Isocentre , pxlSize , DoseSpot)|
%
%
%% Description
% |DoseMap = getSlice(PtsG , GantryAngle , table_angle , Isocentre , pxlSize , DoseSpot)| Description
%
%
%% Input arguments
% |PtsG| -_SCALAR MATRIX_- |PtsG(i,:)= [x,y,z]| Coordinate (mm) of the pixels belonging to the 2D slice, expressed in the gantry CS. Z is a constant
%
% |GantryAngle| -_SCALAR_- Gantry angle (deg)
%
% |table_angle| - _SCALAR_ -  Yaw angle of the PPS table |degree|
%
% |ImagePositionPatient| - _SCALAR VECTOR_ - Coordinate (in |mm|) of the first pixel of the image in the coordinate system of the image
%
% |isocenter| -_SCALAR VECTOR_- [x,y,z] Coordiantes (mm) of the isocentre in the DICOM CS
%
% |pxlSize| -_SCALAR_- Size (mm) of the square pixels of the 2D slice
%
% |DoseSpot| -_SCALAR MATRIX_- |DoseSpot(x,y,z)| Dose at the pixel coordinate (x,y,z) in the CT scan
%
%
%% Output arguments
%
% |DoseMap| -_SCLAR MATRIX_- |DoseMap(x,y)| Interpolated dose value at the |PtsG| points
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)


function DoseMap = getSlice(PtsG , GantryAngle , table_angle , ImagePositionPatient , isocenter , pxlSize , DoseSpot)

    PtsG = [PtsG , ones(size(PtsG,1),1)]; %Make a 4 vector to work with matrices
    Beam.GantryAngle = GantryAngle;  %| -_SCALAR_- Gantry Angle (deg)
    Beam.PatientSupportAngle = table_angle; %| -_SCALAR_- Couch Angle (deg)
    Beam.isocenter = isocenter; %| -_SCALAR VECTOR_- [x,y,z] Coordiantes (mm) of the isocentre in the planning CT scan
    [~ , PTsDICOM] = IECgantry2CTindex(PtsG, Beam , pxlSize , ImagePositionPatient , size(DoseSpot));
    if ~isempty(PTsDICOM)
      DoseMap = interp3( DoseSpot , PTsDICOM(:,2),PTsDICOM(:,1),PTsDICOM(:,3),'linear',0); %Note that X is the second index and Y is the first index in interp3
    else
      DoseMap = [];
    end

end
