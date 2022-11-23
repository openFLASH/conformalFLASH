%% getIDD
% Build integrated depth dose profile
%
%% Syntax
% |[IDD , Zg] = getIDD(Xg , Yg, GantryAngle , table_angle, Isocentre , pxlSize , DoseSpot)|
%
%
%% Description
% |[IDD , Zg] = getIDD(Xg , Yg, GantryAngle , table_angle, Isocentre , pxlSize , DoseSpot)| Description
%
%
%% Input arguments
% |Xg| -_SCALAR MATRIX_- |Xg(i,j)| X coordinate in IEC Gantry at which DoseSpot(i,j) is evaluated
%
% |Yg| -_SCALAR MATRIX_- |Yg(i,j)| X coordinate in IEC Gantry at which DoseSpot(i,j) is evaluated
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
% |IDD| - _SCALAR VECTOR_ - |IDD(j)| Integrated depth dose at depth  |Zg(i)| in IEC gantry CS
%
% |Zg| - _SCALAR VECTOR_ - |Zg(j)| Depth (mm) in IEC gantry at which the IDD is evaluated
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)


function [IDD , Zg] = getIDD(Xg , Yg, GantryAngle , table_angle, ImagePositionPatient , isocenter , pxlSize , DoseSpot , showGraph)

  NbPts = numel(Xg);
  Zg = -320:1:320; %mm along the IEC Z axis.
  for idxZg = 1:numel(Zg)
    PtsG = [Xg(:) , Yg(:) , Zg(idxZg) .* ones(NbPts,1) ];
    DoseMap = getSlice(PtsG , GantryAngle , table_angle , ImagePositionPatient , isocenter , pxlSize , DoseSpot);
    if ~isempty(DoseMap)
      IDD(idxZg) = sum(DoseMap,'all');
    else
      IDD(idxZg) = 0; %We are outside of the CT scan. Arbitrarily set to 0
    end
  end

  if showGraph
    figure(210)
    plot(Zg,IDD)
    xlabel('Zg (mm)')
    ylabel('IDD')
    grid on
  end

end
