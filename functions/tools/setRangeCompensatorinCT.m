%% setRangeCompensatorinCT
% Add the range compensator in the CT scan for all beams with tag "RangeCompensator"
% The HU of the range compensator is the same as the HU of the material for the CEM
% The CT scan is the image with name |CTname| stored in |handles|
% Expand the CT if it is too small.
% Save the updated CT in the |handles| with name |CTname|
% Update the |Plan| structure with the update CT size
%
%% Syntax
% |[Plan , handles ] = setRangeCompensatorinCT(handles , Plan , CTname)|
%
%
%% Description
% |[Plan , handles ] = setRangeCompensatorinCT(handles , Plan , CTname)| Description
%
%
%% Input arguments
% |handles| -_STRUCTURE_- REggui data handle. The CT scan is stored in |handles| in the image with name |Plan.CTname|.
%
% |Plan| -_STRUCTURE_- Information about the treament plan
%  * |Plan.DoseGrid| - _struct_ - Structure containing the information about the dose grid. The following data must be present:
%  * |Plan.DoseGrid.nvoxels| - _scalar_ - Total number of voxels of the dose grid, i.e., number of voxels in the CT image.
%  * |Plan.CTinfo| -_STRUCTURE_- DICOM header of the CT scan
% * |Plan.Spike.MaterialID| - _STRING_ - Name of the CEM material, as defined in the file "plugins\openMCsquare\lib\Materials\list.dat"
% * |Plan.ScannerDirectory| - _STRING_ - Name of the folder containing the definition of the CT scanner properties in MCsquare in folder "plugins\openMCsquare\lib\Scanners"
% * |Plan.Beams| -_VECTOR of STRUCTURES_- Information about the different beams in the plan
%     * |Plan.Beams(b).GantryAngle| -_SCALAR_- Gantry Angle (deg)
%     * |Plan.Beams(b).PatientSupportAngle| -_SCALAR_- Couch Angle (deg)
%     * |Plan.Beams(b).isocenter| -_SCALAR VECTOR_- [x,y,z] Coordiantes (mm) of the isocentre in the planning CT scan
%     * |Plan.Plan.Beams(b).RangeCompensator| -_STRUCT_- Structure with description of the range compensator
%        * |Plan.Beams(b).RangeCompensator.IsocentertoCompensatorTrayDistance| -_SCALAR_- Distance (mm) from isocentre to upstream side of the range compensaot
%        * |Plan.Beams(b).RangeCompensator.CompensatorColumns| -_SCALAR_- Number of pixels along the X dimension
%        * |Plan.Beams(b).RangeCompensator.CompensatorRows| -_SCALAR_- Number of pixels along the Y dimension
%        * |Plan.Beams(b).RangeCompensator.CompensatorPixelSpacing| -_SCALAR VECTOR_- [x,y] size (mm) of the pixels projected into the isocentre plane
%        * |Plan.Beams(b).RangeCompensator.CompensatorPosition| -_SCALAR VECTOR_- [x,y] Coordinate (mm) of the voxels in the IEC gantry CS
%        * |Plan.Beams(b).RangeCompensator.CompensatorThicknessData| -_SCALAR MATRIX_- |CompensatorThicknessData(x,y)| Height (mm) of the range compensator at pixel (x,y)
%
% |CTname| -_STRING_- Name of the CT image in handles.images into which the device is to be inserted
%
%% Output arguments
%
% |handles| -_STRUCTURE_- Updated REggui data handle. The CT scan is stored in |handles| in the image with name |Plan.CTname|.
%
% |Plan| -_STRUCTURE_- Updated structure defining a multi energy layer PBS plan
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [Plan , handles ] = setRangeCompensatorinCT(handles , Plan , CTname)

  CT = Get_reggui_data(handles,CTname,'images'); %Update the CT scan with the aperture block in handles
  HUcem = getMaterialPropCT(Plan.Spike.MaterialID , Plan.ScannerDirectory) + 1 ; %Hounsfield unit associated to CEM in the material file
  HUair =  getMaterialPropCT('Schneider_Air' , Plan.ScannerDirectory) + 1; %Hounsfield unit associated to air in the material file
  Zres = min(handles.spacing) ./3;
  ImagePositionPatient = handles.origin;

  for b = 1: size(Plan.Beams,2) %Loop for each beam

      if isfield(Plan.Beams(b) , 'RangeCompensator')

        mag = magnification(Plan.Beams(b).RangeCompensator.IsocentertoCompensatorTrayDistance , Plan.BDL);
        fprintf('Magnification factor X : %f \n', mag(1))
        fprintf('Magnification factor Y : %f \n', mag(2))

          %Define the coordinates of the voxels to ad to CT
          X = (1:Plan.Beams(b).RangeCompensator.CompensatorColumns) .*  Plan.Beams(b).RangeCompensator.CompensatorPixelSpacing(1) + Plan.Beams(b).RangeCompensator.CompensatorPosition(1); %Coordinate (mm) of the voxels in the IEC gantry CS
          Y = (1:Plan.Beams(b).RangeCompensator.CompensatorRows) .*  Plan.Beams(b).RangeCompensator.CompensatorPixelSpacing(2) + Plan.Beams(b).RangeCompensator.CompensatorPosition(2);

          minThick = min(Plan.Beams(b).RangeCompensator.CompensatorThicknessData(Plan.Beams(b).RangeCompensator.CompensatorThicknessData > 0),[],'all');
          maxThick = max(Plan.Beams(b).RangeCompensator.CompensatorThicknessData,[],'all');

          Acem = []; %Indices of the voxels containing brass

          for Z = minThick : Zres : maxThick
            [Xi , Yi] = find((Plan.Beams(b).RangeCompensator.CompensatorThicknessData >= Z) & (Plan.Beams(b).RangeCompensator.CompensatorThicknessData > 0)); %Coordinate of the voxels containing range compensator
            Xcem = X(Xi) .* mag(1); %Project the CEM from the isocentre plane into the CEM plane
            Ycem = Y(Yi) .* mag(2);
            Zp = Plan.Beams(b).RangeCompensator.IsocentertoCompensatorTrayDistance - Z ;
            Ap = [Xcem',Ycem']; %X,Y coordinate of voxels containing brass
            Acem = [Acem ; Ap , ones(size(Ap,1),1).* Zp , ones(size(Ap,1),1)];
          end

          %Insert the range compensator in the CT
          [CT , ImagePositionPatient] = insertDeviceInCT(CT , Acem, HUcem , Plan.Beams(b) , handles.spacing , ImagePositionPatient, HUair);
      end

  end %for b


  %Update the handles and plan
  [handles , Plan] = updateAllImages(handles , Plan , CT , ImagePositionPatient , HUair  , CTname);


end
