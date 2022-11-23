%% setApertureInCT
% Add an aperture in the CT scan for all beams with tag "ApertureBlock"
%
%% Syntax
% |[Plan , handles ] = setApertureInCT(handles , Plan , CTname , PTV , Body)|
%
%
%% Description
% |[Plan , handles ] = setApertureInCT(handles , Plan , CTname , PTV , Body)| Description
%
%
%% Input arguments
% |handles| -_STRUCTURE_- REggui data handle. The CT scan is stored in |handles| in the image with name |Plan.CTname|.
%
% |Plan| -_STRUCTURE_- Information about the treament plan
%  * |Plan.DoseGrid| - _struct_ - Structure containing the information about the dose grid. The following data must be present:
%  * |Plan.DoseGrid.nvoxels| - _scalar_ - Total number of voxels of the dose grid, i.e., number of voxels in the CT image.
%  * |Plan.CTinfo| -_STRUCTURE_- DICOM header of the CT scan
% * |Plan.ScannerDirectory| - _STRING_ - Name of the folder containing the definition of the CT scanner properties in MCsquare in folder "plugins\openMCsquare\lib\Scanners"
% * |Plan.Beams| -_VECTOR of STRUCTURES_- Information about the different beams in the plan
%     * |Plan.Beams(b).GantryAngle| -_SCALAR_- Angle (deg) of the b-th beam in the plan
%     * |Plan.Beams(b).isocenter| -_SCALAR VECTOR_- [x,y,z] Coordiantes (mm) of the isocentre in the planning CT scan for b-th beam
%     * |Plan.Beams(b).PatientSupportAngle| -_SCALAR_- Angle (deg) of the couch for the b-th beam in the plan
%     * |Plan.Beams(b).Iso2Skin| -_SCALAR_- Isocentre (mm) To Skin Distance along the proton beam axis
%     * |Plan.Beams(b).Layers(l).Energy| -_SCALAR_- Energy (MeV) of the l-th layer of the b-th beam
%
% |CTname| -_STRING_- Name of the CT image in handles.images into which the device is to be inserted
%
%
%% Output arguments
%
% |Plan| -_STRUCTURE_- Updated treatment plan.
%
% |handles| -_STRUCTURE_- Updated REggui data handle. The CT scan is stored in |handles| in the image with name |Plan.CTname|.
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [Plan , handles ] = setApertureInCT(handles , Plan , CTname)

  CT = Get_reggui_data(handles,CTname,'images'); %Update the CT scan with the aperture block in handles
  apertureFlag = false; %Is there any aperture in this plan ?
  CTsize = size(CT);
  ImagePositionPatient = handles.origin;

  for b = 1: size(Plan.Beams,2) %Loop for each beam

    if Plan.Beams(b).ApertureBlock

        Abrass  = apertureContour2Block(Plan.Beams(b) ,  handles.spacing  , Plan.BDL);

        %Insert the aperture block inside the CT scan
        %-------------------------------------------
        % Get the Hu value of the material in the table HU_Material_Conversion.txt
        % In the table, if the boundary value is HU, then the material start at HU + 1
        % For example: 2000  2     # carbon
        % means that the carbon (ID =2) starts at HU = 2001 inclu
        % Therefore, when setting the value of the voxel, add 1 to the HU from the calibrationtable so that the correct mateiral ID is identified.
        % The density of the material is lineraly interpolated. This is why we added +1 to the values in HU_Density_Conversion.txt
        HUbrass = getMaterialSPR('Brass' , Plan.ScannerDirectory) + 1 ; %Hounsfield unit associated to brass in the material file
        HUair =  getMaterialSPR('Schneider_Air' , Plan.ScannerDirectory) + 1; %Hounsfield unit associated to air in the material file

        [CT , ImagePositionPatient] = insertDeviceInCT(CT , Abrass, HUbrass , Plan.Beams(b) , handles.spacing , ImagePositionPatient , HUair);
        CTsize = size(CT); %The CT size was increased. Update the parmaeter.
        apertureFlag = true;

  end %for b

  %Update the handles and plan
  [handles , Plan] = updateAllImages(handles , Plan , CT , ImagePositionPatient , HUair  , CTname);
  handles = Set_reggui_data(handles,'CTwithAperture',CT,Plan.CTinfo,'images',1); %Create a new image with only the aperture


end
