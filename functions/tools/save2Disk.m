%% save2Disk
% Save a 3D matrix to disk as a DICOM image. as a dose map or a dose rate map
% If a filename of a DICOM plan is provided, the dose map is reference to this plan
%
%% Syntax
% |handles = save2Disk(handles , Imatrix , imgsize , CTinfo , imgName , output_path , planFullPath , doseType)|
%
%
%% Description
% |handles = save2Disk(handles , Imatrix , imgsize , CTinfo , imgName , output_path , planFullPath , doseType)| Description
%
%
%% Input arguments
% |handles| - STRUCTURE_ -  RREGGUI handles
%
% |Imatrix| - SCLAR MATRIX_ - Matrix to be saved as a DICOM file. If ndims(Imatrix)==1, then the vector is reshaped as a 3D matrix (with |imgsize|) and the Y axis is flipped.
%                             If ndims(Imatrix)==3, the matrix is saved without fipping
%
% |imgsize| - _SCALAR VECTOR_ -  [x,y,z] Number of pixel along each dimension of the DICOM CS
%
% |CTinfo| - STRUCTURE_ - DICOM tags to be saved with the image
%
% |imgName| - _STRING_ -  Name of the DICOM file
%
% |output_path| - _STRING_ - Directory in which to save the image
%
% |planFullPath| -_STRING_- Filen name (including full path to the folder) of the DICOM file with the mono-layer plan. Used to link the dose DICOM files
%
% |doseType| -_STRING_- [OPTIONAL. Default = 'RTDOSE'] Define the DICOM tag (3004,0004) Dose Type and (0008,0060) Modality Options are :
%                           * RTDOSE for a dose map : Modality = RTDOSE ; DoseType = RTDOSE ; DoseType = PHYSICAL ;
%                           * PERCENTILE for a percentile dose rate map: Modality = RTDOSERATE ; DoseType = PERCENTILE ;
%                           * DADR for a dose averaged dose rate map: Modality = RTDOSERATE ; DoseType = DADR ;
%
%% Output arguments
%
% |handles| - STRUCTURE_ -  RREGGUI handles
%
%
%% Contributors
% Authors : R. Labarbe, L. Hotoiu (open.reggui@gmail.com)

function [handles, CorrectedName] = save2Disk(handles , Imatrix , imgsize , CTinfo , imgName , output_path , planFullPath , doseType)

  if nargin < 7
    planFullPath = []; % don't use plan info to save the file
  end

  if nargin < 8
    doseType = 'RTDOSE'; % By default, we save a RTDOSE DICOM object
  end

  if (ndims(Imatrix) == 3)
    reshapeFLAG = false;
  else
    reshapeFLAG = true;
  end

  if reshapeFLAG
    fprintf('Reshaping and flipping the Y axis \n')
    Imatrixa = reshape(full(Imatrix),imgsize);
    Imatrixa = flip(Imatrixa,3); %The Zaxis (patient) is inverted when created by the Pij matrices. The Pij matrice create a left handed CS
  else
    Imatrixa = full(Imatrix);
  end

  %Create an info structure to save with the DICOM image
  if (isempty(planFullPath))
      comparison = strcmp('plan_miropt',handles.plans.name); %Get the list of name of all plans in handles
      [~,index] = find(comparison); %Get the index of the plan related to this dose rate map
      info_out = Create_patient_dose_info( handles, CTinfo, handles.plans.info{index});

  elseif ~isempty(planFullPath)
      %Load the info from the mono energy layer plan
      [~,PlanInfo] = load_Plan(planFullPath,'dcm',0); %Load quitely the mono layer plan to get the header info to save in the dose files

      %Create data structure for dose file
      info_out = Create_patient_dose_info( handles, CTinfo, PlanInfo);
  end

  %If this is a dose rate image, change the DICOM tags
  switch doseType
    case 'RTDOSE' % for a dose map : Modality = RTDOSE ; DoseType = RTDOSE ; DoseType = PHYSICAL ;
        info_out.OriginalHeader.DoseType = 'PHYSICAL';
        info_out.OriginalHeader.DoseComment = 'DOSE/PHYSICAL';

    case 'DRAD' %  for a dose rate averaged dose map : Modality = RTDOSE ; DoseType = EFFECTIVE ;
        info_out.OriginalHeader.DoseType = 'EFFECTIVE';
        info_out.OriginalHeader.DoseComment = 'DOSE/DRAD';

    case 'PERCENTILE' % for a percentile dose rate map: Modality = RTDOSERATE ; DoseType = PERCENTILE ;
        info_out.OriginalHeader.Modality = 'RTDOSERATE'; %TODO This crashes in RaysTation and is therefore des-activated. Should be activated again when supported by RS
        info_out.OriginalHeader.DoseType = 'PERCENTILE';
        info_out.OriginalHeader.DoseUnits = 'GY/S'; %This is a dose rate map. Change the default units (Gy) to Gy/s
        info_out.OriginalHeader.DoseComment = 'RTDOSERATE/PERCENTILE';

    case 'MPDR' % for a percentile dose rate map: Modality = RTDOSERATE ; DoseType = PERCENTILE ;
        info_out.OriginalHeader.Modality = 'RTDOSERATE'; %TODO This crashes in RaysTation and is therefore des-activated. Should be activated again when supported by RS
        info_out.OriginalHeader.DoseType = 'MPDR';
        info_out.OriginalHeader.DoseUnits = 'GY/S'; %This is a dose rate map. Change the default units (Gy) to Gy/s
        info_out.OriginalHeader.DoseComment = 'RTDOSERATE/MPDR';

    case 'DADR' % for a dose averaged dose rate map: Modality = RTDOSERATE ; DoseType = DADR ;
        info_out.OriginalHeader.Modality = 'RTDOSERATE'; %TODO This crashes in RaysTation and is therefore des-activated. Should be activated again when supported by RS
        info_out.OriginalHeader.DoseType = 'DADR';
        info_out.OriginalHeader.DoseUnits = 'GY/S'; %This is a dose rate map. Change the default units (Gy) to Gy/s
        info_out.OriginalHeader.DoseComment = 'RTDOSERATE/DADR';
  end

  %Add the info structure to the REGGUI image
  info_out.OriginalHeader.SeriesDescription = strcat(imgName); %strcat removes the trailing white spaces
  [handles,CorrectedName] = Set_reggui_data(handles,imgName,single(Imatrixa),info_out,'images',0);

  % Export Dicom Dose or dose rate
  if (~exist(output_path,'dir'))
    %If the output directory does not exist, create it
    mkdir(output_path)
  end
  fprintf('Exporting image \n')
  handles = Export_image(CorrectedName,fullfile(output_path, imgName) , 'dcm' , handles);

end
