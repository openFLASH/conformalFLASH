%% getRS_doses
% Build the dose influence matrix
% using the dose per beamlets (in the CS of the original CT scan) computed with computeDoseWithCEF
%
% If |Plan.SpotTrajectoryInfo| is not defined, the field is added to the Plan structure in order to
% define a spot trajectory that matches the order fo the spots defined  in the plan.
%
%% Syntax
% |[handles, Plan] = getRS_doses(Plan, handles)|
%
%
%% Description
% |[handles, Plan] = getRS_doses(Plan, handles)| Description
%
%
%% Input arguments
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are stored.
%               This can be a mono or multi energy layer plan
%
% |handles| - _STRUCT_ - REGGUI data structure.
%
%% Output arguments
%
% |PlanMono| - _struct_ - MIROpt structure where all the plan parameters are stored.
%                   This is a mono energy layer plan
%
% |handles| - _STRUCT_ - REGGUI data structure.
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [handles, PlanMono ] = getRS_doses(Plan, handles)

  if numel(Plan.Beams) ~= 1
    numel(Plan.Beams)
    error('There are more than one beam in the plan. FLASH plans should have one single beam')
  end

  CTimageName = Plan.CTname;

  %If this is a multilayer plan, convert into a single layer plan
  if (numel(Plan.Beams.Layers) ~= 1)
      %Create amonolayer plan
      fprintf('Computing a monolayer plan \n')
      PlanMono = CreatePlanMonoLayer(Plan , 'Plan' , Plan.protonsFullDose); %Create the monolayer plan so we have the spot list and their position
  else
      %this is already a monolayer plan
      fprintf('Plan is monolayer \n')
      PlanMono = Plan;
  end

  %Load the dose map of ech beamlet and store it in a sparse matrix
  %Keep only the voxels of the selected structures to save memory
  fprintf('Loading the influence matrix \n')
  OptConfig.BeamletsMatrixPrecision = 'd';
  PlanMono.PlanExistentFile = PlanMono.output_path;
  NbBeamlets = numel(PlanMono.Beams.Layers(1).SpotWeights); %Number of PBS spots in plan with energy monolayer
  PlanMono.Scenario4D(1).RandomScenario(PlanMono.rr_nominal).RangeScenario(PlanMono.rs_nominal).P = sparse(PlanMono.DoseGrid.nvoxels , NbBeamlets); %|Pij(vox,spot)|. Create  zeros sparse matrix
  w = plan2weight(PlanMono); %The beamlet weight in the monolayer plan. This weight is the sum over all fractions

  %Load the dose map computed by MCsquare by computeDoseWithCEF and
  %build the dose influence matrix
  b=1; %FLASH plan have one single beam
  nvoxels = prod(handles.size);
  PlanMono.Scenario4D(1).RandomScenario(PlanMono.rr_nominal).RangeScenario(PlanMono.rs_nominal).P = sparse(nvoxels , NbBeamlets);
          %Create empty sparse matrix with proper size
          %This will make reading faster because it avoids resizing of matrices

  for spt = 1:NbBeamlets

      SptPos = PlanMono.Beams.Layers(1).SpotPositions(spt,:);
      fprintf('Loading dose for PBS spots %d of %d at [%f mm, %f mm] \n',spt,NbBeamlets,SptPos(1),SptPos(2))
      DoseFileName = ['dose_UserGrid_spt' , num2str(spt) , '_spt_X_',num2str(round(SptPos(1))),'_Y_',num2str(round(SptPos(2)))]; %%Define the name of the output dose file
      DosePath = fullfile(PlanMono.output_path,'Outputs','Outputs_beam1','CEF_beam','Outputs', 'SpotDoseInCT'); %Path to the file with the same orientation as the original CT scan

      switch Plan.SaveDoseBeamlets
        case 'dcm'
            handles = Import_image(DosePath, DoseFileName ,'dcm','MCsquare_Dose',handles); %Load the dose map
            Dose = Get_reggui_data(handles,'MCsquare_Dose');
                          %The dose map is oriented on the same grid as the original CT scan
                          %The Matlab importer has rescaled the dose map on the same grid as the CT scan
                          %If the MCsquare dose map was tallied on a larger grid than the high res CT, then the dose map saved on disk  has the spatial resolution of |CEFDoseGrid|
                          %If handles.auto_mode=true, and handles.spatialpropsettled = true then Import_image.m will resample the low resolution dose map into the spatial resolution of the hi-res CT scan


        case 'sparse'
            data = load(fullfile(DosePath , [DoseFileName '.mat']));
            Dose = data.doseBeamlet;
        end

      Dose = Dose ./ (w(spt) .* PlanMono.fractions); %Normalise the dose for w=1 and for 1 fraction
                  %In SpotWeightsOptimization at line 85, the weight are divided by the number of fractions.
      temp = flip(Dose,3);
      dose1D = temp(:); %spare matrix with the dose influence of the spot spt
      dose1D = sparse(double(full(PlanMono.OptROIVoxels_nominal) .* dose1D)); %Apply the mask to force to zero the voxels outside of the RT struct of interrest. This saves memory
              %the .* product should use full matrices. Product .* with sparse matrices seems buggy in Matlab
      PlanMono.Scenario4D(1).RandomScenario(PlanMono.rr_nominal).RangeScenario(PlanMono.rs_nominal).P(:,spt) =  dose1D; %|Pij(vox,spot)|

      switch Plan.SaveDoseBeamlets
        case 'dcm'
          handles = Remove_image('MCsquare_Dose', handles); %Remove the spot dose map from local handles. We will load a new dose map with the smae name for the next spot
      end

  end % for spt

end
