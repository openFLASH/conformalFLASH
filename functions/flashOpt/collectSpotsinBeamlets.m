%% collectSpotsinBeamlets
% Collect the spots (individual Bragg peaks) into sets (=beamlets = SOBP) of spots located at the same [x,y] position.
% The beamlets is a SOBP at the position [x,y]
% The matrix |weight2spot| contains the pointer linking the spots to their beamlets
%
%
%% Syntax
% |[beamletPosition , weight2spot , DRcritical , Plan] =  collectSpotsinBeamlets(Plan , ROI)|
%
%
%% Description
% |[beamletPosition , weight2spot , DRcritical , Plan] =  collectSpotsinBeamlets(Plan , ROI)| Description
%
%
%% Input arguments
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are
% stored. The following data must be present in the structure:
%
% * |Plan.DoseGrid| - _struct_ - Structure containing the information about the dose grid. The following data must be present:
% * |Plan.DoseGrid.nvoxels| - _scalar_ - Total number of voxels of the dose grid, i.e., number of voxels in the CT image.
% * |Plan.Beams| -_VECTOR of STRUCTURES_- Information about the different beams in the plan
%     * |Plan.Beams(i).GantryAngle| -_SCALAR_- Angle (deg) of the i-th beam in the plan
%     * |Plan.Beams(i).PatientSupportAngle| -_SCALAR_- Angle (deg) of the couch for the i-th beam in the plan
%     * |Plan.Beams(i).isocenter| -_SCALAR VECTOR_- [x,y,z] Coordiantes (mm) of the isocentre in the planning CT scan for i-th beam
%
% * |Plan.optFunction| - _struct_ - Structure containing the information about the objective functions set by the user on the dose to the target volume and organs at risk. The following data must be present for each objective function with index i:
% * |Plan.optFunction(i).ROIname| - _string_ - Name of the ROI to which the objective function must be applied.
% * |Plan.optFunction(i).name| - _string_ - Name of the objective function type. The possible types supported so far are: 'min', 'max', 'min_mean' and 'max_mean'.
% * |Plan.optFunction(i).Dref| - _scalar_ -  Reference dose in Gy
% * |Plan.optFunction(i).Vref| - _scalar_ -  Reference volume in 1% only needed for minDVH and maxDVH functions
% * |Plan.optFunction(i).impw| - _scalar_ -  Importance weight for the current objective function (impw > 0)
% * |Plan.optFunction(i).robust| - _scalar_ - Boolean value indicating if the objective function must be considered in all scenarios (if set as robust, equal to 1), or only in the nominal case (if set as non-robust, equal to zero).
%
% * |Plan.NbrRandomScenarios| - _scalar_ - Number of scenarios used to simulate the random setup errors.
% * |Plan.NbrRangeScenarios| - _scalar_ - Number of scenarios used to simulate the range errors.
% * |Plan.Nbr4DScenarios| - _scalar_ - Number of scenarios used to simulate the errors related to breathing motion.
% * |Plan.NbrSystSetUpScenarios| - _scalar_ - Number of scenarios used to simulate the systematic setup errors.
% * |Plan.PlanExistentFile| - _string_ - Path pointing to the location of the file containing the spot doses. If empty, the spot doses (dose influence matrix) will be computed from scratch.
%
% |ROI| - _struct_ - [OTPIONAL] MIROpt structure containing information about all
% volumes in the RTSTRUCT file. The following data must be present in the
% structure:
% * |ROI(i).mask3D.value|- _array_ - |ROI(i).mask3D.value(x,y,z)=1| if the voxel at (x,y,z) is located inside the RT struct
%
%% Output arguments
% |beamletPosition| - CELL VECTOR_ - beamletPosition{b}(i,:) = [x,y] The i-th spot of the b-th beam is located at position [x,y] in the BEV coordinates
%
% |weight2spot| - SCALAR MATRIX_ - |weight2spot(weightIdx,:) = [b,BeamLetNb,l]| The spot at index |weightIdx| is related to the b-th beam and the SOBP # |BeamLetNb| in layer # l
%
% |DRcritical| -_CELL VECTOR_- [If  |ROI| is provided] |DRcritical{b}| gives a vector with the |bemalet| number of all the SOBPs (for beam |b|) that hit a structure for which there is a dose rate constraint
%
% |Plan| -_STRUCTURE_- Updated data about plan
%     * |Plan.bevPTV| - _CELL VECTOR OF SCALAR MATRIX_ - |Plan.bevPTV{b}(x,y)| Value of the paralell projection of the |target| on a a surface perpendicular to the |beam| axis. The surface is centered on |beam| axis. It is a square with size equal to |bev_size|*|bev_size|
%     * |Plan.bevOAR| - _CELL VECTOR OF SCALAR MATRIX_ - |Plan.bevOAR{b}(x,y)| Value of the paralell projection of the union of all OAR on a a surface perpendicular to the |beam| axis. The surface is centered on |beam| axis. It is a square with size equal to |bev_size|*|bev_size|
%     * |Plan.bev_x| -_CELL VECTOR OF SCLAR VECTOR_- |Plan.bev_x{b}(i)| X coordinate of Plan.bevPTV{b}(i,:) or Y coordinate of Plan.bevPTV{b}(:,i) for b-th beam
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [beamletPosition , weight2spot , DRcritical , Plan] =  collectSpotsinBeamlets(Plan, ROI)

flagROI = 0;
if (nargin > 1) && ~isempty(ROI)
  flagROI = 1;
end

beamletPosition = [];
weight2spot = [];
weightIdx = 1;
DRcritical{size(Plan.Beams,2)} = [];

if (~isfield(Plan,'bevPTV') | ~isfield(Plan,'bevOAR'))
  %There is no beam's eye view projections available
  %set the flag to compute them this time
  if flagROI
    fprintf('Projection of BEV \n')
  end
  Project = true;
  Plan.bevOAR = [];
else
  if flagROI
    fprintf('Reuse BEV \n')
  end
  Project = false;
end

  for b = 1: size(Plan.Beams,2) %Loop for each beam
      AllSpots = [];
      ProjStruc = [];
      indx = 1;

      if(Project)
        Plan.bevOAR{b} = [];
      end


      if flagROI
        %If ROI is provided, prepare the size of the BEV
        ExternalROI_ID = getROIByName(ROI, Plan.ExternalROI);
        beam.gantry_angle = Plan.Beams(b).GantryAngle; %- _SCALAR_ - Gantry angle (in degree) of the treatment beam beam
        if (isfield(Plan.Beams(b),'PatientSupportAngle'))
          beam.table_angle = Plan.Beams(b).PatientSupportAngle;
        else
          beam.table_angle =0; % Miropt only support yaw = 0 - _SCALAR_ - Table top yaw angle (degree) of the treatment beam
        end
        beam.isocenter = Plan.Beams(b).isocenter; %_SCALAR VECTOR_ - |beam.isocenter= [x,y,z]| Cooridnate (in mm) of the isocentre in the CT scan for the treatment beam

        %Collect all the ROI from dose rate objectives together to get a BEV with all ROIs
        Ind = find(ROI(Plan.TargetROI_ID).mask3D.value); %Put the PTV in the list of ROI
        for optFidx = 1:length(Plan.optFunction)
          if (Plan.optFunction(optFidx).ID == 8 || Plan.optFunction(optFidx).ID == 9 || Plan.optFunction(optFidx).ID == 10)
            Ind = [Ind; find(ROI(Plan.optFunction(optFidx).ROIindex).mask3D.value)];
          end
        end
        Ind = unique(Ind);
        [x,y,z] = ind2sub(size(ROI(ExternalROI_ID).mask3D.value),Ind);
        bev_size = getBEVsize(x,y,z ,beam);
        bev_size = round(2.5 .* max(bev_size)); %BEV window is equal to the longest diameter of the RT struc + a margin
        x = linspace(-bev_size/2 , bev_size/2 , bev_size); %Make a linear spacing every mm
        SAD = 1000; %mm Arbitrary as the rays are paralell. Place the plane of the source outside of CT scan

        %Project the PTV onto the BEV
        handles.spacing = Plan.CTinfo.Spacing; %- _SCALAR VECTOR_ - Pixel size (|mm|) of the displayed images in GUI
        handles.origin = Plan.CTinfo.ImagePositionPatient; %- _SCALAR VECTOR_ - Coordinate (in |mm|) of the first pixel of the image in the coordinate system of the image

        if (Project)
          [Plan.bevPTV{b} , x] = Project_on_isoplane(handles,ROI(Plan.TargetROI_ID).mask3D.value,beam,SAD,bev_size); %Linear spacing every mm
          Plan.bevPTV{b} = Plan.bevPTV{b} > 0;
          Plan.bevPTV{b} = Plan.bevPTV{b} * 20; %Make PTV bright
          Plan.bev_x{b} = x;
        else
          %Retrieve the x value from previous iteration
          x = Plan.bev_x{b};
        end

        for optFidx = 1:length(Plan.optFunction)
          if (Plan.optFunction(optFidx).ID == 8 || Plan.optFunction(optFidx).ID == 9 || Plan.optFunction(optFidx).ID == 10)

            %Project the RT strcut requiring DR computation onto the BEV
            if (Project)
              %dilatedROI = DilateROI(ROI(Plan.optFunction(optFidx).ROIindex).mask3D.value,Plan.Beams(b).SpotSpacing,Plan.CTinfo.Spacing);
              fprintf('Dilating %s \n' , Plan.optFunction(optFidx).ROIname )
              dilatedROI = DilateROI(ROI(Plan.optFunction(optFidx).ROIindex).mask3D.value , Plan.Beams(b).spotSigma .*2 , Plan.CTinfo.Spacing);
                      %Dilate the ROI by 2 spot sigma. Beyond that limit, the spots will no longer significantly contribute to dose in the ROI and will not impact the percentile dose rate
              fprintf('Projecting %s \n',Plan.optFunction(optFidx).ROIname)
              proj = Project_on_isoplane(handles,dilatedROI,beam,SAD,bev_size);
              proj = (proj > 0); %proj(x,y)= 1 indicates that the BEV pixel is included inside the projection of the structure
              if (isempty(Plan.bevOAR{b}))
                %Project the first OAR on the beam's eye view
                Plan.bevOAR{b} = proj;
              else
                %Add the projection of a new OAR on the BEV
                Plan.bevOAR{b} =  Plan.bevOAR{b} + proj;
                Plan.bevOAR{b} = Plan.bevOAR{b} > 0;
              end
            end

            if Plan.showGraph
                if ishghandle(b)
                  close (b)
                end
                figure(b)
                hold off

                disStr = Plan.bevPTV{b} + 50 .* Plan.bevOAR{b};
                image(x,x,(disStr).*255 ./ max(disStr,[],'all'))
                xlabel('X (mm)')
                ylabel('Y (mm)')
                title(['Beam ',num2str(b)])
                grid on
            end
         end %if (Plan.optFunction(optFidx).ID == 8)
       end %for optFidx
     end %if(flagROI)

      for l = 1:size(Plan.Beams(b).Layers,2) %Loop for each layer
          for idxSpt = 1:numel(Plan.Beams(b).Layers(l).SpotWeights) %Loop for each spot in the layer
            BeamLetNb = 0; %reset Beamlet number, to be sure that a proper number is found
             if (isempty(AllSpots))
               %There is no spot in the list. Add this first spot
               %There is no spot in any previous layer at this position. Add spot to the list
               AllSpots(end+1,:)= Plan.Beams(b).Layers(l).SpotPositions(idxSpt ,:);
               BeamLetNb = size(AllSpots,1);
               if flagROI && ~isempty(Plan.bevOAR{b})
                 %There is a projected structure to check
                 %Check whether this spot is inside the projection of a structure with DR constraints
                 if Plan.showGraph
                   hold on
                   plot(AllSpots(BeamLetNb,1),AllSpots(BeamLetNb,2), 'or')
                 end
                 DRcritical = checkSpotInStructure(b , x , AllSpots(BeamLetNb,:) , BeamLetNb , Plan.bevOAR{b} , DRcritical , Plan.showGraph);
               end

             elseif (isempty(find(getDistances2(AllSpots , Plan.Beams(b).Layers(l).SpotPositions(idxSpt ,:)) == 0)))
               %There is no spot in any previous layer at this position. Add spot to the list
               AllSpots(end+1,:)= Plan.Beams(b).Layers(l).SpotPositions(idxSpt ,:);
               BeamLetNb = size(AllSpots,1);

               %Check whether this spot is inside the projection of a structure with DR constraints
               if flagROI
                 %There is a projected structure to check
                 if Plan.showGraph
                   figure(b)
                   hold on
                   plot(AllSpots(BeamLetNb,1),AllSpots(BeamLetNb,2), 'or')
                 end
                 if ~isempty(Plan.bevOAR{b})
                   DRcritical = checkSpotInStructure(b,x,AllSpots(BeamLetNb,:),BeamLetNb,Plan.bevOAR{b},DRcritical , Plan.showGraph);
                 else
                   DRcritical{b}= [];
                 end
               end %if flagROI

             else
               %There is already a spot at that position in another layer
               BeamLetNb = find(getDistances2(AllSpots , Plan.Beams(b).Layers(l).SpotPositions(idxSpt ,:)) == 0);
             end
             if(BeamLetNb==0)
               error('Did not find matching beamlet number')
             end
             weight2spot(weightIdx,:) = [b , BeamLetNb , l];%The weight w(weightIdx) is included in beam b and is in the beamlet number |BeamLetNb|. It belongs to layer l. It is one of the Bragg peaks composing beamlet # |BeamLetNb|
             weightIdx = weightIdx + 1;
            % The spots weights (w) are ordered in the following manner:
            % [beam1-layer1-spot1, B1L1S2, ... B1L1Sn, B1L2S1,...B1LmSn, B2L1S1,....,BwLmSn]
            % All layers of the same beam have the same (x,y) spot position on the lattice
            % In FLASH, we use Conformal Energy Filter. So we only need to consider ONE trajectory per beam: all layers are created by the Conformal Energy Filter
            % For all layers in the beam, take all the spots covering the largest surface
          end
      end
      beamletPosition{b} = AllSpots; %beamlet{b}(i,:)= [x,y] coordinates of the i-th beamlet to be delivred in the beam b with Conformal Energy Filter
      if flagROI
        DRcritical{b} = unique(DRcritical{b}); %Keep one single copy of all beamlet number in structures
      end %if(flagROI)
      drawnow
  end

end




%-------------------------------
% Check whether the spot is inside the structure
% If yes, add it to DRcritical
%-------------------------------
function DRcritical = checkSpotInStructure(b,x,Spot,BeamLetNb,ProjStruc,DRcritical, showGraph)

    [X,Y] = meshgrid(x,x);

    %Convert the physical dimension into pixel indices
    D = (X-Spot(1)).^2 + (Y - Spot(2)).^2;
    v = min(D,[],'all');
    [xind , yind] = find(D==v);

    if (ProjStruc(xind , yind))
      %The axis of this beamlet hits the projection of the structure
      %record the bemlet index
      if (~isempty(DRcritical{b}))
        DRcritical{b}=[DRcritical{b}(:)' , BeamLetNb];
      else
        DRcritical{b}= BeamLetNb;
      end %if (~isempty

      if showGraph
        figure(b)
        hold on
        plot(Spot(1),Spot(2), '+r')
      end
  end %if (x > 0

end
