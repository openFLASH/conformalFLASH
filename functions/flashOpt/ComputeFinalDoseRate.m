%% ComputeFinalDoseRate
% Compute the dose rate in a region of interest when the spot weight and the spot trajectory is known
%
% The dose rate requires the spot duration and the sweep time between spots. The spot duration is computed using the spot |weightIN|.
% The sweep time between spot depends on the sequence of spots (depending on the wpot weight, some spot may be removed from the sequence)
% The spot duration and the sweep time can be computed using one of the following options:
%  * OPTION 'scanAlgo': Using scanAlgo gateway
%  * OPTION 'Record': Using the provided |Plan.SpotTrajectoryInfo.beam{b}.dT(s)| and |Plan.SpotTrajectoryInfo.beam{b}.TimePerSpot(s)|
%  * OPTION 'Model': Using a machine model defined in |Plan.SpotTrajectoryInfo.beam{b}.Nmaps|
%
%% Syntax
% |[handles] = ComputeFinalDoseRate(Plan , w , handles , ROI)|
%
%
%% Description
% |[handles] = ComputeFinalDoseRate(Plan , w , handles , ROI)| Description
%
%
%% Input arguments
%
% |Plan| - _STRUCTURE_ - The optimised treatment plan
% * |Plan.scenario4D|
% * |Plan.rr_nominal| - _scalar_ - Index for the nominal case with respect to random setup errors.
% * |Plan.rs_nominal| - _scalar_ - Index for the nominal case with respect to range errors.
% * |Plan.Scenario4D(i_4D).RandomScenario(rr).RangeScenario(rs).P| -_SCALAR VECTOR_- Dose influence matrix: |P(vox,spot)| The dose contribution to voxel |vox| of the spot number |spot| for the i_4D-th breathing scenario, rr-th random setup sceanrio and rs-th range error scenario
% *|Plan.SpotTrajectoryInfo| -_STRUCTURE_- Spot trajectory information. |SpotTrajectoryInfo{b}| Parameter for beam |b|.
%    * |SpotTrajectoryInfo.TimingMode| -_STRING_- Option used to evaluate the spot duration. The options are:
%                     * OPTION 'scanAlgo' : Using scanAlgo gateway
%                     * OPTION 'Record' : Using the provided |Plan.SpotTrajectoryInfo.beam{b}.dT(s)| and |Plan.SpotTrajectoryInfo.beam{b}.TimePerSpot(s)|
%                     * OPTION 'Model'  : Using a machine model defined in |Plan.SpotTrajectoryInfo.beam{b}.Nmaps|
%    * |SpotTrajectoryInfo.beam{b}.sobpSequence| -_SCALAR VECTOR_-  Order of the indices of |spot| to sort the spots. |OrderedSOBP = spot(sobpSequence,:)|
%    * |SpotTrajectoryInfo.sobpPosition{b}| - CELL VECTOR_ - beamletPosition{b}(i,:) = [x,y] The i-th spot of the b-th beam is located at position [x,y] in the BEV coordinates
%    * |SpotTrajectoryInfo.weight2spot| - SCALAR MATRIX_ - |weight2spot(weightIdx,:) = [b,BeamLetNb,l]| The spot at index |weightIdx| is related to the b-th beam and the SOBP # |BeamLetNb| in layer # l
%    * |SpotTrajectoryInfo.scanAlgoGW| -_STRUCTURE_- [OPTIONAL. Only required for |TimingMode='scanAlgo'|] ] Information about the scanAlgo gateway
%           * |scanAlgoGW.scanalgoGateway_IP| -_STRING_- IP address, including port, to the scanAlkgo gatewat
%           * |scanAlgoGW.room_id| -_STRING_- Room ID as defined  inthe gateway
%           * |scanAlgoGW.snout_id|  -_STRING_- snout ID as defined in the gateway
%           * |scanAlgoGW.spot_id| -_STRING_- Spot tune ID as defined in the gateway
%    * |SpotTrajectoryInfo.beam{b}.dT| -_SCALAR VECTOR_- [OPTIONAL. Only required for |TimingMode='Record'|] |dT(s)| Sweep (ms) to move from the s-1-th spot to the s-th spot in |sobpSequence(s)|
%    * |SpotTrajectoryInfo.beam{b}.TimePerSpot| -_SCALAR VECTOR_- [OPTIONAL. Only required for |TimingMode='Record'|] |TimePerSpot(s)| Duration (ms) of the s-th spot in |sobpSequence(s)|
%    * |SpotTrajectoryInfo.beam{b}.Nmaps| -_STRUCTURE_-  [OPTIONAL. Only required for |TimingMode='Model'|] Topological information on the initial spot sequence
%         * |Nmaps.NeighbourghMap| -_SCALAR VECTOR_- |NeighbourghMap(d,i)| d=# of delivered spot; i=# of impacted spot; |NeighbourghMap(d,i)|  = 0 if there is an impact
%         * |Nmaps.NeighbourghWeightMap| -_SCALAR VECTOR_- |NeighbourghWeightMap(d,i)| fraction of the dose of spot d (# of delivered spot) that is also delivered at spot i (i=# of impacted spot)
%         * |Nmaps.NeighbourghTimeMap| -_SCALAR MATRIX_- |NeighbourghTimeMap(d,i)| Time (ms) required to go from spot d to spot i
%
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are stored. The following data must be present in the structure:
%  * |Plan.fractions| -_SCALAR_- Number of fraction to deliver the dose. Used to compute the dose per fraction (and hence dose rate per fraction)
%  * |Plan.Inozzle| -_SCALAR_- Nozzle current (nA) during spot delivery
%  * |Plan.Beams(b).spotSigma| -_SCALAR_- Spot lateral sigma (mm) at maximum of deepest Bragg peak along the optical axis
%  * |Plan.DoseGrid.size| -_SCALAR VECTOR_- [Nx, Ny, Nz] Number of pixels of the Dose map along each dimension
%  * |Plan.output_path|-_STRING_- Folder in which the results are saved
%  * |Plan.BDL| -_STRUCTURE_- Spot information from the Beam Data Library. See |getSpotFromBDL| for more information
%  * |Plan.SaveDoseBeamlets| -_BOOL_- save the dose of each beamlet in the reference frame of the CT with aperture and with the resolution of that CT: dcm (DICOM format) , sparse (sparse matrix) , false (not saved)
%
%
% |handles| - _STRUCT_ - REGGUI data structure.
%
% |ROI| - _struct_ - MIROpt structure containing information about all volumes in the RTSTRUCT file. The following data must be present in the structure
%     If empty, then the dose rate is computed in the whole dose map
% * |ROI(i).mask1D| - _array_ - Logical column vector storing a binary mask for ROI i (voxels inside the volume of interest are equal to 1, and those outside are equal to zero).
%
%% Output arguments
%
% |handles| - _STRUCT_ - REGGUI data structure.
%
% |doseRatesCreated| -_CELL VECTOR of STRINGS_- List of names of dose rate maps saved to disk
%
%
%% Contributors
% Authors : R. Labarbe, Lucian Hotoiu (open.reggui@gmail.com)

function [handles, doseRatesCreated] = ComputeFinalDoseRate(Plan, handles, ROI, collectSpotsInBeamlets)

    if nargin < 4
        collectSpotsInBeamlets = true;
    end

    if collectSpotsInBeamlets
        [~ , ~, ~, Plan] = collectSpotsinBeamlets(Plan, ROI); %display the BEV images
    end

    Pij = Plan.Scenario4D(1).RandomScenario(Plan.rr_nominal).RangeScenario(Plan.rs_nominal).P ; %The dose influence matrix for the nominal case of the first breathing phase
    if Plan.showGraph
      plots_BEV = 1:numel(Plan.Beams); %Figure # where to plot the spot trajectory
    else
      plots_BEV = []; %User has turned off display
    end


    %Construct the vector of SOBP weights from the content of the plan
    w = plan2weight(Plan); %This is the spot weight PER FRACTION

    %In SpotWeightsOptimization at line 85, the weight are divided by the number of fractions.
    %But getDRa expects the total dose, for whole plan
    %So we must rescale the weight for whole treatment
    w = w .* Plan.fractions;

    %Compute the total dose delivered to a voxel.
    %This will be used for the dose threshold in the dose rate computation
    D = sparse(Pij * w');

    %For info: display dose stat in PTV
    if ~isempty(ROI)
      targetROIid = getROIByName(ROI,  Plan.TargetROI);
      Dose = D(ROI(targetROIid).mask1D);
      fprintf('Dose in %s : %f <= D(Gy)= %f <= %f \n',Plan.TargetROI,full(min(sum(Dose,2))),full(mean(sum(Dose,2))),full(max(sum(Dose,2))));
    end

    for optFidx = 1:length(Plan.optFunction)
        %Loop for each constraints
        if (Plan.optFunction(optFidx).ID == 8 || Plan.optFunction(optFidx).ID == 9)
            %This is a constraint on dose rate
            if ~isfield(Plan.optFunction(optFidx),'Vref')
              %Use the default percentile
              percentile = 0.01; %Default percentile for dose rate
            else
              percentile = (1 - Plan.optFunction(optFidx).Vref) ./2;
            end

            plotID.plots_BEV = plots_BEV;
            plotID.plot_DR = 79 + optFidx;

            if ~isfield(Plan, 'version')
              Plan.version = 'default';
            end
            switch Plan.version
                case 'Upenn'
                  plotID.pDR_D = 499 + optFidx ;
                  plotID.SaveHisto = true; %Save the dose rate histograms as YAML and JSON
                otherwise
                  plotID.pDR_D = [];
                  plotID.SaveHisto = false;
            end

            if isempty(ROI)
              %No mask provided. Compute dose rate in full volume
              ROImask = [];
            else
              %Mask for the structre were provided. Compute dose rate only in defined mask
              ROImask = ROI(Plan.optFunction(optFidx).ROIindex).mask1D;
            end

            if isfield(Plan.optFunction(optFidx), 'DMF') & isfield(Plan.optFunction(optFidx) , 'DR50') & ~isempty(plotID.pDR_D)
              [~, DRaStru, ~, DADR, ~, DRAD, ~, ~ , MPDR]  = getDRa(Plan.SpotTrajectoryInfo, w, Pij, Plan, Plan.optFunction(optFidx).Dref,  ROImask, D, Plan.optFunction(optFidx).DMF, Plan.optFunction(optFidx).DR50, percentile, Plan.optFunction(optFidx).ROIname, plotID ) ;
            else
              DRAD = [];
              [~, DRaStru, ~, DADR, ~, ~, ~, ~ , MPDR]  = getDRa(Plan.SpotTrajectoryInfo, w, Pij, Plan, Plan.optFunction(optFidx).Dref,  ROImask, D, [], [] , percentile, Plan.optFunction(optFidx).ROIname, plotID ) ;
            end


            doseRatesCreated = {};
            for b = 1:length(DRaStru)
              %Loop for each beam
              path2beamResults = getOutputDir(Plan.output_path , b);
              planFullPath = fullfile(path2beamResults,'Plan');

              %Save the percentile dose rate
              [handles, doseRateName] = save2Disk(handles, DRaStru{b}, Plan.DoseGrid.size, Plan.CTinfo, ['DRprct_beam_',num2str(b),'_in_' , Plan.optFunction(optFidx).ROIname], path2beamResults , planFullPath ,'PERCENTILE');
              doseRatesCreated{1} = doseRateName;

              %Save the maximum percentile dose rate
              [handles, doseRateName] = save2Disk(handles, MPDR{b}, Plan.DoseGrid.size, Plan.CTinfo, ['maxDRprct_beam_',num2str(b),'_in_' , Plan.optFunction(optFidx).ROIname], path2beamResults , planFullPath ,'MPDR');
              doseRatesCreated{4} = doseRateName;

              %Save the dose averaged dose rate
              [handles, doseRateName] = save2Disk(handles, DADR{b}, Plan.DoseGrid.size, Plan.CTinfo, ['DADR_beam_',num2str(b),'_in_' , Plan.optFunction(optFidx).ROIname], path2beamResults , planFullPath ,'DADR');
              doseRatesCreated{2} = doseRateName;

              %Save the dose-rate averaged dose
              if ~isempty(DRAD)
                [handles, doseRateName] = save2Disk(handles, DRAD{b}, Plan.DoseGrid.size, Plan.CTinfo, ['DRAD_beam_',num2str(b),'_in_' , Plan.optFunction(optFidx).ROIname], path2beamResults , planFullPath , 'DRAD');
                doseRatesCreated{3} = doseRateName;
              end

              %NB: the file name must be less than 64 characters, otherwise (0008,103E)	SeriesDescription will rise a warning at export

            end %for b
        end %if(Plan.optFunction(optFidx)
    end %for optFidx
end
