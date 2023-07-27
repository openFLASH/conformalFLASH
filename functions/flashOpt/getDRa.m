%% getDRa
% Compute different dose rate indicators for the specified structures:
%  * average (_a) and median (_m) of dose rate in the OAR for all the voxels in the organ that receives a dose larger than a specified threshold. The follwoind dose rate are used:
%       * percentile dose rate |DRa| |DRm| of each voxel inside the organ defined by |ROImask| for a specified spot sequence. This is the dose rate at which 98% of the dose is delivered by the beam. The dose rate is compute using the dose per fraction and per beam.
%       * dose averaged dose rate |DADR| [3]
%
% The dose rate indicators are computed per beam.
% Note that the spot weights provided to the function are spot weight for the TOTAL dose (sum of all fractions)
%
% The dose rate requires the spot duration and the sweep time between spots. The spot duration is computed using the spot |weightIN|.
% The sweep time between spot depends on the sequence of spots (depending on the wpot weight, some spot may be removed from the sequence)
% The spot duration and the sweep time can be computed using one of the following options:
%  * OPTION 'scanAlgo': Using scanAlgo gateway
%  * OPTION 'Record': Using the provided |Plan.SpotTrajectoryInfo.beam{b}.dT(s)| and |Plan.SpotTrajectoryInfo.beam{b}.TimePerSpot(s)|
%  * OPTION 'Model': Using a machine model defined in |Plan.SpotTrajectoryInfo.beam{b}.Nmaps|
%
%
%% Syntax
% |[dr, DRa, SpotTrajectoryInfo, drm, DADR, DADRm, DRAD, DRADm, SpotTiming] = getDRa(SpotTrajectoryInfo , weightIN , Pij , Plan , Dref , ROImask , DoseAtPxl , DMF, DR50, percentile , ROIName , plotID)|
%
%
%% Description
%
% |[dr, DRa, SpotTrajectoryInfo, drm, DADR, DADRm, DRAD, DRADm, SpotTiming] = getDRa(SpotTrajectoryInfo , weightIN , Pij , Plan , Dref , ROImask , DoseAtPxl , DMF, DR50, percentile , ROIName , plotID)| Use the provided (mono-energy) spot trajectory and compute the dose rate
%
%% Input arguments
% *|SpotTrajectoryInfo| -_STRUCTURE_- Spot trajectory information. |SpotTrajectoryInfo{b}| Parameter for beam |b|.
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
% |weightIN| -_SCLAR VECTOR_- |weightIN(j)| the weight of the j-th spot. D(i) = P*weightIN gives the TOTAL dose at the i-th voxel.
%
% |Pij| -_SCALAR MATRIX_- dose influence matrix: |Pij(vox,spot)| The dose contribution to voxel |vox| of the spot number |spot|
%
% |ROImask| -_SCLAR VECTOR_- [If empty, compute dose rate for all voxels] MAsk identifying the position of the voxel contained in the organ. |ROImask(i)=1| if the i-th voxel is inside the organ
%
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are
% stored. The following data must be present in the structure:
%  * |Plan.fractions| -_SCALAR_- Number of fraction to deliver the dose. Used to compute the dose per fraction (and hence dose rate per fraction)
%  * |Plan.Inozzle| -_SCALAR_- Nozzle current (nA) during spot delivery
%  * |Plan.Beams(b).spotSigma| -_SCALAR_- Spot lateral sigma (mm) at maximum of deepest Bragg peak along the optical axis
%  * |Plan.DoseGrid.size| -_SCALAR VECTOR_- [Nx, Ny, Nz] Number of pixels of the Dose map along each dimension
%  * |Plan.output_path|-_STRING_- Folder in which the results are saved
%  * |Plan.BDL| -_STRUCTURE_- Spot information from the Beam Data Library. See |getSpotFromBDL| for more information
%  * |Plan.SaveDoseBeamlets| -_BOOL_- save the dose of each beamlet in the reference frame of the CT with aperture and with the resolution of that CT: dcm (DICOM format) , sparse (sparse matrix) , false (not saved)
%
% |Dref| -_SCALAR_- Threshold dose (Gy / FRACTION) in OAR above which the dose rate condition must be respected. Voxels below the threshold dose are not included in DR condition
%
% |DoseAtPxl| -_SCALAR VECTOR_- |DoseAtPxl(pxl)| Total dose (Gy, over all fractions) delivered at pixel |pxl|
%
% |percentile| -_SCALAR_- [OPTIONAL. Default : 0.01] |P = 1 - 2.*percentile| is the percentile dose that is included in the dose rate computation
%
% |ROIName| -_STRING_- [OPTIONAL] Name of the ROI for which DRa is computed. If given, then the function displays results
%
% |plotID| -_STRUCTURE_- [OPTIONAL. Default = []] If present, display different figures with dose rate data
%    * |plotID.plots_BEV|  -_SCALAR_-  |plots_BEV(b)| is the number of the figure in which to plot the scanning trajectory of beam |b|
%    * |plotID.plot_DR|  -_INTEGER_-  Figure number where to plot the peak dose rate vs time graph
%    * |plotID.SaveHisto| -_BOOL_- TRUE = save a JSON and YAML file of the dose rate histograms
%
%
%% Output arguments
% |dr| - _SCLAR VECTOR_ - |dr(b)| Average value of the percentile dose rate (Gy/s) of all pixel in |ROImask| receiving a dose > |Dref| Gy by the beam |b|. If no voxel in the volume receives a dose > |Dref| by beam |b|, then |dr(b)| = -1
%
% |DRa| -_SCALAR_- |DRa{b}(i)| Percentile dose rate (Gy/s) at the i-th voxel of |ROImask| for the dose delivered by beam |b|. This is similar to the dose rate defined in [2].
%
% |drm| - _SCLAR VECTOR_ - |drm(b)| Median value of the percentile dose rate (Gy/s) of all pixel in |ROImask| receiving a dose > |Dref| Gy by the beam |b|. If no voxel in the volume receives a dose > |Dref| by beam |b|, then |drm(b)| = -1
%
% |DADR| -_CELL VECTOR of SCLAR VECTOR_-|DADR{b}(pxl)| Dose averaged dose rate (Gy/s) as defined in [3] at pixel |pxl| for the b-th beam
%
% |DADRm| - _SCLAR VECTOR_ - |DADR(b)| Median value of the dose averaged dose rate (Gy/s) of all pixel in |ROImask| receiving a dose > |Dref| Gy by the beam |b|. If no voxel in the volume receives a dose > |Dref| by beam |b|, then |DADRm(b)| = -1
%
% |SpotTiming| -_SCALAR VECTOR_- |SpotTiming{b}(i)| Time (ms) taken at the i-th pixel to deliver the |percentile| of the dose for the b-th beam
%
%% REFERENCE
% [1] Charlton, N., & Vukadinovi, D. (n.d.). On minimum cost local permutation problems and their application to smart meter data, 1–24. Retrieved from https://www.reading.ac.uk/web/files/maths/TechReport_2_13.pdf
% [2] Folkerts, M., Abel, E., Busold, S., Perez, J., & Vidhya Krishnamurthi, C. C. L. (2020). A framework for defining FLASH dose rate for pencil beam scanning. Medical Physics. https://doi.org/10.1002/MP.14456
% [3] van de Water, S., Safai, S., Schippers, J. M., Weber, D. C., & Lomax, A. J. (2019). Towards FLASH proton therapy: the impact of treatment planning and machine characteristics on achievable dose rates. Acta Oncologica, 58(10), 1463–1469. https://doi.org/10.1080/0284186X.2019.1627416
%
%% Contributors
% Authors : R. Labarbe, Lucian Hotoiu (open.reggui@gmail.com)

function [dr, DRa, drm, DADR, DADRm, DRAD, DRADm, SpotTiming] = getDRa(SpotTrajectoryInfo , weightIN , Pij , Plan , Dref , ROImask , DoseAtPxl , DMF, DR50, percentile , ROIName , plotID )

if nargin < 10
  percentile = [];
end
if isempty(percentile)
  percentile = 0.01; %Compute dose rate for the 98-percentile
end

if nargin >= 11
  verbose = 1;
else
  verbose =0;
end

if nargin < 12
  plotID.plots_BEV = [];
  plotID.plot_DR = [];
  plotID.pDR_D = [];
  plotID.SaveHisto = false;
end
if isempty(plotID)
  plotID.plots_BEV = [];
  plotID.plot_DR = [];
  plotID.pDR_D = [];
  plotID.SaveHisto = false;
end


sobp = SpotTrajectoryInfo.sobpPosition;
weight2spot = SpotTrajectoryInfo.weight2spot;
fractions = Plan.fractions;
Inozzle = Plan.Inozzle .* 1e-9; % nA -> A Current in nozzle
weightIN = weightIN' ./ fractions; %Compute MU PER FRACTION because the dose rate is computed per fraction. Make it a sparse matrix


dr = zeros(length(sobp),1); %dr(b) average of the percentile dose rate for the b-th beam
drm = zeros(length(sobp),1); %drm(b) Median of the percentile dose rate for the b-th beam
DADRm = zeros(length(sobp),1); %DADRm(b) Median of the dose averaged dose rate for the b-th beam

if ~isempty(ROImask)
    Tested =  DoseAtPxl > Dref; %find the voxels receiving a total dose  > than the threshold and that are inside the ROI. The DR computation will occur only in those pixels
    Tested =  Tested .* ROImask;
else
    Tested =  DoseAtPxl > Dref; %find the voxels receiving a total dose  > than the threshold and that are inside the ROI. The DR computation will occur only in those pixels
end

NbPixels = numel(find(Tested));

fprintf('-- getDRa() --- Nb of pixels tested = %d \n', NbPixels)

for b = 1:length(sobp) %Loop for each beam
    %Find an ordering of the sobp based on their neighbourhood
    fprintf('Beam %d : Using provided spot trajectory \n',b)
    SpotIdxNoProton = []; %Contains the indices of sobp{b}(i) where no proton was delivered

    OutputDir = getOutputDir(Plan.output_path , b);

    sobpSequence = SpotTrajectoryInfo.beam{b}.sobpSequence;

    if (NbPixels)
        %If the dose is higher than the threshold in the OAR, then
        %Compute the dose given by each beamlet to each pixel of ROI
        NbBeamlets = max(weight2spot(:,2) .* (weight2spot(:,1) == b)); %Number of beamlets in beam |b|
        fprintf('Beam %d : Nb beamlets : %d \n' , b , NbBeamlets)

        Dose = sparse(NbPixels,NbBeamlets); %Dose(pxl,spot) dose given at pixel |pxl| by spot number |spot|
        wSOBP = [];
        tailOfScan = [];

        for beamletIndex = 1:NbBeamlets

          %weight2spot(weightIdx,:) = [b,BeamLetNb]
          spotIndices = find((weight2spot(:,1) == b) .* (weight2spot(:,2) == beamletIndex)); %These are all the Bragg peaks contributing to beamlet |beamletIndex| of beam |b|
          T = find(Tested); %Indices of the pxl that shall be included in the computation

          %NB : There is a strange behaviour wit hthe .* multiplication of sparse matrix.
          % In some cases full(spA) .* full(spB) ~= full (spA .* spB)
          % For some spare matrices spA and spB is works fine but not for some other. I do not understand the origin of this problem.
          % To avoid the problem, we selet out of the sparse Pij matrix only the elements of the selected voxels and selected beamlets.
          % Then convert the Pij matrix into a full matrix and carry out the matrix multiplication using the full matrix
          % The result is then converted back into a sparse matrix

          if (size(weightIN(spotIndices),1)==1)
            %transpose weights
            tmp = full(Pij(:,spotIndices)) * full(weightIN(spotIndices)'); %Compute dose only for voxels in mask and above dose threshold
          else
            tmp = full(Pij(:,spotIndices)) * full(weightIN(spotIndices)); %Compute dose only for voxels in mask and above dose threshold
          end

          Dose(:,beamletIndex) = sparse(tmp(T)); %Select only the pixel above the minimum dose |Dref|


          %TODO This code does not work robustly. With some dataset, it works. With other dataset it does not work
          % It is difficult to understand the origin of the problem
          % if (size(weightIN(spotIndices),1)==1)
          %   %transpose weights
          %   Dose(:,beamletIndex) = sparse(full(Pij(T,spotIndices)) * weightIN(spotIndices)'); %Compute dose only for voxels in mask and above dose threshold
          %   max(Dose,[],'all')
          % else
          %   Dose(:,beamletIndex) = sparse(full(Pij(T,spotIndices)) * weightIN(spotIndices)); %Compute dose only for voxels in mask and above dose threshold
          % end

          if (~sum(weightIN(spotIndices)))
            %This spot has no charge. Do not include it in the sequence
            WsobpSequence = find(sobpSequence == beamletIndex); %Find the index of the spot in the sequence
            sobpSequence(WsobpSequence)=[]; %Remove the spot from the sequence
            SpotIdxNoProton(end+1) = beamletIndex; %Add one spot in the list of the spot with no protons delivered in OAR
          end

          wSOBP(beamletIndex) = sum(weightIN(spotIndices)); %Weight of the SOBP spot # beamletIndex

        end %for beamletIndex

        Dose = Dose'; %Dose(spot,pxl)

        if (verbose)
          fprintf('Beam %d : Dose per fraction per beamlet in %s : %f <= D(Gy) <= %f in %d pixels (out of %d) \n',b,ROIName,full(min(Dose,[],'all'))  ,full(max(Dose,[],'all')) , numel(find(sum(Dose,1))), NbPixels);
          fprintf('Beam %d : Dose per fraction             in %s : %f <= D(Gy) <= %f \n',b,ROIName,full(min(sum(Dose,1))),full(max(sum(Dose,1))));
        end

        switch SpotTrajectoryInfo.TimingMode
        case 'scanAlgo'
            %We received information about a scanAlgo gateway
            %Connect to the gateway to get the correct timing information
            fprintf('getDRa : Getting timing from scanAlgo gateway: %s \n',SpotTrajectoryInfo.scanAlgoGW.scanalgoGateway_IP)
            param = getMachineParam(Plan.BDL);
            [dT , TimePerSpot] = getScanalgoTiming(sobp{b} , param.MAXenergy , SpotTrajectoryInfo.scanAlgoGW , [wSOBP(sobpSequence)] , true);

        case 'Record'
            fprintf('getDRa : Getting timing from records \n')
            dT = SpotTrajectoryInfo.beam{b}.dT;
            TimePerSpot = SpotTrajectoryInfo.beam{b}.TimePerSpot;

        case 'Model'
            %Order the spot timing for the provided sequence
            %Make an estimation of the timing
            fprintf('getDRa : Getting timing from simplified machine model --- I nozzle = %f nA : --\n', Plan.Inozzle)
            index = 2:length(sobpSequence);
            dT = round(diag(SpotTrajectoryInfo.beam{b}.Nmaps.NeighbourghTimeMap(sobpSequence(index-1),sobpSequence(index))),1); %Time (ms) to go from spot spotSequence(i) to spotSequence(i+1). Round to nearest 0.1ms
            wSOBPNotNull = wSOBP;
            wSOBPNotNull(wSOBPNotNull==0)=[]; %REmove the spot with zero weight as we will not deliver them
            TimePerSpot = getSpotDeliveryTime(wSOBP , Inozzle , Plan.BDL); %Time (s) to deliver each spot. The dose is defined **PER FRACTION**
            TimePerSpot = TimePerSpot(sobpSequence) .* 1000; %Order the spot delivery time according to the delivery sequence. Convert into ms
        otherwise
            error('Unknown timing mode')
        end

        if (~isempty(plotID.plots_BEV))

          if isfield(plotID, 'pDR_D') & ~isempty(plotID.pDR_D)
              [tmp , dr(b) , DRmin, DRmax, drm(b), Tstart , Tend , DADRtmp , DADRm(b), DRADtmp, DRADm(b), SpotTiming{b}, pxlSelected , DRhisto] = DRaEstimate(sobpSequence , dT , TimePerSpot , Dose , DMF, DR50, plotID.plot_DR , percentile); %Average dose rate at several measurement point (MP) in ROI. The MP are located at the centre of the spot. Same order than |weight|
                  % If the BEV plot is displayed, then we also compute the  |DRhisto| dose rate vs dose histogram
          else
              [tmp , dr(b) , DRmin, DRmax, drm(b), Tstart , Tend , DADRtmp , DADRm(b), DRADtmp, DRADm(b), SpotTiming{b}, pxlSelected ]          = DRaEstimate(sobpSequence , dT , TimePerSpot , Dose , DMF, DR50, plotID.plot_DR , percentile); %Average dose rate at several measurement point (MP) in ROI. The MP are located at the centre of the spot. Same order than |weight|
              DRhisto = [];
          end

          %Add a star in the BEV at the projection of the voxel at which the peak dose rate vs time was measured
          figure(plotID.plot_DR)
          [i,j,k] = ind2sub(Plan.DoseGrid.size , T(pxlSelected)); %Get the x,y,z pixel index of the voxel with the dose rate plot
          k = Plan.DoseGrid.size(3)-k; %Flip the 3rd dimension
          Pmax = PXLindex2DICOM([i,j,k] , Plan.DoseGrid.resolution , Plan.CTinfo.ImagePositionPatient , Plan.Beams(b).isocenter); %Convert into physical coordinate in DICOM CS
          M = matDICOM2IECgantry(Plan.Beams(b).GantryAngle , Plan.Beams(b).PatientSupportAngle);  %4x4 rotation matrix to go from IEC gantry to DICOM : dicom = R * IEC_G
          Pgantry = M * Pmax'; %Get the cooridnate of the voxle in IEC gnatry CS

          %Add the title to the dose rate graph
          titleStr = [ROIName,'-- Pxl @ [',num2str(Pgantry(1)),' ,',num2str(Pgantry(2)),'] mm'];
          title(titleStr)
          if ~exist(OutputDir,'dir')
            mkdir (OutputDir)
          end
          saveas(plotID.plot_DR,fullfile(OutputDir,[ROIName,'_DRvtime.fig']),'fig')

          %Add a marker in the BEV graph
          figure(plotID.plots_BEV(b))
          hold on
          plot(Pgantry(1) ,Pgantry(2) , 'p', 'MarkerFaceColor','red', 'MarkerSize',15)
          saveas(plotID.plots_BEV(b),fullfile(OutputDir,[ROIName,'_BEV.fig']),'fig')

          if ~isempty(DRhisto) & isfield(plotID, 'pDR_D') & ~isempty(plotID.pDR_D)
              plot_D_DR_histo(plotID.pDR_D , DRhisto , OutputDir , ROIName , titleStr);
          end

        else
          if isfield(plotID , 'SaveHisto') & plotID.SaveHisto
            [tmp , dr(b) , DRmin, DRmax, drm(b), Tstart , Tend , DADRtmp , DADRm(b), DRADtmp, DRADm(b), SpotTiming{b} , ~ , DRhisto] = DRaEstimate(sobpSequence , dT , TimePerSpot , Dose , DMF, DR50, [] , percentile); %Average dose rate at several measurement point (MP) in ROI. The MP are located at the centre of the spot. Same order than |weight|
          else
            [tmp , dr(b) , DRmin, DRmax, drm(b), Tstart , Tend , DADRtmp , DADRm(b), DRADtmp, DRADm(b), SpotTiming{b} ] = DRaEstimate(sobpSequence , dT , TimePerSpot , Dose , DMF, DR50, [] , percentile); %Average dose rate at several measurement point (MP) in ROI. The MP are located at the centre of the spot. Same order than |weight|
          end
        end

        if isfield(plotID , 'SaveHisto') & plotID.SaveHisto
            if ~exist(OutputDir,'dir')
              mkdir (OutputDir)
            end
            %WriteYaml(fullfile(OutputDir , [ROIName,'_DRvsD_hist.yaml']) , DRhisto); %Save the dose rate histogram data into a Yaml file
            savejson('',DRhisto,fullfile(OutputDir , [ROIName,'_DRvsD_hist.json']));
        end


        if (~isempty(plotID.plots_BEV))
          %Overlay the spot trajectory on the BEV
          plotTrajectory(sobp{b} , sobpSequence , plotID.plots_BEV(b) , Tstart , SpotIdxNoProton, b , TimePerSpot);
        end

        if (Plan.SaveDoseBeamlets & ~isempty(plotID.plots_BEV))
         %Save the text file only if we need to save the bemalet dose map
         %AND if the plots should be displayed. Otherwise the file will be saved at
         %each iteration of the optimiser

          %Check tha tthe folder exists
          folderName = fullfile(OutputDir,'CEF_beam','Outputs','SpotDoseInCT');
          if ~exist(OutputDir,'dir')
            mkdir (OutputDir)
          end
          %Save a text file with the spot timing information
          fileName = fullfile(folderName,'SpotTiming.txt');
          saveSpotTiming(fileName , sobpSequence , sobp{b} , Tstart , Tend);

        end

        [i,j]= find(Tested);
        %[i,j]= ind2sub(size(DoseAtPxl),T);
        DRa{b}  = sparse(full(i),full(j),tmp,size(DoseAtPxl,1),size(DoseAtPxl,2));
        DADR{b} = sparse(full(i),full(j),full(DADRtmp),size(DoseAtPxl,1),size(DoseAtPxl,2));

        if numel(DRADtmp > 1)
          %If the DRAD has been computed, then process it
          DRAD{b} = sparse(full(i),full(j),DRADtmp,size(DoseAtPxl,1),size(DoseAtPxl,2));
        else
          DRAD{b}= 0;
        end

        if (verbose)
          fprintf('Beam %d : Percentile Dose Rate in %s : %3.3g <= DRa = %3.3g (Gy/s) // DRm = %3.3g (Gy/s) <= %3.3g \n',b,ROIName,full(DRmin), full(dr(b)), full(drm(b)), full(DRmax));
          fprintf('Beam %d : DADR                 in %s : %3.3g <= DADR(Gy/s)= %3.3g // DADRm = %3.3g (Gy/s)  <= %3.3g \n',b,ROIName,round(full(min(DADRtmp))),round(full(mean(DADRtmp))),full(DADRm(b)),round(full(max(DADRtmp))));
          if isfield(plotID, 'pDR_D') & ~isempty(plotID.pDR_D)
            fprintf('Beam %d : DRAD                 in %s : %3.3g <= DRAD(Gy)= %3.3g // DRADm = %3.3g (Gy)  <= %3.3g \n',b,ROIName,round(full(min(DRADtmp))),round(full(mean(DRADtmp))),full(DRADm(b)),round(full(max(DRADtmp))));
          end
        end
    else
      % There no dose delivered in the OAR for any beam
      DRa{b} = sparse(size(DoseAtPxl,1),size(DoseAtPxl,2));
      DADR{b} = sparse(size(DoseAtPxl,1),size(DoseAtPxl,2));
      DRAD{b} = sparse(size(DoseAtPxl,1),size(DoseAtPxl,2));
      SpotTiming{b}= sparse(size(DoseAtPxl,1),size(DoseAtPxl,2));

      dr(b)  = -1; %There is no voxel receiving a dose above the threshold. Define DRa < 0 to indicate that Dra is not defined in this case
      drm(b) = -1;
      DADRm(b) = -1;
      DRADm(b) = -1;
      if (verbose)
        fprintf('Beam %d : Percentile Dose Rate %s : %3.3g <= DRa = %3.3g (Gy/s) // DRm = %3.3g (Gy/s) <= %3.3g \n',b,ROIName,full(dr(b)), full(dr(b)), full(drm(b)), full(dr(b)));
      end
    end

  end %for b

end
