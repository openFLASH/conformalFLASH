%% DRaEstimate
% Compute the different dose rates inside the organ receiving |Dose|.
%
%% Syntax
% |[DoseRate , DRa , DRmin, DRmax , DRm , Tstart , Tend , DADR , DADRm, DRAD, DRADm, SpotTiming , pxlSelected , DRhisto] = DRaEstimate(spotSequence , dT , TimePerSpot , Dose , DMF, DR50, fig , percentile)|
%
%
%% Description
% |[DoseRate , DRa , DRmin, DRmax , DRm , Tstart , Tend , DADR , DADRm, DRAD, DRADm, SpotTiming , pxlSelected , DRhisto] = DRaEstimate(spotSequence , dT , TimePerSpot , Dose , DMF, DR50, fig , percentile)| Description
%
%
%% Input arguments
%
% |spotSequence| -_SCALAR VECTOR_- The s-th spot of the trajectory  corresponds to |Dose(spotSequence(s))|
%
% |dT| -_SCALAR VECTOR_- |dT(s)| Time (ms) required to move from the sth-1 spot to the s-th spot in the ordered sequence
%
% |TimePerSpot| -_SCALAR VECTOR_- |TimePerSpot(s)| Time (ms) required to deliver all the protons in s-th spot of the seqeunce
%
% |Dose| -_SCALAR MATRIX_- |Dose(spot,pxl)| Dose (Gy) delivered by spot number |spot| to the pixel |pxl|. Only the pixel receiving a dose above the threshold dose are provided to the function
%
% |DMF| -_SCALAR_- [OPTIONAL. Default = []] Dose modifying factor of flash for this organ.
%
% |DR50| -_SCALAR_- [OPTIONAL. Default = []] Gy/s 50% dose rate of the sigmoid probability function.
%
% |fig| -_INTEGER_- Figure number where to plot the peak dose rate vs time graph
%
% |percentile| -_SCALAR_- [OPTIONAL. Default : 0.01] |P = 1 - 2.*percentile| is the percentile dose that is included in the dose rate computation
%
%% Output arguments
% |DoseRate| - _SCLAR VECTOR_ - |DoseRate(pxl)| Average dose rate (Gy/s) of the 98-percentileat the pxl-th measurement point.
%               The dose rate = 0 Gy/s at pixel were no dose was delivered (0Gy).
%               The 98-percentile dose rate is the ratio of 98% of the dose delivered to a pixel divided by the time taken to deliver this 98% of dose
%
% |DRa| -_SCALAR_- Average dose rate (Gy/s) of all pixel.
%
% |DRmin| -_SCALAR_- Minimum dose rate (Gy/s) of all pixel.
%
% |DRmax|  -_SCALAR_- MAximum dose rate (Gy/s) of all pixel.
%
% |DRm| -_SCALAR_- Median dose rate (Gy/s) of all pixel
%
% |Tstart| -_SCALAR VECTOR_- |Tstart(i)| Timing (ms) at the begining of the i-th spot delivery in |spotSequence|
%
% |Tend| -_SCALAR VECTOR_- |Tend(i)| Timing (ms) at the end of the i-th spot delivery in |spotSequence|
%
% |DADR(pxl)| Dose averaged dose rate (Gy/s) as defined in [3] at pixel |pxl|
%
% |DADRm| -_SCALAR_- Median of the dose averaged dose rate (Gy/s) of all pixel.
%
% |SpotTiming| -_SCALAR VECTOR_- |SpotTiming(i)| Time (ms) taken at the i-th pixel to deliver the |percentile| of the dose
%
% |pxlSelected| -_INTEGER_- Index of the pixel in |Dose(: , pxlSelected)| for which the dose rate vs time plot was drawn
%
% |DRhisto| -_STRUCT_- Structure with information to reconstruct the dose vs dose rate histogram
%
%% REFERENCE
% [1] Charlton, N., & Vukadinovi, D. (n.d.). On minimum cost local permutation problems and their application to smart meter data, 1–24. Retrieved from https://www.reading.ac.uk/web/files/maths/TechReport_2_13.pdf
% [2] Folkerts, M., Abel, E., Busold, S., Perez, J., & Vidhya Krishnamurthi, C. C. L. (2020). A framework for defining FLASH dose rate for pencil beam scanning. Medical Physics. https://doi.org/10.1002/MP.14456
% [3] van de Water, S., Safai, S., Schippers, J. M., Weber, D. C., & Lomax, A. J. (2019). Towards FLASH proton therapy: the impact of treatment planning and machine characteristics on achievable dose rates. Acta Oncologica, 58(10), 1463–1469. https://doi.org/10.1080/0284186X.2019.1627416
%
%% Contributors
% Authors : R. Labarbe, Lucian Hotoiu (open.reggui@gmail.com)

function [DoseRate , DRa , DRmin, DRmax , DRm , Tstart , Tend , DADR , DADRm, DRAD, DRADm, SpotTiming , pxlSelected , DRhisto] = DRaEstimate(spotSequence , dT , TimePerSpot , Dose , DMF, DR50, fig , percentile)

  if nargin < 5
    DMF = [];
  end
  if nargin < 6
    DR50 = [];
  end

    if nargin < 7
      fig = 80;
    end
    if nargin < 8
      percentile = [];
    end
    if isempty(percentile)
      percentile = 0.01; %Compute dose rate for the 98-percentile
    end
    if (percentile < 0) || (percentile > 1)
      percentile = 0.01; %Default value of percentile if an invalid value is given
    end

    NbSpots = length(spotSequence); %# of spots
    NbPxl = size(Dose , 2); %# Nb de measurement points

    SweepTime = [0 , cumsum(dT)']; %SweepTime(i) is time at which spot spotSequence(i) was delivered (still need to add the beam delivery time)
    SweepTime = SweepTime + 1; %the clock start at '1', so that the 0 is reserved for "no dose delivered"
    DeliveryTime =  cumsum(TimePerSpot);  %Cumulated time (ms) required to deliver the protons for each spot
    spotTimingStart = double(SweepTime + [0 , DeliveryTime(1:end-1)]); %Time at the begining of the spot delivery
    spotTimingStop = double(SweepTime + DeliveryTime); %Time at the end of the spot delivery
                  %Cast to double to allow division of sparse matrix by time

    if ~isempty(fig)
      %ROIname is provided. Therefore, we ant to see a plot
      [~ , pxlSelected] = max(sum(Dose,1)); %Plot the DR vs time profile for the pixel with highest dose
      plotDoseRate(spotTimingStart , Dose(spotSequence , pxlSelected) , TimePerSpot , fig );
    else
      pxlSelected = [];
    end

    %Re-order the dose in the same order as the spotTimingStart and spotTimingStop vectors
    [DoseRate  , SpotTiming] = GetPercentileDR(Dose(spotSequence,:) , spotTimingStart, spotTimingStop , percentile);

    %Compute dose averaged dose rate
    DADR = doseAveragedDoseRate(Dose(spotSequence,:) , spotTimingStart, spotTimingStop);

    if ~isempty(DMF) & ~isempty(DR50)
      DRAD = dRaDs(Dose(spotSequence,:) , spotTimingStart, spotTimingStop, DMF, DR50);
    else
      DRAD = [];
    end


    %Sort the start and stop time of the voxel delivery
    Nspot = numel(TimePerSpot);
    Tstart = zeros(1,Nspot);
    Tstart(spotSequence) = spotTimingStart - 1;
    Tend  = zeros(1,Nspot);
    Tend(spotSequence) = spotTimingStop - 1;

    fprintf('Percentile for dose rate computation : %d %% \n',round(100.*(1-2.*percentile)))
    fprintf('Spot  duration         (ms): %3.3g <= (Tav=%3.3g) <= %3.3g \n',min(TimePerSpot),mean(TimePerSpot),max(TimePerSpot))
    wNotnan = find(~isnan(SpotTiming));
    fprintf('Pixel Irradiation time (ms): %3.3g <= (Tav=%3.3g) <= %3.3g \n',min(SpotTiming),mean(SpotTiming(wNotnan)),max(SpotTiming))
    fprintf('Irradiation time of OAR(s): %3.3g \n', spotTimingStop(end)*1e-3) %includes sweep time in between spots and spot irradiation time


    %Compute average dose rate at each spot
    %=========================================
    DRa = mean(DoseRate); %Average dose rate of pixels with dose above threshold (Dref) dose
    DRmin = min(DoseRate);
    DRmax = max(DoseRate);
    DRm = median(DoseRate);
    DADRm = median(DADR);
    if ~isempty(DRAD)
      DRADm = median(DRAD);
    else
      DRADm = 0;
    end

    if nargout >= 14
      %Compute the dose rate histogram only if required because it takes time
      DRhisto = peakDoseRateHistogram(Dose(spotSequence,:) , spotTimingStart, spotTimingStop);
    end

end


% INPUT
% |Dose| -_SCALAR MATRIX_- |Dose(spot,pxl)| Dose (Gy) delivered by spot number |spot| to the pixel |pxl|
%
% |spotTimingStart| -_SCALAR VECTOR_- |spotTimingStart(i)| Timing (ms) at the begining of the i-th spot delivery
%
% |spotTimingStop| -_SCALAR VECTOR_- |spotTimingStop(i)| Timing (ms) at the end of the i-th spot delivery
%
% |percentile| -_SCALAR_-  Define the lower and higher dose percentile to ignore to compute dose rate and delivery time (eg. 0.01 for the 98-percentile)
%
%
% OUTPUT
%
% |DR| -_SCALAR VECTOR_- |DR(pxl)| Dose rate (Gy/s) to deliver the |percentile| of the dose to the pixel |pxl|
%
% |SpotTiming| -_SCALAR VECTOR_- |SpotTiming(i)| Time (ms) taken at the i-th pixel to deliver the |percentile| of the dose

function [DR ,  SpotTiming] = GetPercentileDR(Dose , spotTimingStart, spotTimingStop , percentile)

    cumDose = cumsum(Dose,1); % cumDose(spot,pxl) cumulated dose at the end of delivery of spot |spot| to pixel # |pxl|
    NbSpots = size(Dose,1);
    DoseMAX = repmat(cumDose(end,:) , NbSpots , 1 ); %|DoseMAX(spot,pxl)| Total dose delivered at pixel # |pxl|

    [~, startTimeIndex] = max( (cumDose ./ DoseMAX) >= percentile,[],1); % index of the first spot for which the cumulated delivered dose is above the lower percentile
    startTime = spotTimingStart(startTimeIndex);

    [~, endTimeIndex] = max( (cumDose ./ DoseMAX) >= (1-percentile),[],1); % index of the first spot for which the cumulated delivered dose is above the higher percentile
    endTime = spotTimingStop(endTimeIndex);

    SpotTiming = endTime - startTime;

    wZero = find(SpotTiming == 0);
    DR = 1000 .*  (1 - 2.* percentile) .* cumDose(end,:) ./ SpotTiming ; %Dose rate (Gy/s) to deliver the percentile dose
    DR(wZero) = 0;

end


%====================================
% Compute the dose averaged dose rate [3]
% INPUT
% |Dose| -_SCALAR MATRIX_- |Dose(spot,pxl)| Dose (Gy) delivered by spot number |spot| to the pixel |pxl|
%
% |spotTimingStart| -_SCALAR VECTOR_- |spotTimingStart(i)| Timing (ms) at the begining of the i-th spot delivery
%
% |spotTimingStop| -_SCALAR VECTOR_- |spotTimingStop(i)| Timing (ms) at the end of the i-th spot delivery
%
% |percentile| -_SCALAR_-  Define the lower and higher dose percentile to ignore to compute dose rate and delivery time (eg. 0.01 for the 98-percentile)
%
% OUTPUT
% |DADR| -_SCALAR VECTOR_- |DADR(pxl)| Dose averaged dose rate (Gy/s) as defined in [3] delivered to the pxl-th pixel
%====================================
function DADR = doseAveragedDoseRate(Dose , spotTimingStart, spotTimingStop)
    nBPxl = size(Dose,2);
    dT = spotTimingStop - spotTimingStart; %dT(spot) Time to deliver i-th spot
    DR = Dose ./ repmat(dT', 1,nBPxl); % DR(spot,pxl) Dose rate at pxl-th voxel and for delivery of spot
    Dtot = sum(Dose,1); %Dtot(pxl) Total dose delivered at pixel |pxl|
    DADR = 1000 .* sum(DR .* Dose ,1) ./ Dtot; %Convert into Gy/s from ms
end
