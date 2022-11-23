%% plotDoseRate
% Plot the dose and dose rate as a function of time
%
%% Syntax
% |res = help_header(im1,im2)|
%
%
%% Description
% |res = help_header(im1,im2)| Description
%
%
%% Input arguments
% |spotTimingStart| -_SCALAR MATRIX_-  |timeAxis(spot)| time (ms) at which |spot| is delivered
%
% |Dose| -_SCALAR MATRIX_- |Dose(spot)| Dose (Gy) delivered by spot number |spot|
%
% |TimePerSpot| -_SCALAR VECTOR_- |TimePerSpot(spot)| Time (ms) required to deliver all the protons in spot # |spot|
%
% |fig| -_INTEGER_- Figure number where to plot
%
%
%% Output arguments
%
% None
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function plotDoseRate(spotTimingStart , Dose , TimePerSpot , figNb)

  timeAxis = spotTimingStart;
  wzeroDose = find(Dose==0);
  Dose(wzeroDose) = []; %Remove beamlet delivering no dose
  TimePerSpot(wzeroDose) = [];
  spotTimingStart(wzeroDose) = [];
  DoseRatePoint = 1000 .* Dose' ./ double(TimePerSpot); %Convert time from ms to s in order to get dose rate in Gy/s
          %|Dose| is sparse and |TimePerSpot| is single. Need to cast |TimePerSpot| to doble to be able to divide the sparse matrix
  DosePoint = cumsum(Dose); %The cumulated dose at each time point |spotTimingStart|

  %Take into account the duration of pulses for proper display
  timeAxis = 1:min(TimePerSpot)./5:round(spotTimingStart(end)+TimePerSpot(end)); %Create a time vector covering the whole delivery
  [DosePoint , DoseRatePoint , timeAxis] = createHistory(timeAxis , spotTimingStart, TimePerSpot , DosePoint , DoseRatePoint);

  figure(figNb)
  yyaxis left
  plot(timeAxis,DosePoint)
  xlabel('Time (ms)')
  ylabel('Dose (Gy)')
  axis([min(timeAxis) , max(timeAxis) , 0 , 1.1 .* max(DosePoint,[],'all')])
  grid on

  yyaxis right
  plot(timeAxis, DoseRatePoint)
  ylabel('Dose rate (Gy/s)')
  axis([min(timeAxis) , max(timeAxis) , 0 , 1.1 .* max(DoseRatePoint,[],'all')])
  grid on

end
