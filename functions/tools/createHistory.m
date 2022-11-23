%% createHistory
% Create the time history for dose and dose rate from the snapshots
%
%% Syntax
% | [DosePoint , DoseRatePoint , timeAxis] = createHistory(timeAxis , SpotTimeStamp, SpotTimeDuration , doseSnaps , DoseRateSnaps)|
%
%
%% Description
% | [DosePoint , DoseRatePoint , timeAxis] = createHistory(timeAxis , SpotTimeStamp, SpotTimeDuration , doseSnaps , DoseRateSnaps)| Description
%
%
%% Input arguments
%  |timeAxis| -_SCALAR VECTOR_- |timeAxis(i)| the i-th time point (ms) at which the dose rate profile is defined
%
% |SpotTimeStamp| - _SCALAR VECTOR_ - |SpotTimeStamp(t)| The time (ms) at the begining of the delivery of the t-th spot
%
% |SpotTimeDuration| - _SCALAR VECTOR_ - |SpotTimeStamp(t)| duration (ms) of the delivery of the t-th spot
%
% |Snapshot| - _SCALAR MATRIX_ - |Snapshot(t)| the cummulated dose (Gy) delivered at pixel at time t
%
% |DoseRateSnaps| - _SCALAR MATRIX_ - |DoseRateSnaps(t)| the dose rate (Gy/s) delivered at pixel (x,y,z) at time t
%
%
%% Output arguments
%
%  |timeAxis| -_SCALAR VECTOR_- |timeAxis(i)| the i-th time point (ms) at which the dose rate profile is defined
%
%  |DoseRatePoint| -_SCALAR VECTOR_- |DoseRatePoint(i)| instantaneous dose rate (Gy/s) at time |timeAxis(i)|
%
%  |DosePoint| -_SCALAR VECTOR_- |DosePoint(i)| Cumulative dose (Gy)  up to time |timeAxis(i)|
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [DosePoint , DoseRatePoint , timeAxis] = createHistory(timeAxis , SpotTimeStamp, SpotTimeDuration , doseSnaps , DoseRateSnaps)

DosePoint = zeros(1,length(timeAxis));
DoseRatePoint = zeros(1,length(timeAxis));

for i = 1:length(SpotTimeStamp)
  %Define dose rate during beam ON window
  wBeamOn = find((timeAxis >= SpotTimeStamp(i)) .* (timeAxis <= SpotTimeStamp(i) + SpotTimeDuration(i))); %The multiplication is a AND logical operator
  if (~isempty(wBeamOn))
      DoseRatePoint(wBeamOn) =  DoseRateSnaps(i); %Constant dose rate

      %Dose increase linearly during beam on
      if(wBeamOn(end) - wBeamOn(1) > 0)
        %The dose is delivered over several points
        DosePoint(wBeamOn) = DosePoint(wBeamOn(1))+ (doseSnaps(i)-DosePoint(wBeamOn(1))).*(wBeamOn-wBeamOn(1))./(wBeamOn(end)-wBeamOn(1));
      else
        %The dose is delivered over one time point
        DosePoint(wBeamOn) = DosePoint(wBeamOn(1))+ doseSnaps(i);
      end

      %Accumulate dose after complete delivery of the spot
      wBeamOn = find(timeAxis > (SpotTimeStamp(i)+SpotTimeDuration(i)));
      DosePoint(wBeamOn) = doseSnaps(i).*ones(1,length(wBeamOn));
    end
end

end
