%% getSweepTime
% Compute the sweep time between sucessive spot from the spot time stamps and the spot delivery duration
%
%% Syntax
% |Tsweep = getSweepTime(time, duration)|
%
%
%% Description
% |Tsweep = getSweepTime(time, duration)| Description
%
%
%% Input arguments
%|time| -_SCALAR VECTOR_- |time(i)| Time stamp (s) of the i-th spot
%
%|duration| -_SCALAR VECTOR_- |duration(i)| Duration (ms) of the delivery of the i-th spot
%
%
%% Output arguments
%
% |Tsweep| -_SCALAR VECTOR_- |Tsweep(i)| Sweep time (ms) to go from the (i-1)-th spot to the i-th spot
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)


%-----------------------------

function [Tsweep , TimePerSpot] = getSweepTime(time, duration)

  Tsweep = diff(time).*1000 - duration(2:end);
  Tsweep = round(Tsweep);
  Tsweep(Tsweep<0) = 0;

end
