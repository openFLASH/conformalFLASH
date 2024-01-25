%% mergeSeqSpots
% Reduce the number of spots by merging SEQUENTIALS spots that are delivered at the same [X,Y] location (within 0.5mm)
%
%% Syntax
% |spots = mergeSeqSpots(spots)|
%
%
%% Description
% |spots = mergeSeqSpots(spots)| Description
%
%
%% Input arguments
% |spots| -_STRUCTURE_- [OPTIONAL. If absent, the spot information is obtained from plan]
% * |spots.charge(s)| - _SCALAR_ - Measured electrical charge of the s-th spot of the l-th energy layer
% * |spots.timeStart(s)| - _SCALAR_ - Time (in [s]) at the begining of the delivery of the s-th spot of the l-th energy layer. Time resolution in multiple of the log period
% * |spots.timeStop(s)| - _SCALAR_ - Time (in [s]) at the end of the delivery of the s-th spot of the l-th energy layer. Time resolution in multiple of the log period
% * |spots.time(s)| - _SCALAR_ - Time (in [s]) of the 'center of mass' of the spot (dose weighted time average)
% * |spots.duration(s)| - _SCALAR_ - Duration (in [s]) of the spot delivery. Time resolution in multiple of the log period
% * |spots.metersetRate(s)| - _SCALAR_ - Dose rate (MU/s) of the s-th spot
% * |spots.effectiveDuration(s)| - _SCALAR_ - Durartion (in [s]) of the spot delivery with a higher resolution than the log period. Estimated from the average current during spot delivery
% * |spots.xyIC(s,:)| - _SCALAR VECTOR_ - Average spot position (x,y) on the IC over the delivery of the s-th spot of the l-th energy layer. The coordinate system is IEC-GANTRY. The IC_Offset has already been corrected.
% * |spots.nb_protons(s)| - _SCALAR_ - Number of proton in the s-th spot of the l-th energy layer
% * |spots.IC_Offset(s,:)| - _SCALAR VECTOR_ - Offset (mm) of the IC coordinate system with respect to IECgantry for the s-th spot of the l-th energy layer
% * |spots.SM_Offset(s,:)| - _SCALAR VECTOR_ - Offset (mm) of the SM coordinate system with respect to IECgantry for the s-th spot of the l-th energy layer
% * |spots.xy(s,:)| - _SCALAR VECTOR_ - Average spot position (x,y) at isocenter over the delivery of the s-th spot of the l-th energy layer. The coordinate system is IEC-GANTRY.
% * |spots.weight(s)| - _SCALAR_ - MU of the s-th spot of the l-th energy layer
%
%
%% Output arguments
%
% |res| - _STRUCTURE_ -  Description
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function spots = mergeSeqSpots(spots)

  tol = 0.5; %mm If two spots are more than 0.5mm apart, they are probably different spot

  if isfield(spots.spots , 'metersetRate') && ~isfield(spots.spots , 'effectiveDuration')
    %Compute effective duration from metersetRate
    fprintf('Computing effective duration from metersetrate \n')
    spots.spots.effectiveDuration = spots.spots.weight ./ spots.spots.metersetRate; %s
  end


  spots2 = spots;
  spots2.spots.xy = spots.spots.xy(1,:);
  %spots2.spots.xyIC = spots.spots.xyIC(1,:);
  spots2.spots.weight = spots.spots.weight(1);
  %spots2.spots.charge = spots.spots.charge(1);
  %spots2.spots.nb_protons = spots.spots.nb_protons(1);


  spots2.spots.time = spots.spots.time(1);
  spots2.spots.timeStart = spots.spots.timeStart(1);
  spots2.spots.timeStop = spots.spots.timeStop(1);
  spots2.spots.duration = spots.spots.duration(1);
  spots2.spots.effectiveDuration = spots.spots.effectiveDuration(1);

  %spots2.spots.painting = spots.spots.painting(1);
  spots2.spots.spot_id = spots.spots.spot_id(1);
  %spots2.spots.SM_Offset = spots.spots.SM_Offset(1);
  %spots2.spots.IC_Offset = spots.spots.IC_Offset(1);
  %spots2.spots.tuning = spots.spots.tuning(1);

  realStart = spots.spots.time(1) - spots.spots.effectiveDuration(1) /2; %Real start time (us) of the first spot

for idx = 2:numel(spots.spots.weight)

    if norm(spots2.spots.xy(end,:) - spots.spots.xy(idx,:)) < tol
      %These are the same spots. Merge them
      spots2.spots.weight(end) = spots2.spots.weight(end) + spots.spots.weight(idx);
      spots2.spots.duration(end) = spots2.spots.duration(end) + spots.spots.duration(idx);
      spots2.spots.timeStop(end) =  spots.spots.timeStop(idx);
      %spots2.spots.charge(end) = spots2.spots.charge(end) + spots.spots.charge(idx);
      %spots2.spots.nb_protons(end) = spots2.spots.nb_protons(end) + spots.spots.nb_protons(idx);

      realStop = spots.spots.time(idx) + spots.spots.effectiveDuration(idx) /2; %Real start time of the idx-th spot
      spots2.spots.time(end) = (realStart + realStop) ./2;
      spots2.spots.effectiveDuration(end) = realStop - realStart;

    else
      %This is a new spot. Add it to the list
      spots2.spots.xy(end+1,:)     = spots.spots.xy(idx,:);
      %spots2.spots.xyIC(end+1,:)     = spots.spots.xyIC(idx,:);
      spots2.spots.weight(end+1,1)   = spots.spots.weight(idx);
      %spots2.spots.nb_protons(end+1,1)   = spots.spots.nb_protons(idx);
      %spots2.spots.charge(end+1,1)   = spots.spots.charge(idx);

      spots2.spots.timeStart(end+1,1) = spots.spots.timeStart(idx);
      spots2.spots.timeStop(end+1,1) = spots.spots.timeStop(idx);
      spots2.spots.duration(end+1,1) = spots.spots.duration(idx);

      realStart = spots.spots.time(idx) - spots.spots.effectiveDuration(idx) /2; %Real start time of the spot
      spots2.spots.time(end+1,1)     = spots.spots.time(idx);
      spots2.spots.effectiveDuration(end+1,1) = spots.spots.effectiveDuration(idx);

      %spots2.spots.painting(end+1,1) = spots.spots.painting(idx);
      spots2.spots.spot_id(end+1,1) = spots.spots.spot_id(idx);
      %spots2.spots.SM_Offset(end+1,1) = spots.spots.SM_Offset(idx);
      %spots2.spots.IC_Offset(end+1,1) = spots.spots.IC_Offset(idx);
      %spots2.spots.tuning(end+1,1) = spots.spots.tuning(idx);
    end

  end

%Compute meterset rate
spots2.spots.metersetRate = spots2.spots.weight ./ (spots2.spots.effectiveDuration/1e6); % in [MU/s]
spots = spots2;

end
