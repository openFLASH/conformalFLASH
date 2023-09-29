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
%     * |spots.name| -_STRING_- Name of the beam to which the log record refers
%     * |spots.spots.xy| - _SCALAR VECTOR_ - Average spot position (x,y) at isocenter over the delivery of the s-th spot
%     * |spots.spots.weight| - _SCALAR_ - Monitor unit of the s-th spot
%     * |spots.spots.time| - _SCALAR_ - Time at the begining of the delivery of the s-th spot
%     * |spots.spots.duration| - _SCALAR_ - Time at the end of the delivery of the s-th spot of the l-th energy layer
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

  spots2 = spots;
  spots2.spots.xy = spots2.spots.xy(1,:);
  spots2.spots.weight = spots2.spots.weight(1);
  spots2.spots.time = spots2.spots.time(1);
  spots2.spots.duration = spots2.spots.duration(1);

  for idx = 2:numel(spots.spots.weight)

    if norm(spots2.spots.xy(end,:) - spots.spots.xy(idx,:)) < tol
      %These are the same spots. Merge them
      spots2.spots.weight(end) = spots2.spots.weight(end) + spots.spots.weight(idx);
      spots2.spots.duration(end) = spots2.spots.duration(end) + spots.spots.duration(idx);
    else
      %This is a new spot. Add it to the list
      spots2.spots.xy(end+1,:)     = spots.spots.xy(idx,:);
      spots2.spots.weight(end+1,1)   = spots.spots.weight(idx);
      spots2.spots.time(end+1,1)     = spots.spots.time(idx);
      spots2.spots.duration(end+1,1) = spots.spots.duration(idx);
    end

  end

spots = spots2;

end
