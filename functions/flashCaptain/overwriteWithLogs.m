%% overwriteWithLogs
% Overwrite the spot info in the |Plan| with the info read from |spots| (usually loaded from the irradiation logs).
% Use |spots.name| to identify the beam to update in |Plan.Beams|
%
%% Syntax
% |Plan = overwriteWithLogs(Plan, spots)|
%
%
%% Description
% |Plan = overwriteWithLogs(Plan, spots)| Description
%
%
%% Input arguments
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are stored.
%
% |spots| -_STRUCTURE_- Record from the treatment plan
%     * |spots.name| -_STRING_- Name of the beam to which the log record refers
%     * |spots.spots.energy| -_SCALAR_- Energy (Mev) of the layer
%     * |spots.spots.xy| - _SCALAR VECTOR_ - Average spot position (x,y) at isocenter over the delivery of the s-th spot
%     * |spots.spots.weight| - _SCALAR VECTOR_ - Monitor unit of the s-th spot
%     * |spots.spots.time| - _SCALAR VECTOR_ - Time (in ms) at the begining of the delivery of the s-th spot
%     * |spots.spots.duration| - _SCALARVECTOR_ - Time (in ms) at the end of the delivery of the s-th spot
%
%% Output arguments
%
% |Plan| - _struct_ - Updated plan information.
%   * |Plan.Beams(idxBeam).Layers.Energy| -_SCALAR_- Energy (Mev) of the layer
%   * |Plan.Beams(idxBeam).Layers.nominalSpotPosition| - _SCALAR VECTOR_ - Average spot position (x,y) at isocenter over the delivery of the s-th spot
%   * |Plan.Beams(idxBeam).Layers.SpotPositions| - _SCALAR VECTOR_ - Average spot position (x,y) at isocenter over the delivery of the s-th spot
%   * |Plan.Beams(idxBeam).Layers.SpotWeights| - _SCALAR VECTOR_ - Monitor unit of the s-th spot
%   * |Plan.Beams(idxBeam).Layers.time| - _SCALAR VECTOR_ - Time (in ms) at the begining of the delivery of the s-th spot
%   * |Plan.Beams(idxBeam).Layers.duration| - _SCALARVECTOR_ - Time (in ms) at the end of the delivery of the s-th spot
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function Plan = overwriteWithLogs(Plan, spots)

  %Identify the beam in |Plan.Beams| with a name matching |spots.name|
  idxBeam = find(strcmp([Plan.Beams(:).name] , spots.name));

% figure(20)
% plot(Plan.Beams(idxBeam).Layers.nominalSpotPosition(:,1),Plan.Beams(idxBeam).Layers.nominalSpotPosition(:,2),'ob')

  %Update the spot info in |Plan.Beams|
  Plan.Beams(idxBeam).Layers.Energy = spots.spots(1).energy;
  Plan.Beams(idxBeam).Layers.nominalSpotPosition = spots.spots(1).xy;
  Plan.Beams(idxBeam).Layers.SpotPositions = spots.spots(1).xy;
  Plan.Beams(idxBeam).Layers.SpotWeights = spots.spots(1).weight';
  Plan.Beams(idxBeam).Layers.time = spots.spots(1).time .* 1000; %convert into ms
  Plan.Beams(idxBeam).Layers.duration = spots.spots(1).duration .* 1000; %convert into ms

% figure(20)
% hold on
% plot(Plan.Beams(idxBeam).Layers.nominalSpotPosition(:,1),Plan.Beams(idxBeam).Layers.nominalSpotPosition(:,2),'+r')
% pause


end
