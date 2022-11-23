%% clip2SpikeResolution
% Group the steps with small width so that the new step has a width at least equal to the resolution
%
%% Syntax
% |[r_min , r_max , h ] = clip2SpikeResolution(a_max , a_min , h_step , resolution)|
%
%
%% Description
% |[r_min , r_max , h ] = clip2SpikeResolution(a_max , a_min , h_step , resolution)| Description
%
%
%% Input arguments
% |a_max| - _SCALAR VECTOR_ - Max apothem (mm) of the L-th polygon (apothem of the L-th hexagon/square)
%
% |a_min| - _SCALAR VECTOR_ - Min apothem (mm) of the L-th polygon (apothem of the (L+1)-th hexagon/square)
%
% |h_step| - _SCALAR VECTOR_ - Height (mm) of the L-th polygon from the base
%
% |resolution| -_SCALAR_- Minimum size (mm) of a terrace
%
%% Output arguments
%
% |r_max| - _SCALAR VECTOR_ - Max apothem (mm) of the L-th polygon (apothem of the L-th hexagon/square)
%
% |r_min| - _SCALAR VECTOR_ - Min apothem (mm) of the L-th polygon (apothem of the (L+1)-th hexagon/square)
%
% |h| - _SCALAR VECTOR_ - Height (mm) of the L-th polygon from the base
%
% |hidx| -_INTEGER VECTOR_- |hidx(i)| Index in h_step of the step |h(i)|
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [r_min , r_max , h , hidx] = clip2SpikeResolution(a_max , a_min , h_step , resolution)

  width = a_max  - a_min; %Width of each step
  cWidth = cumsum(width); %cumulative distance from edge to the current step
  Z = floor(cWidth ./ resolution);
  wWidthJmp = find(diff([0 , Z])); %Index of the steps that are at multiple of the resolution
            %Add zero as the first term of the cumsum series. The first step occurs after a cumsum of zero
            %diff reduces the vector length by 1. Adding the zero allows keeping correct track of the indices

  r_min  = zeros(1 , numel(wWidthJmp));
  r_max   = zeros(1 , numel(wWidthJmp));
  h  = zeros(1 , numel(wWidthJmp));
  hidx = zeros(1 , numel(wWidthJmp));
  wWidthJmp = [0 , wWidthJmp];

  for idx = 2:numel(wWidthJmp)

    idxSteps = wWidthJmp(idx-1)+1:wWidthJmp(idx); %Indices of the steps to be gathered into a single step with width |resolution|
    H_weighted = sum(h_step(idxSteps) .* width(idxSteps) ./ sum(width(idxSteps)) );  %Weighted height of the step
    [~, Ti] = min((H_weighted - h_step(idxSteps)).^2); %Index of the step with a height closezr to |H_weighted|

    h(idx-1) = h_step(idxSteps(Ti));
    hidx(idx-1) = idxSteps(Ti);
    r_min(idx-1) = min(a_min(idxSteps));
    r_max(idx-1) = max(a_max(idxSteps));
  end

end
