%% plan2weight
% Convert the MIROPT plan structure into a vector of weight that can be used by the MIROPT optimizer
% The spots weights (w) are ordered in the following manner:
% [beam1-layer1-spot1, B1L1S2, ... B1L1Sn, B1L2S1,...B1LmSn, B2L1S1,....,BwLmSn]

%
%% Syntax
% |w = plan2weight(Plan)|
%
%
%% Description
% |w = plan2weight(Plan)| Description
%
%
%% Input arguments
% |Plan| - _STRUCTURE_ - The optimised treatment plan
%  * |Plan.Beams(b).Layers(L).SpotWeights(k)| -_SCALAR VECTOR_-  Weight per fraction of the k-th spot in the layer
%  * |Plan.fractions| -_SCALAR VECTOR_- |Plan.fractions(b)| Number of fraction to deliver the b-th beam.
%
%
%% Output arguments
%
% |wF| - _SCALAR VECTOR_ - |wF(i)| optimised wieght of the i-th spot. This is the weight per fraction of each beam
%
% |wT| - _SCALAR VECTOR_ - |wF(i)| optimised wieght of the i-th spot. This is the total weight, sum of all fractions
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [wF , wT] = plan2weight(Plan)

  idx_start = 1;
  idx_end = 0;
  wF = [];
  wT = [];

  %Construct the vector of BP weights from the content of the plan
  for b = 1: size(Plan.Beams,2)
      for l = 1:size(Plan.Beams(b).Layers,2)
          idx_end = idx_end + size(Plan.Beams(b).Layers(l).SpotWeights , 2);
          wF = [wF , Plan.Beams(b).Layers(l).SpotWeights]; %This is the weight PER FRACTION. In SpotWeightsOptimization at line 85, the weight are divided by the number of fractions.
          wT = [wT , Plan.Beams(b).Layers(l).SpotWeights .* Plan.fractions(b)]; %This is the total weight for all fractions
          %w(idx_start:idx_end) = Plan.Beams(b).Layers(l).SpotWeights;
          idx_start = idx_end + 1;
      end
  end

  wF = double(wF); %cast to double
  wT = double(wT); %cast to double

end
