%% packRgofFluence
% Extract from the optimised parameter vector x:
%  * the radius into a_min and a_max of the steps of the spike
%  * the multiplicative factor of the fluence
%
%% Used to compute the gofFluence.m function
%
%% Syntax
% |[RidgeFilter , A] = packRgofFluence(x , RidgeFilter)|
%
%
%% Description
% |[RidgeFilter , A] = packRgofFluence(x , RidgeFilter)| Description
%
%
%% Input arguments
% |x| -_SCALAR VECTOR_- The optimised parameters |x = [maxR , h , dr1 , dr2, ... , drn]|
%                 where the relative width of the lowers step is dr1 and the relative width of the highest step is drn
%                 h is the index in the vector |RidgeFilter(SpikeIdx).h_step(x0(2))| defining the height the base is a fitting parameter
%                 maxR is the maximum apothem of the spike
%
% |RidgeFilter| -_VECTOR of STRUCT_- Structure describing the shape of the Conformal Energy Filter filer
%     * |Plan.Beams(b).RidgeFilter(k).a_max| -_SCALAR VECTOR_- |a_max(L)| Largest apothem (mm) of the L-th step
%     * |Plan.Beams(b).RidgeFilter(k).a_min| -_SCALAR VECTOR_- |a_min(L)| Smalest apotherm of the L-th step
%     * |Plan.Beams(b).RidgeFilter(k).h_step| -_SCALAR VECTOR_- |RF(k).h_step(L)| Height (mm) of the L-th polygon from the base
%
% |pxlSize| -_SCALAR_- [OPTIONAL. Default : 0] Minimum step size (mm) of the width of a step. Corresponds to the lateral resolution of the 3D printer or the voxel resolution of CT scan
%
%% Output arguments
% |RidgeFilter| -_VECTOR of STRUCT_- Updated structure describing the shape of the Conformal Energy Filter filer
%     * |Plan.Beams(b).RidgeFilter(k).a_max| -_SCALAR VECTOR_- |a_max(L)| Largest apothem (mm) of the L-th step
%     * |Plan.Beams(b).RidgeFilter(k).a_min| -_SCALAR VECTOR_- |a_min(L)| Smalest apotherm of the L-th step
%     * |Plan.Beams(b).RidgeFilter(k).h_outside| -_SCALAR_- Height (mm) of the base of the k-th spike
%
% |Xout|  -_SCALAR VECTOR_- Modified values of of |x| so as to meet the constraints.
%
% |Xother| -_SCALAR VECTOR_- component of |x| that are not included in |RidgeFilter|
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [RidgeFilter , Xout , Xother]  = packRgofFluence(x , RidgeFilter, pxlSize)

  if nargin < 3
    pxlSize =0;
  end

  Xother = [];
  NbSpikes = numel(RidgeFilter); %Number of spikes in the CEF

  switch RidgeFilter(1).SpikeType %All spike have the same shape. We can switch based on type of 1st spike
  case 'smooth'
        %The spike has a smooth shape defined by an analytical expression
        for k = 1:NbSpikes
            RidgeFilter(k).ShapeParam = x;
            % r = 0:0.1:RidgeFilter(k).apothem;
            % f = spikeShapeFnc(x,r);
            % a_cent = interp1(f,r,RidgeFilter(k).h_step); %Coordinate of centre of step of height h_step
            a_cent = spikeShapeH2R(RidgeFilter(k).ShapeParam, RidgeFilter(k).h_step, RidgeFilter(k).apothem);
            Nstp = numel(RidgeFilter(k).h_step);

            figure(600)
            hold on
            plot(a_cent,RidgeFilter(k).h_step,'+r')

            % |a_cent| gives the coordinate of the centre of the steps. Unfortunately, not all sequence of position can
            % be converted into |a_min| and |a_max| such that (|a_min| + |a_max|) / 2 = |a_cent|
            % Some sequence of |a_cent| give position that cannot be located at the centre of rings for ALL steps.
            % We will therefore give up the idea of building rings with |a_cent| forced to be at the middle of the ring
            % Instead, the divide between 2 rings will be chosen at the middle point between successive values in |a_cent|
            seq1 = [a_cent(2:end)];
            seq2 = [a_cent(1:end-1)];
            middlePts = (seq1 + seq2) ./2;
            RidgeFilter(k).a_min = [middlePts , 0];
            RidgeFilter(k).a_max = [RidgeFilter(k).apothem , middlePts];
        end

  case 'fractal'
      allHeights = RidgeFilter(1).t_RF'; %Vector with all the step heights available for all cells
      x = allHeights(round(x));
      Xout = round(x);

      %TODO This will not work when we optimise multiple cell at the same time
      for k = 1:NbSpikes
          RidgeFilter(k).h_grid = reshape(x, RidgeFilter(k).NbColumns , RidgeFilter(k).NbColumns);
      end

  case 'multiStairs'
      x = abs(x); %CONSTRAINT: We do not want negative weight
      x = x ./ sum(x);
      Xout = x;
      %TODO This will not work when we optimise multiple cell at the same time
      for k = 1:NbSpikes
          RidgeFilter(k).h_grid = x; %normalise the weights
      end

  otherwise
        %The spike is made of a sequence of steps
        x = abs(x); %CONSTRAINT: We do not want negative radius or intensity
        Xout = x;
        idx0 = 1; %x(idx0) is the first step of the k-th spike

        for k = 1:NbSpikes

          %FIRST ELEMENT OF VECTOR is the maximum diameter of the chimney
          %Compute the maximum apothem of the spike.
          maxR = x(idx0);
          if (maxR > RidgeFilter(k).latticePeriod./2)
            maxR = RidgeFilter(k).latticePeriod./2; %CONSTRAINT: clip the radius to half the distance between spikes
            Xout(idx0) = maxR;
          end
          idx0 = idx0 + 1;

          %SECOND ELEMENT OF VECTOR is the outside height
          %clip the height of the outside are to the available step heights
          StpIdx = round(x(idx0));
          if (StpIdx <= 1)
              StpIdx = 1; %CONSTRAINT: The index cannot be smaler than 1
              Xout(idx0)= 1;
          end
          if (StpIdx > numel(RidgeFilter(k).h_step))
              StpIdx = numel(RidgeFilter(k).h_step); %CONSTRAINT: The index cannot be larger than the length of the vector
              Xout(idx0) = StpIdx;
          end
          RidgeFilter(k).h_outside = RidgeFilter(k).h_step(StpIdx); %Define the height of the base outside the spike
          idx0 = idx0 + 1;

          NbLayers = numel(RidgeFilter(k).a_max); %Number of energy layers
          r = dR2R(x(idx0:idx0+NbLayers-1) , maxR); %The lateral coordinate of each step
          if pxlSize > 0
            r = rounding(r,pxlSize);
          end

          RidgeFilter(k).a_min = [r(2:end),0];
          RidgeFilter(k).a_max = r;
          idx0 = idx0 + NbLayers; %Index of the first index (of X) for next spike
        end

    end %switch case
end



%=============================
% Convert the relative width of each steps into the smallest and largest apothem for each step
%
% INPUT
% |dr| -_SCALAR VECTOR_- |dr = [dr1 , dr2, ... , drn]| the relative width of the lower step is dr1 and the relative width of the highest step is drn
%               The width of the i-th step is dr(i) .* Rtot
% |Rtot| -_SCALAR_- Total width of the base of a spike
%=============================
function r = dR2R(dr , Rtot)
  dr = abs(dr); %only positive distance
  dr = dr ./ sum(dr,'all'); %Normalise the delta so that each width is expressed as a fraction of the base width
  r = Rtot .* flipdim(cumsum(flipdim(dr,2)),2); %Determine the proportion of each steps
end

%=============================
% Round the value of X to the closest available step height
%
% INPUT
% |x| -_SCALAR VECTOR_- The vector of values to be rounded
% |allHeights| -_SCALAR VECTOR_- List of all availbe heights
%
% OUTPUT
% |x| -_SCALAR VECTOR_- The value of the initial |x| are rounded to the closest value listed in |allHeights|
%=============================
function x = roundX2stepsHeights(x , allHeights)
  [t_RF , xr] = meshgrid(allHeights , x);
  [~ , roundIdx] = min((t_RF - xr).^2,[],2);
  x = allHeights(roundIdx);
  x = x';

end
