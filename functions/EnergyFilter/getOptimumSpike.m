%% getOptimumSpike
% Optimise the shape of one spike of the CEF
% The optimisation minimizes the difference between the fluence map at each energy between
%  * a beam going through the spike
%  * an IMPT plan comprising a beam composed of several Bragg peaks with the required weights
% Both beams then go through a 150mm water phantom placed at isocentre.
% The fluence is compared 50mm downstream to the water phantom
%
% The oiptimisation varies:
%   * width of each steps of the spike (a_max, a_min) are varied
%   * multiplies the weight of the SOBP by a factor A to rescale the fluence map
%
%% Syntax
% |[Plan , FluenceRef] = getOptimumSpike(Plan , b)|
%
%
%% Description
% |[Plan , FluenceRef] = getOptimumSpike(Plan , b)| Description
%
%
%% Input arguments
% |Plan| - _struct_ - Information about the initial guess of the spike shape:
%   * |Plan.showGraph| -_LOGICAL_- If true, display the graphs ofthe fluence at different steps i nthe computations
%   * |Plan.Spike.WET| -_SCALAR_- Relative water equivalent thickness of the Conformal Energy Filter material
%   * |Plan.Spike.intrpCTpxlSize| -_SCALAR_- Lateral spatial resolution (mm) of the 3D printer
%   * |Plan.Spike.min_thickness| -_SCALAR_- Thickness (mm) in of the base on which the spikes are built. The base has the same |R_WET| as the spikes
%   * |Plan.Beams(b)| -_STRUCTURE_- Information about the b-th beam
%       * |Plan.Beams(b).spotSigmaAtSkin| -_SCALAR_- Spot lateral sigma (mm) at skin level
%       * |Plan.Beams(b).spotSigma| -_SCALAR_- Spot lateral sigma (mm) at maximum of deepest Bragg peak along the optical axis
%       * |Plan.Beams(b).Layers(L).Energy| -_SCALAR_- Energy (MeV) of the L-th layer
%       * |Plan.Beams(b).GridLayout| -_STRING_- Layout of the PBS spots on the grid. Options: HEXAGONAL (default), SQUARE
%       * |Plan.Beams(b).RidgeFilter| -_VECTOR of STRUCT_- Structure describing the shape of the Conformal Energy Filter filer
%          * |Plan.Beams(b).RidgeFilter(k).x_centre| -_SCALAR_- X coordinate (mm in IEC Gantry CS) of central axis of the k-th spike
%          * |Plan.Beams(b).RidgeFilter(k).y_centre| -_SCALAR_- Y coordinate (mm in IEC Gantry CS) of central axis of the k-th spike
%          * |Plan.Beams(b).RidgeFilter(k).apothem| -_SCALAR_- Apothem (mm) of the base of the spike
%          * |Plan.Beams(b).RidgeFilter(k).w| -_SCALAR_- Weight of the PBS spot of energy |max(Layers(:).Energy)|  to deliver at the base of spike = weight of the SOBP
%          * |Plan.Beams(b).RidgeFilter(k).h_outside| -_SCALAR_- Height (mm) of the base of the k-th spike
%          * |Plan.Beams(b).RidgeFilter(k).a_max| -_SCALAR VECTOR_- |RF(k).a_max(L)| Max apothem (mm) of the L-th polygon (apothem of the L-th hexagon/square)
%          * |Plan.Beams(b).RidgeFilter(k).a_min| -_SCALAR VECTOR_- |RF(k).a_min(L)| Min apothem (mm) of the L-th polygon (apothem of the (L+1)-th hexagon/square)
%          * |Plan.Beams(b).RidgeFilter(k).h_step| -_SCALAR VECTOR_- |RF(k).h_step(L)| Height (mm) of the L-th polygon from the base
%          * |Plan.Beams(b).RidgeFilter(k).w_step = w_RF(kspt,k_layer)| Weight of the step k_layer at location k_RF, expressed in unit proportional to number of protons
%
% |b| -_INTEGER_- Index of the beam for which the spike is to be optimised
%
%
%% Output arguments
%
% |Plan| - _STRUCT_ - Updated treatment plan with the width of the spikes optimised:
%
%
% REFERENCES
% [1] PBS 10s Concept Performance Assessment by Jonathan Hubeau MID 29367
% [2] https://nl.mathworks.com/help/optim/ug/fmincon.html
%
%
%% Contributors
% Authors : Rudi Labarbe (open.reggui@gmail.com)

function Plan = getOptimumSpike(Plan , b)

  Meas_Zg = 0; %Place water phantom at isoncetre
  Iso2RS = Plan.Beams(b).RSinfo.IsocenterToRangeShifterDistance; %Distance from the top of the CEF spike to the isocentre

  %Compute fluence map with the CEF
  %================================
  fprintf('Computing the fluence maps with initial CEF \n')
  pencil_x = [Plan.Beams(b).RidgeFilter(:).x_centre]; %Position of the centre of the spikes is also the centre of the PBS spots
  pencil_y = [Plan.Beams(b).RidgeFilter(:).y_centre];
  weight_table = [Plan.Beams(b).RidgeFilter(:).w]; % w(spt) is the weight of the spt-th PBS spot. The value is proportional to the number of proton at **maximum** energy

  % Get the model parameters by calling fluenceWithCEF
  [Fluence , X_flu , Y_flu , E_flu, sigmas , CEF_thickness] = fluenceWithCEF(Plan , b , pencil_x , pencil_y  , weight_table , [] ,  Meas_Zg , Plan.showGraph , [], 0, 'config_RS_CEM');

  if (isfield(Plan.Beams(b), 'BlockData'))
      %Get the information about the aperture block
      %The GOF function shall only take into account the fluence inside the aperture block
      for cntIdx = 1:numel(Plan.Beams(b).BlockData)
        aperture = Plan.Beams(b).BlockData{cntIdx};
        pxlSize = [min(diff(X_flu)) ,  min(diff(Y_flu))]; %Size of the pixel in fluence map
        orig = [min(X_flu) , min(Y_flu)]; %coordinate of the first pixel of the fluence map
        xi = double(round(1+(aperture(:,1)-orig(1)) ./ pxlSize(1))); %Pixel indices of the contour of the pareture
        yi = double(round(1+(aperture(:,2)-orig(2)) ./ pxlSize(1)));
        apertureMask = poly2mask(xi,yi,double(numel(X_flu)),double(numel(Y_flu))); %mask defining the paerture in the fluence image
      end
  else
      %There is no aperture block
      apertureMask = [];
  end

  if Plan.showGraph
        figure(50)
        [yg , xg] = meshgrid(Y_flu , X_flu); %meshgrid puts Y to the first index of xg and yg
        [~,imx] = max(Fluence,[],'all','linear');
        [~, ~ , iZ] = ind2sub(size(Fluence),imx);
        Fcef = squeeze(Fluence(:,:,iZ));
        Fcef = Fcef ./ max(Fcef,[],'all');
        contour(xg' , yg' ,Fcef' ,'ShowText','on')
        if(isfield(Plan.Beams(b), 'BlockData'))
          hold on
          plot(aperture(:,1),aperture(:,2),'-r') %Overlay aperture shape on the contour
        end
        xlabel('X (mm)')
        ylabel('Y (mm)')
        title(['Fluence with unoptimised CEM at ' num2str(E_flu(iZ)) 'MeV'])
        grid on
  end


  %Iterative optimisation of the fluence
  %This is done one spike at a time
  %=====================================
  [~ , step , nrSides] = getConvGridParam(Plan , b);

  Xglobal = [];

  for SpikeIdx = 1:numel(Plan.Beams(b).RidgeFilter)

      nElayers = numel(unique([Plan.Beams(b).RidgeFilter(SpikeIdx).Energy]));
      if (nElayers > 1)
            %The optimisation is done only if there is more than one step in the spike

            %Compute the nominal fluence map around the selected spike
            fprintf('Computing the reference fluence maps for spike %d at (%3.2f , %3.2f) mm \n',SpikeIdx , Plan.Beams(b).RidgeFilter(SpikeIdx).x_centre , Plan.Beams(b).RidgeFilter(SpikeIdx).y_centre)
            [X_flu , Y_flu ] =  makeConvGrid(Plan.Beams(b).spotSigma , step , Plan.Beams(b).RidgeFilter(SpikeIdx).x_centre , Plan.Beams(b).RidgeFilter(SpikeIdx).y_centre); %REcompute the X and Y convolution grid around this spike
            %This is the new cartesian grid around the spike. The convolution will take place on that grid only
            FluenceRef = computeFluenceMap(X_flu , Y_flu , E_flu , Plan.Beams(b) , [] ); %The fuence map is computed at the depth of the max BP
            FluenceRef = FluenceRef ./  max(FluenceRef,[],'all'); %Normalise the reference fluence

            %Define inital guess
            %====================
            switch Plan.Beams(b).RidgeFilter(SpikeIdx).SpikeType
              case 'smooth'
                x0 = Plan.Beams(b).RidgeFilter(SpikeIdx).ShapeParam;

              case 'fractal'
                x1 = [ reshape(Plan.Beams(b).RidgeFilter(SpikeIdx).h_grid,1,numel(Plan.Beams(b).RidgeFilter(SpikeIdx).h_grid))]; %Initial guess.
                allHeights = Plan.Beams(b).RidgeFilter(SpikeIdx).t_RF; %VEctor with all the steps heights available for all cells
                for i = 1:numel(x1)
                  x0(i) = find(x1(i)==allHeights); %x0 is the index of the possible energy layers
                end

              case 'multiStairs'
                x0 = [Plan.Beams(b).RidgeFilter(SpikeIdx).h_grid];

              otherwise
                r0 = abs(diff([Plan.Beams(b).RidgeFilter(SpikeIdx).a_max 0])); %Radial distance of each jump to the next step
                x0 = r0 ./ sum(r0,'all'); %Initial guess. This is the radial distance normalised to 1
                %One additional degree of freedom is the height of the base outside of the spike.
                %This will be the 1st element of the vector.
                [~ , idxMinHeight] = min(Plan.Beams(b).RidgeFilter(SpikeIdx).h_step);
                x0 = [Plan.Beams(b).RidgeFilter(SpikeIdx).apothem , idxMinHeight , x0]; %The index in the vector |RidgeFilter(SpikeIdx).h_step(x0(2))| defining the height the base is a fitting parameter

            end

            %Define the options for the optimizer
            switch Plan.Beams(b).RidgeFilter(SpikeIdx).SpikeType
              case 'fractal'
                  %The Xs are indices of the height vector. They must be integer. Setp size should be close to integer value
                  % [2] 	 FinDiffRelStep : Vector step size factor for finite differences. The forward finite differences delta are
                  %             delta = v.*sign′(x).*max(abs(x),TypicalX);
                  TypX = ones(size(x0)) .* numel(Plan.Beams(b).RidgeFilter(SpikeIdx).t_RF) ./ 2; %Typical value of X
                  options = optimset('Display','iter-detailed','TolX',1e-6,'TolFun', 1e-6, 'FinDiffRelStep' , ones(size(x0)) ./ (numel(allHeights)) , 'TypicalX' , TypX);
                  if Plan.showGraph
                    options = optimset(options , 'PlotFcns', {@optimplotfval});
                  end

              otherwise
                    if (isfield(Plan.Spike,'OptMaxIter') && Plan.Spike.OptMaxIter ~= 0)
                        options = optimset('Display','iter','MaxIter', Plan.Spike.OptMaxIter);
                    else
                        options = optimset('Display','iter','TolX',1e-2,'TolFun', 1e-4);
                    end

                    if Plan.showGraph
                      options = optimset(options , 'PlotFcns', {@plotTheGraphs,@optimplotfval});
                    end
              end

            %Optimise spike shape
            %=======================
            %gofType = 'SOBP'
            gofType = 'fluence'

            switch Plan.Beams(b).RidgeFilter(SpikeIdx).SpikeType
              case 'fractal'
                  %Constrained optimisation
                  A = [];
                  bs = [];
                  Aeq = [];
                  beq = [];
                  nonlcon = [];
                  HeightMax = max(find(Plan.Beams(b).RidgeFilter(SpikeIdx).w_step)); %There is no tower higher than this value in the IMPT plan. Constrain the optimiser to stay below this limit
                  lb = zeros(size(x0)) + 0.6;  %Lower bound
                  ub = ones(size(x0)) .* HeightMax + 0.4; %Upper bound
                  %ub = ones(size(x0)) .* numel(allHeights) + 0.4; %Upper bound
                    %NB: contrarily to what is suggested on the MAtlab documentation, the boundary conditions are defined as STRICTLY larger or STRICTLY smaller
                    % So we defined the boundary as ]0.6 , numel(allHeights)+0.4] so that the rounding will bring back to the correct index
                  [x,fval,exitflag] = fmincon(@gofSelector , x0 ,A,bs,Aeq,beq,lb,ub,nonlcon , options,  gofType , Plan , b , SpikeIdx, Meas_Zg , sigmas , FluenceRef , [] , Iso2RS , E_flu , Plan.Beams(b).RidgeFilter(SpikeIdx).SpikeType);

              otherwise
                  %Unconstrained optimisation
                  [x,fval,exitflag] = fminsearch(@gofSelector, x0, options, gofType , Plan , b , SpikeIdx, Meas_Zg , sigmas , FluenceRef , [] , Iso2RS , E_flu , Plan.Beams(b).RidgeFilter(SpikeIdx).SpikeType);
            end


            %Update the width of the spikes for the exported plan
            Plan.Beams(b).RidgeFilter(SpikeIdx) = packRgofFluence(x , Plan.Beams(b).RidgeFilter(SpikeIdx), Plan.Spike.intrpCTpxlSize);

        end

 end %for SpikeIdx

end %function


%===============================
% Plot the evolution of the maximum diameter and the height of the external area as a function of
%
% INPUT
% |varargin| -_CELL VECTOR_- Same parameters than in gofSelector.m
%===============================
function stop = plotTheGraphs(x,optimvalues,state,varargin)

      SpikeType = varargin{11};
      stop = false;

      switch SpikeType
       case {'fractal' , 'multiStairs'}
           %do nothing

         otherwise
             plot(optimvalues.iteration,x(1),'ob')
             hold on
             plot(optimvalues.iteration,x(2),'+r')
             legend({'Max radius (mm)','External height index'})
      end

      xlabel('Iteration')
      grid on
   end
