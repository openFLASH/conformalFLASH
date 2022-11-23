%% CEFsingleSpkPerSpot
% Compute the energy filter modulator parameter in the case of a single spike of chimney per spot
%
%% Syntax
% |RF = CEFsingleSpkPerSpot(RF , k_RF , kspt , SpotPositions , t_RF , tw_RF , w_RF , w_RF0 , Beams  , Layers , Spike , nrSides , mag , apothem , BaseArea)|
%
%
%% Description
% |RF = CEFsingleSpkPerSpot(RF , k_RF , kspt , SpotPositions , t_RF , tw_RF , w_RF , w_RF0 , Beams  , Layers , Spike , nrSides , mag , apothem , BaseArea)| Description
%
%
%% Input arguments
% |RF| -_VECTOR of STRUCT_- Structure describing the shape of the Conformal Energy Filter filer
%
% |k_RF| -_SCALAR_- Index in |RF| to which the spike information is to be added
%
% |kspt| -_SCALAR_- Index of the PBS spot to which the |k_RF| spike corresponds
%
% |spotPosition| - SCALAR MATRIX_ - |spotPosition(i,:) = [x,y]| The i-th spot is located at position [x,y] in the BEV coordinates
%
% |t_RF| -_SCALAR VECTOR_- Thickness (mm) of the different steps of the spikes, corresponding to each energy layer
%
% |tw_RF| -_SCALAR VECTOR_- |tw_RF(i)| Weight of the SOBP allocated to the i-th spike
%
% |w_RF| -_SCALAR MATRIX- w_RF(i,L) weight of the hexagon/square at SpotPositions(i,:) for energy layer L. The weight is rescaled for proton of the maximum evergy
%
% |w_RF0| -_SCALAR MATRIX- w_RF(i,L) weight of the hexagon/square at SpotPositions(i,:) for energy layer L. The weight is assuming the ernergy of layer l
%
% |Beam| -_STRUCTURE_- Information about the beam
%
% |Layers| -_STRUCTURE_- Information from  |Beam.Layers| with the weight renormalised at max energy.
%
% |Spike| -_STRUCTURE_- Description of the propoerties of the CEF spikes
%
% |nrSides| - _INTEGER_ - Number of side of the polygon defining the base of the spike
%
% |mag| - _SCALAR VECTOR_ - Magnification factor [Mx,My] to project from the isocentre plane to the projection plane
%
% |apothem| -_SCALAR_- Radius (mm) of the base of the spikes
%
% |BaseArea| -_SCALAR_- Area (mm2) of the base of the spike
%
%% Output arguments
%
% |RF| -_VECTOR of STRUCT_- Structure describing the shape of the Conformal Energy Filter filer with the |k_RF| updated
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function RF = CEFsingleSpkPerSpot(RF , k_RF , kspt , SpotPositions , t_RF , tw_RF , w_RF , w_RF0 , Beams  , Layers , Spike , nrSides , mag , apothem , BaseArea)

  RF(k_RF).t_RF = t_RF; %mm Thickness corrsponding to each energy
  RF(k_RF).latticePeriod = Beams.SpotSpacing .* min(mag); %The lattice period projected into the plane of the CEF. Take the smallest magnification
  RF(k_RF).x_centre = SpotPositions(kspt,1) .* mag(1); %Position of the spike in the CEF plane
  RF(k_RF).y_centre = SpotPositions(kspt,2) .* mag(2);
  RF(k_RF).x_spot = SpotPositions(kspt,1); %Position of the PBS spot in the isocentre plane
  RF(k_RF).y_spot = SpotPositions(kspt,2);
  RF(k_RF).w = tw_RF(kspt); %The total weight of the SOBP: sum the charges at position (x,y) of the BP from all layers
  apothem_prev_layer = apothem; %radius of the circle inscribed in the base
  apothem_base = apothem;
  apothem_layer = 0;
  N_layer = length([Beams.Layers(:).Energy]);


  %Find the first spot with non zero weight.
  % Use that spot to define the thickness of the range compensator
  switch Spike.SpikeType
    case {'up','random'}
        %The centre of the spike goes up and corresponds to the shallowest Bragg peak
        ExternalLayer = min(find(w_RF(kspt,:) > 0)); %This is the index of the deepest Bragg peak in the SOBP
        direct = 1; %Layer index goes up from ExternalLayer
        InternalLayer = N_layer;
        %All deeper BP have a weight of zero. We can ignore them
    case {'down','smooth','ellipse'}
        %The centre of the spike goes down and corresponds to the deepest BP
        ExternalLayer = max(find(w_RF(kspt,:) > 0)); %This is the index of the shallowest Bragg peak in the SOBP
        direct = -1; %Layer index goes down from ExternalLayer
        InternalLayer = 1;
        %All shallower BP have a weight of zero. We can ignore them
    otherwise
      Spike.SpikeType
        error('Unknown spike type')
  end

  %Define the height of zone outside the spike / chimney
  HighestEnergyLayer = min(find(w_RF(kspt,:) > 0)); %This is the index of the deepest Bragg peak in the SOBP
  RF(k_RF).h_outside = t_RF(HighestEnergyLayer); %Height of the base of the spike

  k_count = 0;
  for k_layer = ExternalLayer : direct : InternalLayer
      %Loop for all the layers following the deepest one

      if w_RF(kspt,k_layer) ~= 0 && ~isnan(w_RF(kspt,k_layer))
          %This BP has a weight. We can create a step in the spike for this BP depth
          k_count = k_count+1;

          switch Spike.PreOptimization
            case 'Gaussian'
                % if spots grid is coarse, then apply gaussian weighting to ridge filter element surface
                apothem_layer = preOpt_gaussian(apothem_layer, apothem_prev_layer, apothem_base, Beams, nrSides, w_RF(kspt,k_layer));
            case 'uniform'
                % if spots grid is finer, then no weighting of ridge filter element surface - use "flat" surface
                apothem_layer = preOpt_uniform (apothem_layer, apothem_prev_layer, BaseArea, Beams, nrSides, w_RF(kspt,k_layer));
            otherwise
                Spike.PreOptimization
                error ('Unknown pre-optimisation algorithm')
            end %switch case

          RF(k_RF).a_max(k_count) = apothem_prev_layer;
          RF(k_RF).a_min(k_count) = apothem_layer;
          RF(k_RF).h_step(k_count) = t_RF(k_layer); %Add the height of the layer
          RF(k_RF).Energy(k_count) = Layers(k_layer).Energy; %MeV Energy at exit of the step of the CEF
          RF(k_RF).w_step(k_count) = w_RF(kspt,k_layer); %Record the weight of the spot, expressed in unit proportional to number of protons  at **maximum** energy
          RF(k_RF).w0_step(k_count) = w_RF0(kspt,k_layer); %Record the weight of the spot, expressed in unit proportional to number of protons  at **maximum** energy
          apothem_prev_layer = apothem_layer;
      else
        %There is no weight for this BP.
        %At this time, set the step width =0. However, we will still record that there is a potential step
        %at this depth so as to give one additional degree of freedom to the getOptimumSpike.m function
        k_count = k_count+1;
        RF(k_RF).a_max(k_count) = apothem_prev_layer;
        RF(k_RF).a_min(k_count) = apothem_prev_layer; %The step width is zero
        RF(k_RF).h_step(k_count) = t_RF(k_layer); %Add the range of the layer
        RF(k_RF).Energy(k_count) = Layers(k_layer).Energy; %MeV Energy at exit of the step of the CEF
        RF(k_RF).w_step(k_count) = 0; %In the IMPT plan, this layer had zero weight.
        RF(k_RF).w0_step(k_count) = 0;
      end
  end

  RF(k_RF).SpikeType = Spike.SpikeType; %Record the spike type
  %Cut the radius of the spike when reaching the centre
  w1stZero = min(find(RF(k_RF).a_max==0));
  if ~isempty(w1stZero)
    w1stZero = w1stZero -1;
    RF(k_RF).a_max = RF(k_RF).a_max(1:w1stZero);
    RF(k_RF).a_min = RF(k_RF).a_min(1:w1stZero);
    RF(k_RF).h_step = RF(k_RF).h_step(1:w1stZero);
    RF(k_RF).Energy = RF(k_RF).Energy(1:w1stZero);
    RF(k_RF).w_step = RF(k_RF).w_step(1:w1stZero);
    RF(k_RF).w0_step = RF(k_RF).w0_step(1:w1stZero);
  end

  %If smooth option is selected, compute the parameters of the curve defining the spike edges
  if (strcmp(RF(k_RF).SpikeType,'smooth'))
    fprintf('Smoothing spike edges \n')
    RF(k_RF).ShapeParam = smoothSpike(RF(k_RF).a_max , RF(k_RF).a_min , RF(k_RF).h_step); %Parameters of the function defining the edges of the spike or hole
  end

end



%-------------------------------
function integral = gaussSpotIntegral(apothem, spotSigma, nrSides)

    if apothem == 0
      integral = 0;
      return;
    end

    % define polygon
    switch nrSides
    case 6
        %Hexagon
        theta = 0:60:360;
        phi = 0; % rotation phase
        x_polygon = apothem.*cosd(theta + phi);
        y_polygon = apothem.*sind(theta + phi);

    case 4
        %square
        theta = 0:90:270;
        phi = 45; % rotation phase
        x_polygon = apothem.*cosd(theta + phi);
        y_polygon = apothem.*sind(theta + phi);
    case 1
        %Ellipse
        x_polygon = -apothem:apothem;
        %Compute the main axis of the ellispes
        [Ax , Ay] = getEllipseAxes(apothem , spotSigma.Sx ./ spotSigma.Sy);
        if (Ax < Ay)
          x_polygon = x_polygon .* spotSigma.Sx ./ spotSigma.Sy;
        end
        y_polygon = sqrt(abs(Ay.^2 .* (1 - (x_polygon./Ax).^2)));
        x_polygon = [x_polygon , flipdim(x_polygon,2)]; %Positive root
        y_polygon = [y_polygon , flipdim(-y_polygon,2)]; %negative root
    end

    %gaussFunk = @(x,y) 1/(2*pi()*spotSigma.^2) .*exp(-((x.^2 + y.^2)./(2*spotSigma.^2))); % Integrant gauss function. Normalization (x.^2 + y.^2)/(2*sigma^2) ? Amplitude -> ? 1/(2*pi()*spotSigma.^2)
    gaussFunk = @(x,y) biNorm([x , y] , 1 ,  [0,0] , spotSigma.Sx , spotSigma.Sy , spotSigma.r); % Integrant gauss function. Bi-normal function to account for elliptical spots. Area normalised to 1
    param = struct('method','gauss','tol',1e-6); % Integration method
    domain = struct('type','polygon','x',x_polygon,'y',y_polygon); % Integration domain
    integral = doubleintegral(gaussFunk, domain, param); % compute integral
end

%-------------------------------
function integral = computeGaussWeightedArea(apothem_layer, apothem_prev_layer, apothem_base, spotSigma, nrSides, spotWeight)

    % The integral having as integration limit apothem_layer (of current ridge filter level) and defined as the differentce between the normalized previous level area, the current level area and the spot weight.
    eps = 0.1; % in mm, defined as the spatial resolution of a voxel.
    if (apothem_layer < eps)
        if (gaussSpotIntegral(apothem_prev_layer, spotSigma, nrSides) > spotWeight)
            integral = gaussSpotIntegral(apothem_prev_layer, spotSigma, nrSides) ...
                       /gaussSpotIntegral(apothem_base, spotSigma, nrSides) - spotWeight;
        else
            integral = 0;
        end
    else
        integral = (gaussSpotIntegral(apothem_prev_layer, spotSigma, nrSides) - gaussSpotIntegral(apothem_layer, spotSigma, nrSides)) ...
                   /gaussSpotIntegral(apothem_base, spotSigma, nrSides) - spotWeight;
    end
end


%-------------------------------------------------
% Pre-optimise the spike Shape
% assuming a Gaussian fluence
%-------------------------------------------------
function apothem_layer = preOpt_gaussian(apothem_layer, apothem_prev_layer, apothem_base, Beams, nrSides, w_RF)

  if ((apothem_prev_layer.^2 - apothem_layer.^2) >=0)

      % define the gaussian function integral over square domain at apothem_max and variable apothem_min and subtract spotWeight
      spotGaussWeightedArea = @(apothem_layer) computeGaussWeightedArea(apothem_layer, apothem_prev_layer, apothem_base, Beams.sigmaAtCEF, nrSides, w_RF);
      options = optimset('Display','iter'); % show iterations
      x0 =[0 apothem_prev_layer]; % search interval for solution;
      eps = 0.1; % in mm, defined as the spatial resolution of a voxel.
      if (apothem_prev_layer < eps)
          root_apothem = 0;
      else
          root_apothem = fzero(spotGaussWeightedArea, x0, options); % computes the zeros of a function. Numerical version of inverse
      end

      apothem_layer = root_apothem;

  else
      %The resudial surface is too small to be manufactured
      apothem_layer = 0;
  end
end


%-------------------------------------------------
% Pre-optimise the spike Shape
% assuming an uniform fluence
%-------------------------------------------------
function apothem_layer = preOpt_uniform(apothem_layer, apothem_prev_layer, BaseArea , Beams, nrSides, w_RF)

    switch nrSides
      case 6 %HEXAGONAL
    %if (strcmp(GridLayout,'HEXAGONAL'))
        if (apothem_prev_layer.^2 - BaseArea.*w_RF./(2*sqrt(3)) >=0);
            %fprintf('case good \n')
            apothem_layer = sqrt(apothem_prev_layer.^2 - BaseArea.*w_RF./(2*sqrt(3))); % hexagonal area to draw hexagonal towers
        else
            %The resudial surface is too small to be manufactured
            apothem_layer = 0;
        end

      case 4 %SQUARE
    %elseif (strcmp(GridLayout,'SQUARE'))
        AreaPrevLayer = (2 .* apothem_prev_layer).^2; %Area of the square from previous layer. Apothem is half the base of the square
        if (AreaPrevLayer - BaseArea .* w_RF >=0)
            apothem_layer = sqrt(AreaPrevLayer - BaseArea .* w_RF)./2; % square area to draw square towers
        else
            %The resudial surface is too small to be manufactured
            apothem_layer = 0;
        end

      case 1 %Ellipse
        [Ax , Ay] = getEllipseAxes(apothem_prev_layer , Beams.sigmaAtCEF.Sx ./ Beams.sigmaAtCEF.Sy) ;
        AreaPrevLayer = pi .* Ax .* Ay; %Area of an ellipse
        if (AreaPrevLayer - BaseArea .* w_RF >=0)
            AreaThisLayer = AreaPrevLayer - BaseArea .* w_RF;
            Sx_Sy = Beams.sigmaAtCEF.Sx ./ Beams.sigmaAtCEF.Sy;
            if (Sx_Sy  > 1)
              %Sx is larger. It must be the apothem
              apothem_layer = sqrt(Sx_Sy .* AreaThisLayer ./ pi);
            else
              %Sy is larger. It must be the apothem
              apothem_layer = sqrt(AreaThisLayer ./ (pi .* Sx_Sy));
            end

        else
            %The resudial surface is too small to be manufactured
            apothem_layer = 0;
        end
    end

  end
