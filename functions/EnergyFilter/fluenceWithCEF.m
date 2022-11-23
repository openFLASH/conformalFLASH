%% fluenceWithCEF
% Semi-Analytical Transfer Function (SATF)
% Compute analytically the proton fluence downstream of a Conformal Energy Filter (CEF)
% This transfer function calculates the beam fluence for every energy exiting the CEF in any transversal plane downstream to the filter at |distance|.
% The fluence map is computed for different proton energies, corrsponding to the different step height of the CEM.
% The fluence map are sorted in decreasing nergy order.
%
% Optionally, a water tank can be placed at the measurement plane. The fluence is then computed 50mm downstream from the water WaterTankThckness
% The water tank simulate the additional scattering of the beam introduced in the patient.
%
%% Syntax
% |[Fluence , Xmeas , Ymeas , E_flu , sigmas , CEF_thickness , maskLayer] = fluenceWithCEF(Plan , b , pencil_x , pencil_y  , SAD , Meas_Zg , showGraph , sigmas , dilateFlag)|
%
%
%% Description
% |[Fluence , Xmeas , Ymeas , E_flu , sigmas , CEF_thickness , maskLayer] = fluenceWithCEF(Plan , b , pencil_x , pencil_y  , SAD , Meas_Zg , showGraph , sigmas , dilateFlag)| Description
%
%
%% Input arguments
%
% |Plan| - _struct_ - MIROpt structure with updated information:
%   * |Plan.Spike.WET| -_SCALAR_- Relative water equivalent thickness of the Conformal Energy Filter material
%   * |Plan.Spike.intrpCTpxlSize| -_SCALAR_- Lateral spatial resolution (mm) of the 3D printer
%   * |Plan.Spike.min_thickness| -_SCALAR_- Thickness (mm) in of the base on which the spikes are built. The base has the same |R_WET| as the spikes
%   * |Plan.Spike.SpikeType| -_STRING_- Type of spike to be designed. The centre of the spike corresponds to the BP with smaller range ('up') or the largest range ('down') or randomise pixel column ('random'), apply Gaussian filter to "smear"('smooth'), draw elliptical spike ('ellipse')
%   * |Plan.Beams(b).Layers(L).Energy| -_SCALAR_- Energy (MeV) of the L-th layer
%   * |Plan.Beams(b).GridLayout| -_STRING_- Layout of the PBS spots on the grid. Options: HEXAGONAL (default), SQUARE
%   * |Plan.Beams(b).RidgeFilter| -_VECTOR of STRUCT_- Structure describing the shape of the Conformal Energy Filter filer
%     * |Plan.Beams(b).RidgeFilter(k).x_centre| -_SCALAR_- X coordinate (mm in IEC Gantry CS) of central axis of the k-th spike
%     * |Plan.Beams(b).RidgeFilter(k).y_centre| -_SCALAR_- Y coordinate (mm in IEC Gantry CS) of central axis of the k-th spike
%     * |Plan.Beams(b).RidgeFilter(k).w| -_SCALAR_- Weight of the PBS spot of **maximum** energy to deliver at the base of spike = weight of the SOBP. Proportional to the number of protons of maximum energy to deliver
%     * |Plan.Beams(b).RidgeFilter(k).w_step(k_layer)| Weight of the step k_layer at location k_RF, expressed in unit proportional to number of protons at **maximum** energy
%     * |Plan.Beams(b).RidgeFilter(k).h_step| -_SCALAR VECTOR_- |RF(k).h_step(L)| Height (mm) of the L-th ring
%     * |Plan.Beams(b).RidgeFilter(k).a_max| -_SCALAR VECTOR_- |a_max(L)| Largest apothem (mm) of the L-th step
%     * |Plan.Beams(b).RidgeFilter(k).a_min| -_SCALAR VECTOR_- |a_min(L)| Smalest apotherm of the L-th step
%     * |Plan.Beams(b).RidgeFilter(k).Energy| -_SCALAR VECTOR_- |RF(k).Energy(L)| Energy (MeV) of the protons coming out of the L-th ring from the base
%
% |b| -_SCALAR_- Index of the beam in |Plan.Beams(b)| for which the fluence is to be computed
%
% |pencil_x|  -_SCALAR VECTOR_- |pencil_x(spt)| X coordinate (mm in IEC Gantry CS) of central axis of the spt-th PBS spot
%
% |pencil_y|  -_SCALAR VECTOR_- |pencil_y(spt)| Y coordinate (mm in IEC Gantry CS) of central axis of the spt-th PBS spot
%
% |weight_table| -_SCALAR VECTOR_- Weight of the i-th spot
%
% |SAD| -_SCALAR VECTOR_- |SAD = [x , y]| Proton source to axis distance (mm) for scanning along the X and Y axis of the IEC-gantry CS
%                         Used to rescale Xmeas and Ymeas to take into account the beam divergence.
%                         If SAD=[], then the beam divergence is ignored.
%
% |Meas_Zg| -_SCALAR_- Z coordinate in IEC-Gantry (mm) of the plane in which the fluence is measured
%
% |showGraph| -_INTEGER_- [Optional: defauly = 0] If 1, display the graphs ofthe fluence at different steps i nthe computations
%
% |sigmas| -_STRUCTURE_- [OTPIONAL: if absent, the sigma are computed by the function]. sigma (mm) of the lateral spread of the beam at |distance| of the scatterer
%       * |sigmas.CEF| -_SCALAR_- Lateral spread introduced by the CEF
%       * |sigmas.water| -_SCALAR_- Lateral spread introduced by the water
%
%
%% Output arguments
%
% |Fluence(x,y,E)| -_SCALAR MATRIX_-  Fluence at |distance| downstream from the CEF. |Fluence(x,y,E)|  fluence at position (x,y) for proton of energy |Plan.Beams(b).Layers(E).Energy|
%
% |Xmeas| -_SCALAR VECTOR_- |Xmeas(x)| Cartesian X coordinate (mm) of pixel (x,y) in the fluence map
%
% |Ymeas| -_SCALAR VECTOR_- |Ymeas(y)| Cartesian Y coordinate (mm) of pixel (x,y)  in the fluence map
%
% |E_flu| -_SCALAR VECTOR_- Energy (MeV) for the step(E). The energies are sorted in decreasing order
%
% |sigmas| -_STRUCTURE_- [OTPIONAL: if absent, the sigma are computed by the function]. sigma (mm) of the lateral spread of the beam at |distance| of the scatterer
%       * |sigmas.CEF| -_SCALAR_- Lateral spread introduced by the CEF
%       * |sigmas.water| -_SCALAR_- Lateral spread introduced by the water
%
% |CEF_thickness| -_SCALAR VECTOR_- |CEF_thickness(E)| Thickness (mm) of the l-th step; i.e. the E-th step generating |Fluence(x,y,E)|.
%
% |maskLayer(x,y,E)| -_SCALAR MATRIX_-  Mask defining the shape of the E-th step of the CEF.
%
% REFERENCES
% [1] PBS 10s Concept Performance Assessment by Jonathan Hubeau MID 29367
%
%% Contributors
% Authors : Frauke Roellinghoff , Jonathan Hubeau , Rudi Labarbe, Lucian Hotoiu (open.reggui@gmail.com)

function [Fluence , Xmeas , Ymeas , E_flu , sigmas , CEF_thickness , maskLayer] = fluenceWithCEF(Plan , b , pencil_x , pencil_y , weight_table , SAD , Meas_Zg , showGraph , sigmas , dilateFlag, configObjectsInBeam)

if nargin < 11
    configObjectsInBeam = 'config_CEM_RS';
end

if nargin < 10
dilateFlag = 0;
end

if nargin < 8
  showGraph = 0;
end
if nargin < 9
  sigmas.water = [];
  sigmas.sigmaCEF = [];
end

if isempty(sigmas)
  sigmas.water = [];
  sigmas.sigmaCEF = [];
end

sigmaAtCEF = Plan.Beams(b).sigmaAtCEF; %mm PBS spot before range shifter and CEF
CEF_Zg = Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance;


[~ , pxlSize , nrSides , GridLayout ] = getConvGridParam(Plan , b);
Emax = max([Plan.Beams(b).Layers(:).Energy],[],'all'); %The maximum energy of incoming protons

x_spike = [Plan.Beams(b).RidgeFilter(:).x_centre]; %Position of the centre of the spikes
y_spike = [Plan.Beams(b).RidgeFilter(:).y_centre];
Nb_Spkies = numel(x_spike); %Number of spikes
CEF_x = unique(x_spike); %The coordinate of the centre of the CEF spikes in the grid
CEF_y = unique(y_spike);

%The thicknesses are ordered in ascending order
CEF_thickness  = sort(unique([Plan.Beams(b).RidgeFilter(:).h_step]),'ascend'); %The thickness of all steps of all the spikes. Steps are ordered from lowest to highest
E_flu  = sort(unique([Plan.Beams(b).RidgeFilter(:).Energy]),'descend'); %The energy of the proton emerging from the steps. The energy must be reordered to match the step size
N_thickn = length(E_flu); %Number of different energy layers

%make a cartesian grid around a spike
[X_flu , Y_flu , n_grid_cart , pts_Spike , SpikeSize] =  makeConvGrid(Plan.Beams(b).spotSigma , pxlSize , CEF_x , CEF_y);
[Ymap , Xmap] = meshgrid(Y_flu,X_flu); %2D coordinates maps for the fluence. Xmap(i,j) : i along the X axis, Y along the Y axis
Fluence     = zeros(numel(X_flu),numel(Y_flu), N_thickn);
maskLayer  = zeros(numel(X_flu),numel(Y_flu), N_thickn);

%Compute incident fluence on each spike
%--------------------------------------
if (showGraph == 1)
  fprintf('Computing fluence at entrance of CEF \n')
end

thick = [];
for SpikeIdx = 1:Nb_Spkies
    %loop over spikes

    %Position of centre of the neighbour spots in CS centered on i-th spike
    Xn = x_spike(SpikeIdx) - pencil_x;
    Yn = y_spike(SpikeIdx) - pencil_y;

    %-------------------------------------------------------
    % Compute the weight contribution of neighbourghing spots
    % and take into account the gaussian lateral profile of spots
    %-------------------------------------------------------
    Fbeam = getWofNeighhbourgh(Xn , Yn , weight_table , sigmaAtCEF , pts_Spike , CEF_Zg); %|Fbeam(i)| is the fluence map at spike (SpikeIdx) at coordiante pts_Spike(i,:)=[x,y]
    Fbeam = reshape(Fbeam,SpikeSize);
    thick = [thick , Plan.Beams(b).RidgeFilter(SpikeIdx).h_step(:)']; %Used for display

    %get spike location
    [~ , idx] = min((Xmap - x_spike(SpikeIdx)).^2 + (Ymap - y_spike(SpikeIdx)).^2,[],'all','linear'); %Find the pixel closest to the centre of spike
    [fcx,fcy] = ind2sub(size(Xmap), idx(1)); %the indices of the centre of the spike in coordinate vectors X_flu and Y_flu
    fcx = fcx - round(n_grid_cart/2); %Index of the corner of the spot in the fluence map
    fcy = fcy - round(n_grid_cart/2);

    %Compute the elevation map of the spike
    %Use the same function than the one used in makeCEF so that the fluence map is computed on exactly the same shape
    switch GridLayout
      case 'HEXAGONAL'
        nrSidesLattice = 6;
      case 'SQUARE'
        nrSidesLattice = 4;
    end
    BaseSize = struct('baseHeight', [], 'apothem', []); %Dummy variable
    if isfield(Plan.Beams , 'RangeCompensator')
      %There is a range compensator. Add its shape to the spikes
      %As it may change the height of the spike and hence the weight attributed to each Bragg peak
      %it should be included in the spike design in order to allocate the proper pixels to the proper mask
      RangeCompensator = Plan.Beams.RangeCompensator;
      RangeCompensator.BDL = Plan.BDL;
    else
      RangeCompensator = [];
    end

    [MaskSpk , ~ , RcThickness] = riseOneSpike(Plan , b, SpikeIdx , pts_Spike' , SpikeSize , CEF_thickness ,  RangeCompensator);

    if isfield(Plan.Beams , 'RangeCompensator')
      %REmove the range compensator from the spike height
      minRCthick = min(RcThickness(RcThickness>0),[],'all');
      if isempty(minRCthick)
        %There are no pixel strictly above level 0
        minRCthick = 0;
      end
      MaskSpk = MaskSpk - minRCthick;
    end

    if dilateFlag
      %Make a lateral dilation
      MaskSpk = dilateSpike(MaskSpk);
    end

    %Compute the fluence map for every step of the spike
    if size(CEF_thickness,2) > 1
      %Reshape the vector in the correct direction
      CEF_thickness = CEF_thickness';
    end
    Boundaries = ([0 ; CEF_thickness ] + [CEF_thickness ; CEF_thickness(end).*10]) ./2;
          %Because of the range compensator, there can be heights between two levels
          %The continuous height of the range compensator will be discritised into the discrete energy masks
          % Boundaries defines the limit for this discrtisation.

    for steps = 1:numel(CEF_thickness)
        mask = (MaskSpk >= Boundaries(steps)) .* (MaskSpk < Boundaries(steps+1));
        Fring = mask .* Fbeam; %multiply by the fluence of the PBS spot

        % The masked fluence is copied to the complete fluence map at the correct energy layer AND at the correct X,Y position
        Fluence   (fcx+(1:n_grid_cart) , fcy+(1:n_grid_cart) , steps) = Fluence   (fcx+(1:n_grid_cart) , fcy+(1:n_grid_cart) , steps) + Fring;
        maskLayer (fcx+(1:n_grid_cart) , fcy+(1:n_grid_cart) , steps) = maskLayer(fcx+(1:n_grid_cart) , fcy+(1:n_grid_cart) , steps) + mask;  %record the shape of the mask for this layer
    end
end

%Display the step profile
if (showGraph == 2 && N_thickn > 1)

    figure(150)
    X = pts_Spike(:,1);
    Y = pts_Spike(:,2);
    X = reshape(X,SpikeSize);
    Y = reshape(Y,SpikeSize);
    %contour(X,Y,MaskSpk,unique(thick),'ShowText','on')
    maxM = max(MaskSpk,[],'all');
    minM = min(MaskSpk,[],'all');
    image((MaskSpk-minM) .* 255 ./ (maxM-minM) )
    % Rmax = round(1.3 .* Plan.Beams(b).RidgeFilter(SpikeIdx).latticePeriod./2);
    % axis([-Rmax , Rmax , -Rmax , Rmax]);

    grid on
    % xlabel('Xg (mm)')
    % ylabel('Yg (mm)')
    drawnow

end





switch configObjectsInBeam
    case 'config_CEM_RS'
          
        if (showGraph == 1)
            drawFluenceMaps(X_flu , Y_flu, Fluence , E_flu , 1,'Fluence before spike')
            fprintf('Fluence before spike \n')
        end
        
        %---------------------------------------
        %Compute scattering in CEM
        %---------------------------------------
        if (showGraph == 1)
            fprintf('Computing lateral scattering in the spikes \n')
        end

        RSThick = sum(Plan.Beams(b).RSinfo.RSslabThickness); % Total thickness (mm) of the range shifter
        RSposition = Plan.Beams(b).RSinfo.IsocenterToRangeShifterDistance + RSThick; %Z IEC gantry position of the upstream side of the range shifter (NB IsocenterToRangeShifterDistance is the downstream side of the RS)
        distance = CEF_Zg - max(CEF_thickness) - RSposition;  %Distance (mm) between CEM exit and entrance of the range filter
                                                          %(CEF_Zg - max(CEF_thickness) ) is the position of the downstream side of the CEM on the Z IEC  gantry axis
                                                          %max(CEF_thickness) is the height of the highest peak
                                                          %The scattering in RS and water tank will be added later
        if (distance < 0)
            CEF_Zg
            maxCEFThickness = max(CEF_thickness)
            RSposition
            fprintf('CEF thickness : %f mm \n', max(CEF_thickness))
            fprintf('Missing space : %f mm\n', distance)
            error('Collision between CEM and range shifter')
        end


        if ~isempty(maskLayer)
            %Check that the thickness are ordered in increasing size
            if((prod(diff(CEF_thickness))>0)==0)
              error('thickness must be sorted in ascending order')
            end
        end


        diverg = mean([Plan.Beams(b).sigmaAtCEF.BDL.Divergence1x , Plan.Beams(b).sigmaAtCEF.BDL.Divergence1y]); %Average the beam divergence in X and Y
        scattering = Plan.Spike.Scattering;


        switch scattering
            case 'SAM'
                fprintf('Using semi-analytical model, based on MC, for CEM spike scattering... \n');
                sigmas.sigmaCEF = round(getSpotSigmaSAM(Emax , CEF_thickness ) ./ pxlSize); %Get the spot sigma (in pixels) of the semi analytical model
            case 'Moliere'
                if isempty(sigmas.sigmaCEF)
                    fprintf('Using Moliere angle for CEM spike scattering... \n');
                    [~, ~ , thtM] = getSpotSigmaMoliere(Emax , CEF_thickness , distance, Plan.Spike.MaterialID); %Get the beam divergence from Moliere scattering
                    %Add the variance of each process to obtain the sigma of the Gaussian kernel for the convolution
                    % Moliere scattering + intrinsic beam divergence at the base of the CEF + add broadening due to transmission through |distance|
                    sigmas.sigmaCEF = sqrt( (thtM .* CEF_thickness).^2 + (thtM .* distance ).^2 + (diverg .* CEF_thickness).^2 + (diverg .* distance ).^2); %sigma in mm
                    sigmas.sigmaCEF = round(sigmas.sigmaCEF ./ pxlSize); %sigma in pixels
                end
            case 'User'
                if isempty(sigmas.sigmaCEF)
                    fprintf('Using User-defined angle for CEM spike scattering... \n');
                    thtU = getUserScatteringAngle(CEF_thickness, Plan.Spike.MaterialID);
                    %Add the variance of each process to obtain the sigma of the Gaussian kernel for the convolution
                    % User scattering + intrinsic beam divergence at the base of the CEF + add broadening due to transmission through |distance|
                    sigmas.sigmaCEF = sqrt( (thtU .* CEF_thickness).^2 + (thtU .* distance ).^2 + (diverg .* CEF_thickness).^2 + (diverg .* distance ).^2); %sigma in mm
                    sigmas.sigmaCEF = round(sigmas.sigmaCEF ./ pxlSize); %sigma in pixels
                end
            otherwise
                if isempty(sigmas.sigmaCEF)
                    fprintf('Using Moliere angle for CEM spike scattering... \n');
                    [~, ~ , thtM] = getSpotSigmaMoliere(Emax , CEF_thickness , distance, Plan.Spike.MaterialID); %Get the beam divergence from Moliere scattering
                    %Add the variance of each process to obtain the sigma of the Gaussian kernel for the convolution
                    % Moliere scattering + intrinsic beam divergence at the base of the CEF + add broadening due to transmission through |distance|
                    sigmas.sigmaCEF = sqrt( (thtM .* CEF_thickness).^2 + (thtM .* distance ).^2 + (diverg .* CEF_thickness).^2 + (diverg .* distance ).^2); %sigma in mm
                    sigmas.sigmaCEF = round(sigmas.sigmaCEF ./ pxlSize); %sigma in pixels
                end
        end

        [Fluence , sigmas.sigmaCEF]  = convScatter(Fluence, sigmas.sigmaCEF , maskLayer);


        if (showGraph == 1)
            drawFluenceMaps(X_flu , Y_flu, Fluence , E_flu , 2 , 'Fluence at entrance of RS')
            fprintf('Fluence at entrance of RS \n')
        end

        %---------------------------------------
        %Compute scattering in range shifter and water tank
        % The range shifter is made of water
        %---------------------------------------
        WaterThickeness = Plan.Beams(b).Iso2Skin + Plan.Beams(b).RSinfo.RangeShifterWET;  %WET of the range shifter + patient
        [Fluence , sigmas.water]= fluenceAfterWater(Fluence , Emax , WaterThickeness , pxlSize , sigmas.water);

        if (showGraph == 1)
            drawFluenceMaps(X_flu , Y_flu, Fluence , E_flu , 3 , 'Fluence in water tank')
            fprintf('Fluence in water tank \n')
        end

        %-------------------------------------------------
        % Prediction of energy spread = range straggling
        %-------------------------------------------------
        if (showGraph == 1)
            fprintf('Computing range straggling \n')
        end

        Fluence = convEspread(E_flu , Fluence);

        if (showGraph == 1)
            drawFluenceMaps(X_flu , Y_flu, Fluence , E_flu , 4 ,'Fluence after range straggling')
        end

        %factors for lense effect
        if (~isempty(SAD))
            if (showGraph == 1)
                fprintf('Computing geometrical projection \n')
            end
            %Take intop account the beam divergence
            Xmeas = X_flu .* (SAD(1) - Meas_Zg) ./ (SAD(1) - CEF_Zg);
            Ymeas = Y_flu .* (SAD(2) - Meas_Zg) ./ (SAD(2) - CEF_Zg);
        else
            %Ignore the beam divergence
            Xmeas = X_flu;
            Ymeas = Y_flu;
        end
    case 'config_RS_CEM'
              
        if (showGraph == 1)
            drawFluenceMaps(X_flu , Y_flu, Fluence , E_flu , 1, 'Fluence before RS')
            fprintf('Fluence before RS \n')
        end
        
        %---------------------------------------
        %Compute scattering in range shifter and water tank
        % The range shifter is made of water
        %---------------------------------------
        WaterThickeness = Plan.Beams(b).Iso2Skin + Plan.Beams(b).RSinfo.RangeShifterWET;  %WET of the range shifter + patient
        [Fluence , sigmas.water]= fluenceAfterWater(Fluence , Emax , WaterThickeness , pxlSize , sigmas.water);

        if (showGraph == 1)
            drawFluenceMaps(X_flu , Y_flu, Fluence , E_flu , 3 , 'Fluence in water tank')
            fprintf('Fluence in water tank \n')
        end

        %-------------------------------------------------
        % Prediction of energy spread = range straggling
        %-------------------------------------------------
        if (showGraph == 1)
            fprintf('Computing range straggling \n')
        end

        Fluence = convEspread(E_flu , Fluence);

        if (showGraph == 1)
            drawFluenceMaps(X_flu , Y_flu, Fluence , E_flu , 4 ,'Fluence after range straggling')
        end

        %factors for lense effect
        if (~isempty(SAD))
            if (showGraph == 1)
                fprintf('Computing geometrical projection \n')
            end
            %Take intop account the beam divergence
            Xmeas = X_flu .* (SAD(1) - Meas_Zg) ./ (SAD(1) - CEF_Zg);
            Ymeas = Y_flu .* (SAD(2) - Meas_Zg) ./ (SAD(2) - CEF_Zg);
        else
            %Ignore the beam divergence
            Xmeas = X_flu;
            Ymeas = Y_flu;
        end
        
        
        %---------------------------------------
        %Compute scattering in CEM
        %---------------------------------------
        if (showGraph == 1)
            fprintf('Computing lateral scattering in the spikes \n')
        end

        RSposition = Plan.Beams(b).RSinfo.IsocenterToRangeShifterDistance; %Z IEC gantry position of the upstream side of the range shifter (NB IsocenterToRangeShifterDistance is the downstream side of the RS)
        CEMPosition = CEF_Zg + max(CEF_thickness);
        distance = RSposition - CEMPosition;  %Distance (mm) between RS exit and entrance of the CEM
                                                          %(CEF_Zg + max(CEF_thickness) ) is the position of the downstream side of the CEM on the Z IEC  gantry axis
                                                          %max(CEF_thickness) is the height of the highest peak
                                                          %The scattering in RS and water tank was added before

        if (distance < 0)
            CEF_Zg
            maxCEFThickness = max(CEF_thickness)
            RSposition
            CEMPosition
            fprintf('CEF thickness : %f mm \n', max(CEF_thickness))
            fprintf('Missing space : %f mm\n', distance)
            error('Collision between RS and CEM')
        end


        if ~isempty(maskLayer)
            %Check that the thickness are ordered in increasing size
            if((prod(diff(CEF_thickness))>0)==0)
              error('thickness must be sorted in ascending order')
            end
        end


        diverg = mean([Plan.Beams(b).sigmaAtCEF.BDL.Divergence1x , Plan.Beams(b).sigmaAtCEF.BDL.Divergence1y]); %Average the beam divergence in X and Y
        scattering = Plan.Spike.Scattering;


        switch scattering
            case 'SAM'
                fprintf('Using semi-analytical model, based on MC, for CEM spike scattering... \n');
                sigmas.sigmaCEF = round(getSpotSigmaSAM(Emax , CEF_thickness ) ./ pxlSize); %Get the spot sigma (in pixels) of the semi analytical model
            case 'Moliere'
                if isempty(sigmas.sigmaCEF)
                    fprintf('Using Moliere angle for CEM spike scattering... \n');
                    [~, ~ , thtM] = getSpotSigmaMoliere(Emax , CEF_thickness , distance, Plan.Spike.MaterialID); %Get the beam divergence from Moliere scattering
                    %Add the variance of each process to obtain the sigma of the Gaussian kernel for the convolution
                    % Moliere scattering + intrinsic beam divergence at the base of the CEF + add broadening due to transmission through |distance|
                    sigmas.sigmaCEF = sqrt( (thtM .* CEF_thickness).^2 + (thtM .* distance ).^2 + (diverg .* CEF_thickness).^2 + (diverg .* distance ).^2); %sigma in mm
                    sigmas.sigmaCEF = round(sigmas.sigmaCEF ./ pxlSize); %sigma in pixels
                end
            case 'User'
                if isempty(sigmas.sigmaCEF)
                    fprintf('Using User-defined angle for CEM spike scattering... \n');
                    thtU = getUserScatteringAngle(CEF_thickness, Plan.Spike.MaterialID);
                    %Add the variance of each process to obtain the sigma of the Gaussian kernel for the convolution
                    % User scattering + intrinsic beam divergence at the base of the CEF + add broadening due to transmission through |distance|
                    sigmas.sigmaCEF = sqrt( (thtU .* CEF_thickness).^2 + (thtU .* distance ).^2 + (diverg .* CEF_thickness).^2 + (diverg .* distance ).^2); %sigma in mm
                    sigmas.sigmaCEF = round(sigmas.sigmaCEF ./ pxlSize); %sigma in pixels
                end
            otherwise
                if isempty(sigmas.sigmaCEF)
                    fprintf('Using Moliere angle for CEM spike scattering... \n');
                    [~, ~ , thtM] = getSpotSigmaMoliere(Emax , CEF_thickness , distance, Plan.Spike.MaterialID); %Get the beam divergence from Moliere scattering
                    %Add the variance of each process to obtain the sigma of the Gaussian kernel for the convolution
                    % Moliere scattering + intrinsic beam divergence at the base of the CEF + add broadening due to transmission through |distance|
                    sigmas.sigmaCEF = sqrt( (thtM .* CEF_thickness).^2 + (thtM .* distance ).^2 + (diverg .* CEF_thickness).^2 + (diverg .* distance ).^2); %sigma in mm
                    sigmas.sigmaCEF = round(sigmas.sigmaCEF ./ pxlSize); %sigma in pixels
                end
        end
        
        [Fluence , sigmas.sigmaCEF]  = convScatter(Fluence, sigmas.sigmaCEF , maskLayer);

        if (showGraph == 1)
            drawFluenceMaps(X_flu , Y_flu, Fluence , E_flu , 2 , 'Fluence at entrance of CEM')
            fprintf('Fluence at entrance of CEM \n')
        end
        
end %switch

end %function


%================================================
% Account for the energy spread due to proton range straggling
% using the semi anlytical model
%
% INPUT
% |E_flu| -_SCALAR VECTOR_- |E_flu(k)| Energy (MeV) of the k-th cylinder / layer
% |Fluence(x,y,E)| -_SCALAR MATRIX_-  |Fluence(x,y,E)| Incident proton fluence at posiiton (x,y) for proton of energy E
%
% OUTPUT
% |Fluence(x,y,E)| -_SCALAR MATRIX_-  Fluence with energy spreading |Fluence(x,y,E)| Incident proton fluence at posiiton (x,y) for proton of energy E
%================================================
function FluenceOUT = convEspread(E_flu , Fluence )

    dE = diff(E_flu); %Work out wether the  E_flu vector is sorted in increasing or decreasing energy

    FluenceOUT = zeros(size(Fluence));
    lexan = materialDescription('lexan');
    R = energy2range(E_flu, lexan.alpha,lexan.p); %Range (cm)

    [~ , pE ] = getSAmodelParam(); %Get the polynomial coefficient of the semi nalytical model
    sigma = polyval(pE,R);

    N_thickn = numel(E_flu); %Number of energy filters
    FluenceOUT = zeros(size(Fluence)); %Fluence map after convolution for range straggling

    %Convolve each energy layer to its 2 neighbouring layers
    for k = 1:N_thickn
        G2 = exp( -(R-R(k)).^2 ./ (2*sigma(k).^2) ); %pseudo-Gaussian kernel for range straggling
        G2 = G2 ./ sum(G2,'all'); %Normalise convolution kernel. It is not a real Gaussian (sigma change at each layer) so we need to numerically integrate the function

        %Convolution between enrgy layers. Only use the 1st neighbourgh with lower energy
        %The range straggling is of the order of 2mm. The range distance between layer is ~5mm => no need to go beyond 1 layer down
        if dE(1) > 0
          %Increasing energy. Leak into k-1
          if ((k-1) > 0)
            FluenceOUT(:,:,k) = Fluence(:,:,k) .* G2(k);
            FluenceOUT(:,:,k-1) = Fluence(:,:,k-1) + Fluence(:,:,k) .* G2(k-1);
          end
        else
          %Decreasing energy. Leak into k+1
          if ((k+1) <= N_thickn)
            FluenceOUT(:,:,k) = Fluence(:,:,k) .* G2(k);
            FluenceOUT(:,:,k+1) = Fluence(:,:,k+1) + Fluence(:,:,k) .* G2(k+1);
          end
        end
    end
end



%======================================================
% Get the fluence contribution from neighbouring spots on the current spot
%
% INPUT
% |Xn| -_SCALAR VECTOR_-X coordinate (mm) of the i-th spot with respect to the current spike
% |Yn| -_SCALAR VECTOR_-Y coordinate (mm) of the i-th spot with respect to the current spike
% |w| -_SCALAR VECTOR_- Weight of the i-th spot
% |sigmaAtCEF| -_STRUCTURE_- Information on spot shape
%    * |sigmaAtCEF.Sx| -_SCALAR_- Sigma (mm) of the lateral Gaussian dose distribution in a PBS spot at the base of the CEF along the X axis
%    * |sigmaAtCEF.Sy| -_SCALAR_- Sigma (mm) of the lateral Gaussian dose distribution in a PBS spot at the base of the CEF along the Y axis
% |pts_Spike| -_SCALAR VECTOR_- |pts_Spike(i,:)=[x,y]| Coordinates (as defined by |x_axis| and |y_axis|) of the i-th point in the grid
% |CEF_Zg| -_SCALAR_- Distance (mm) between isocentre and the base of the CEM
%
% OUTPUT
% |C| -_SCALAR VECTOR_- |C(x,y,E)| Fluence at position |Xn(x)|,|Yn(y)| of the proton going throught the E-th step

%======================================================
function C = getWofNeighhbourgh(Xn , Yn , w , sigmaAtCEF , pts_Spike , CEF_Zg)

  C = [];
  for spts = 1:numel(Xn)
    if isempty(C)
      C = w(spts) .* getSpotFluence(pts_Spike , [Xn(spts),Yn(spts)] , sigmaAtCEF.BDL, CEF_Zg);
    else
      C = C + w(spts) .* getSpotFluence(pts_Spike , [Xn(spts),Yn(spts)] , sigmaAtCEF.BDL, CEF_Zg);
    end
  end

end

%=============================================
% Draw the mask for the polygonal rings
%
% |apothem_max| -_SCALAR_- Distance (mm) between centre of polygon and edge of the external polygon of the ring
% |apothem_min| -_SCALAR_- Distance (mm) between centre of polygon and edge of the internal polygon of the ring
% |nrSides| -_INTEGER_- Number of sides to the polygon. Can be 4 or 6
% |pts_Spike| -_SCALAR MATRIX_- |pts_Spike(:,i)= [x,y]| Coordinate (mm) of the i-th pixel in |Zgrid(sub2ind(SpikeSize,i))|
% |SpikeSize| -_SCALAR MATRIX_- Number of pixels [Nx,Ny] in |Zgrid|
% |r| -_SCALAR_- Ratio aX / aY of the ellipse moin axes (used for ellipses)
%
% OUTPUT
% |mask| -_SCALAR MATRIX_- |mask(x,y)| =1 if the pixel (x,y) is inside the ring. O=0 otherwise. The matrix has dimension [Nx,Ny]
% |maskIntern| -_SCALAR MATRIX_- |maskIntern(x,y)| =1 if the pixel (x,y) is in the central holde of the polygon
% |maskExtern| -_SCALAR MATRIX_- |maskExtern(x,y)| =1 if the pixel (x,y) is inside within the controur of the polygon (including the central hole)
%=============================================
function [mask , maskIntern , maskExtern]  = drawPolygonMask(apothem_max,apothem_min, nrSides, pts_Spike , SpikeSize , r)
  maskExtern = drawPolygon(nrSides, pts_Spike' , SpikeSize , apothem_max , 1 , r); %h_step = 1 to define a mask
  maskIntern = drawPolygon(nrSides, pts_Spike' , SpikeSize , apothem_min , 1 , r); %h_step = 1 to define a mask
  mask = ~~(maskExtern - maskIntern); %Remove the central polygon. Keep only the ring.
end

%=============================================
% Simulate printing error
% The width of each spike is dilated by 1 pixel in all direction
% The highest tower override smaller tower
%=============================================
function MaskOut = dilateSpike(MaskSpk)

  NbPxl = 1; %Number of pixel of the dilation. This gives the magnitude of the printer error
  se = strel('disk',1);

  MaskOut = zeros(size(MaskSpk));

  Z = unique(MaskSpk(:));
  Z = sort(Z); %sort in ascending order


  for steps = 1:numel(Z)
        mask = (MaskSpk == Z(steps)); %select the mask of the level Z
        mask = imdilate(mask,se);
        MaskOut(find(mask))=Z(steps);
  end

end
