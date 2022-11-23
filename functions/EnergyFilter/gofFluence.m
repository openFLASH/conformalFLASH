%% gofFluence
% Compute the goodness of fit function between a target fluence map |FluenceRef|
% and the fluence map generated by a spike of a CEF:
%  gof = A .* Fluence([r1,r2,...,rn]) - FluenceRef
%
%% Syntax
% |[gof , FluenceOpt] = gofFluence(x , Plan , b , SpikeIdx , sigma0 , Meas_Zg , sigmas_far , FluenceRef, apertureMask , Iso2RS , E_ref)|
%
%
%% Description
% |[gof , FluenceOpt] = gofFluence(x , Plan , b , SpikeIdx , sigma0 , Meas_Zg , sigmas_far , FluenceRef, apertureMask , Iso2RS , E_ref)| Description
%
%
%% Input arguments
% |x| -_SCALAR VECTOR_- The optimised parameters |x = [A , dr1 , dr2, ... , drn]| where the relative width of the lowers step is dr1 and the relative width of the highest step is drn
%
% |Plan| - _struct_ - MIROpt structure with updated information. See fluenceWithCEF.m for more information
%     * |Plan.Spike.intrpCTpxlSize| -_SCALAR_- Lateral resolution (mm) of the 3D printer. The steps are multiple of this value.
%
% |b| -_INTEGER_- Index of the beam in |Plan.Beams(b)| for which the fluence map is to be computed
%
% |SpikeIdx| -_INTEGER_- Index of the spike in |Plan.Beams(b).RidgeFilter(SpikeIdx)| that should be optimizes
%
% |Meas_Zg| -_SCALAR_- Z coordinate in IEC-Gantry (mm) of the plane in which the fluence is measured
%
% |showGraph| -_LOGICAL_- [Optional: defauly = false] If true, display the graphs ofthe fluence at different steps i nthe computations
%
% |sigmas_far| -_SCALAR VECTOR_- [OTPIONAL: if [], the sigma are computed by the function]. sigma (mm) of the lateral spread of the beam at |distance| of the scatterer
%
% |FluenceRef(x,y,E)| -_SCALAR MATRIX_-  Reference fluence at |distance| downstream from the CEF that the gof tries to reproduce. |Fluence(x,y,E)|  fluence at position (x,y) for proton of energy |Plan.Beams(b).Layers(E).Energy|
%
% |apertureMask| -_SCALAR MATRIX_- |apertureMask(x,y)=1| if the pixel is inside the aperture opening. The GOF is based only on the fluence inside the aperture opening. Empty if no aperture is used
%
% |Iso2RS| -_SCALAR_- Distance in WET (mm) from base of range shifter to isocentre. Includes range shiter and water tank
%
% |conv_coor| -_SCALAR MATRIX_- [OPTIONAL: used to make computation faster. If absent, the parameter is computed] |conv_coor(m,n)=i| the point |[x_grid(m),y_grid(n)]| of the cartesian coordinate system has polar coordinates |R(i)| and |Phi(i)| in polar coordinates
%
%% Output arguments
%
% |gof| - _SCALAR_ -  Goodness of fit. Sum of the square differences
%
% |FluenceOpt|  -_SCALAR MATRIX_-  Fluence of the optimised spike at |distance| downstream from the CEF that the gof tries to reproduce. |Fluence(x,y,E)|  fluence at position (x,y) for proton of energy |Plan.Beams(b).Layers(E).Energy|
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [gof , FluenceOpt] = gofFluence(x , Plan , b , SpikeIdx , sigma0 , Meas_Zg , sigmas_far , FluenceRef, apertureMask , Iso2RS , E_ref)

    sFluenceRef = size(FluenceRef);

    if isempty(apertureMask)
      %If the aperture mask is empty, then use all fluence pixels in the computation of GOF
      apertureMask = ones(sFluenceRef);
    else
      %An aperture mask is provided. Repeat for eack energy layer
      apertureMask = repmat(apertureMask,1,sFluenceRef(3)); %repeat for every energy layer
      apertureMask = reshape(apertureMask,sFluenceRef);
    end

    pencil_x = [Plan.Beams(b).RidgeFilter(:).x_centre]; %Position of the centre of the spikes is also the centre of the PBS spots
    pencil_y = [Plan.Beams(b).RidgeFilter(:).y_centre];
    weight_table = [Plan.Beams(b).RidgeFilter(:).w]; % w(spt) is the weight of the spt-th PBS spot. The value is proportional to the number of proton at **maximum** energy

    W0 = Plan.Beams(b).RidgeFilter(SpikeIdx).w_step; % w_step(l) weight of the l-th energy layer to the spot |SpikeIdx|. Proportional to number of protons of **maximum** energy
    W0 = W0 ./ sum(W0); %Normalise the reference weights of the IMPT plan

    %Plan.Beams(b).RidgeFilter  = packRgofFluence(x,Plan.Beams(b).RidgeFilter(SpikeIdx), Plan.Spike.LateralStep); %Update the width of the steps using the X variable
    Plan.Beams(b).RidgeFilter = packRgofFluence(x,Plan.Beams(b).RidgeFilter(SpikeIdx), 0); %Update the width of the steps using the X variable

    %Plan contains the definition of one single spike after the call to packRgofFluence
    if Plan.showGraph
      showGraph = 2;
    else
      showGraph = 0;
    end
    [FluenceSpk , X_far , Y_far , E_filter , ~ , ~ , maskLayer] = fluenceWithCEF(Plan , b , sigma0 , pencil_x , pencil_y  , weight_table , [] ,  Meas_Zg , 2, sigmas_far, 0, 'config_RS_CEM');

    mask = ~~sum(maskLayer,3); %The gof is computed only in the lattice cell fo the spike

    %If there are missing energy layers in this spike, resort energy layers to match reference
    Fluence = zeros(numel(X_far),numel(Y_far), numel(E_ref));
    W0ext = zeros(numel(E_ref),1); %There are missing energy layers in the W0 reference

    for ind_E = 1:numel(E_filter)
        idx = find(E_ref == E_filter(ind_E)); %Find the index of the energy layer corrsponding to cylinder j
        Fluence    (: , : , idx) = FluenceSpk(:,:,ind_E);
        W0ext(idx) = W0(ind_E); %Create a reference weight vector with all energy layers
    end

    W = squeeze(sum(Fluence .* mask,[1,2])); %The weight of each energy layer is the integral of the fluence of the energy layer
    %W = squeeze(sum(Fluence .* mask .* apertureMask,[1,2])); %The weight of each energy layer is the integral of the fluence of the energy layer
    W = W ./ sum(W); %Normalise the weight
    W0 = W0ext; %Fill the missing energy layers in the reference with zeros

    figure(151)
    hold off
    plot(E_ref , W0 , '-ob')
    hold on
    plot(E_ref , W , '-or')
    xlabel('Energy (MeV)')
    ylabel('Weight (relative)')
    title('GOF : Bragg peaks weights')
    grid on

    gof = sum((W - W0).^2);


end
