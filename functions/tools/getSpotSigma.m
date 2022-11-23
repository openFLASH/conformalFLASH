%% getSpotSigma
% Determine the sigma (lateral spread) of the deepest BP which is part of a SOBP beamlet.
% For the specified beam, find the SOBP beamlet that is closest to a specified coordinate
% Then estimate the spot sigma at the specified depth:
%  * either at the depth of maximum dose of IDD for the BP with highest energy in the SOBP
%  * or at the level of the external contour
%
% Alternatively, the spot sigma can be estiamted form the MCsquare beam data library at ta specified distancefrom isocentre
%
% The sigma is estimated by fitting a bi-normal Gaussian
%
%% Syntax
% |[spotSigma , Diso2Meas ]  = getSpotSigma(Plan , b , RefPoint)|
%
% |[spotSigma , Diso2Meas ]  = getSpotSigma(Plan , b , RefPoint , Location , miscData)|
%
%
%% Description
% |[spotSigma , Diso2Meas ]  = getSpotSigma(Plan , b , RefPoint)| Compute spot sigma at the maximum of the deepest Bragg peak
%
% |[spotSigma , Diso2Meas ]  = getSpotSigma(Plan , b , RefPoint , 'skin' , ROI)| Compute the spot sigma at the surface of the ROI mask
%
% |[spotSigma , Diso2Meas ]  = getSpotSigma(Plan , b , RefPoint , 'PLANE' , Dist2Iso)| Compute the spot sigma in a plane at a specified distance from isocentre using the data from the MCsquare Beam data library
%
% |[spotSigma , Diso2Meas ]  = getSpotSigma(Plan , b , RefPoint , 'NOZZLE')| Compute the spot sigma at the nozzle eexit using data from the MCsquare BDL
%
%
%% Input arguments
% |Plan| - _STRUCT_ - MIROpt structure where all the plan parameters are stored.
%
% |b| -_INTEGER_- Beam number
%
% |RefPoint| -_SCALAR VECTOR_- |RefPoint=[x,y]| Coordinate (mm) in the IEC gantry CS of the point near which the SOBP beamlet must be located
%
% |Location| -_STRING_- [OTPIONAL: defaul 'BP'] Position of the plane in which the spot sigma is computed. 'BP' , 'skin', 'PLANE'
%
% |ROI| -_SCALAR MATRIX_- [OTPIONAL. Only necessary for |Location = 'skin'|]. |ROI(x,y,z)=1| if the point belongs to the body. MAsk defining of the body ocntour
%
% |Dist2Iso| -_SCALAR_- [OTPIONAL. Only necessary for |Location = 'PLANE'|]. Distance from isocentre to measurement plane = Zgantry coordinate of the plane
%
%% Output arguments
%
% |spotSigma| - _SCALAR_ - [For |Location = 'skin'| and |Location = 'BP'|] Sigma (mm) of the lateral Gaussian profile
% |spotSigma| - _STRUCTURE_ - [For |Location = 'PLANE' and 'NOZZLE'|] Sigma of the lateral Gaussian profile
%       * |spotSigma.Sx| -_SCALAR_- Sigma (mm) of the lateral Gaussian dose distribution in a PBS spot at the base of the CEF along the X axis
%       * |spotSigma.Sy| -_SCALAR_- Sigma (mm) of the lateral Gaussian dose distribution in a PBS spot at the base of the CEF along the Y axis
%       * |spotSigma.r| -_SCALAR_- Correlation between X and Y
%       * |spotSigma.sigma| -_SCALAR_-  lateral sigma (mm) of the PBS spot along the main axis of the ellipse
%       * |spotSigma.BDL| -_STRUCTURE_- [only if |Location == 'NOZZLE'| ] Spot information from the Beam Data Library. See |getSpotFromBDL| for more information
%
% |Diso2Meas| -_INTEGER_- Distance (mm) between the isocentre and the centre of the plane in which the spot sigma is reported
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [spotSigma , Diso2Meas ] = getSpotSigma(Plan , b , RefPoint , Location , misData)

  if nargin < 4
    Location = 'BP';
  end

  %Create a coordinate grind on which the Gaussian profile will be estimated
  pxlSize = min(Plan.DoseGrid.resolution); %smallest pixel dimension
  Xaxis = -40:pxlSize:40; %mm
  Yaxis = -40:pxlSize:40;
  [PtsG , Xg , Yg] = makeListOfCoordinates(Xaxis , Yaxis); %Create a meshgrid in the IEC gantry CS. This will cover most of the Gaussian shape
  NbPts = numel(Xg);


  switch Location
    case 'NOZZLE'
        E0 = max([Plan.Beams(b).Layers(:).Energy],[],'all'); %The maximum energy of incoming protons
        Diso2Meas = getNozzle2IsoDistanceFromBDL(Plan.BDL);
        SpotData = getSpotFromBDL(Plan.BDL , E0 );
        spotSigma.BDL = SpotData;
        [DoseMap , Dmax ] = getSpotFluence(PtsG , [0,0] , SpotData , Diso2Meas);

    case 'PLANE'

        E0 = max([Plan.Beams(b).Layers(:).Energy],[],'all'); %The maximum energy of incoming protons
        Diso2Meas = misData;
        SpotData = getSpotFromBDL(Plan.BDL , E0 );
        spotSigma.BDL = SpotData;
        [DoseMap , Dmax ] = getSpotFluence(PtsG , [0,0] , SpotData, Diso2Meas);

        % %Compute spot profile from the sum of the 2 Gaussian defined in BDL
        % DoseMap = biNorm(PtsG , SpotData.Weight1 .* 2 .* pi .* SpotData.SpotSize1x .* SpotData.SpotSize1y,  [0,0] , SpotData.SpotSize1x , SpotData.SpotSize1y , 0) ...
        %         + biNorm(PtsG , SpotData.Weight2 .* 2 .* pi .* SpotData.SpotSize2x .* SpotData.SpotSize2y,  [0,0] , SpotData.SpotSize2x , SpotData.SpotSize2y , 0);
        % Dmax  = max(DoseMap,[],'all'); %Find the highest dose
    otherwise
      %case 'skin' and 'BP'

      %Identify the spots composing the SOBP closest to the specified point
      [spotPosition , weight2spot ] =  collectSpotsinBeamlets(Plan);
      [Dist2 , idxCentre] = min(getDistances2(spotPosition{b} , RefPoint)); %Find the spot closest to specified point for b-th beam
      Pij = Plan.Scenario4D(1).RandomScenario(Plan.rr_nominal).RangeScenario(Plan.rs_nominal).P ; %The dose influence matrix for the nominal case of the first breathin phase

      if (size(Pij,2) ~= size(weight2spot,1))
        %Sanity check: make sure that the Pij matrix was computed with this plan
        %If the test failed, this indicate that the Pij matrix was not generated with the Plan provided
        size(Pij,2)
        size(weight2spot,1)
        error('Matrix Pij does not have same number of spots than defined in the plan')
      end

      %Compute the 3D dose map of the selected BP
      w = zeros(1,size(Pij,2)); %Get the weight for the spot on the optical axis
      idxBPonTarget = find((weight2spot(:,2)==idxCentre).*(weight2spot(:,1)==b)); %Index of the spots forming the SOBP and located near the RefPoint
      [~, idxE] = max([Plan.Beams(b).Layers(weight2spot(idxBPonTarget,3)).Energy]);
      w(idxBPonTarget(idxE)) = 1; %Select the BP with the deepest range in the SOBP
      DoseSpot = Pij2Dose(Pij , w , Plan.DoseGrid.size); %3D dose for the selected BP. The indices of the matrix are running along the DICOM axes

      %Take a 2D slice, orthogonal to beam axis, at the required depth
      [Dmax , wPeak ] = max(DoseSpot,[],'all','linear'); %Find the highest dose
    end


  %Compute the DoseMap for 'BP' and 'skin'. This is already done for 'PLANE'
  switch Location

    case 'BP'
      %Compute spot sigma in the plane of the maximum of the Bragg peak
      [IDD , Zg] = getIDD(Xg , Yg, Plan.Beams(b).GantryAngle , Plan.Beams(b).PatientSupportAngle , Plan.CTinfo.ImagePositionPatient ,  Plan.Beams(b).isocenter , Plan.DoseGrid.resolution , DoseSpot , Plan.showGraph); %Get integrated depth dose profile around isocentre
      [~, Zidx] = max(IDD);
      Zidx = Zidx(1); %In case there are several maxima
      PtsG = [Xg(:) , Yg(:) , ones(NbPts,1) .* Zg(Zidx)]; %Physical coordinate (gantry) of the point in the slice perpendicular to proton axis.
                                            %Select point at max of IDD
      PmaxG = [0 , 0 , Zg(Zidx) ./ Plan.CTinfo.Spacing(3) , 1]'; %Coordinate of the maximum of the PBS spot (in pixel coordiante in IEC gantry)
      M = matDICOM2IECgantry(Plan.Beams(b).GantryAngle , Plan.Beams(b).PatientSupportAngle);  %4x4 rotation matrix to go from IEC gantry to DICOM : dicom = R * IEC_G
      R = inv(M);
      PmaxDCM = R * PmaxG; %in DICOM CS, in pixel unit
      Diso2Meas = abs(Zg(Zidx)); %mm Distance along IEC gantry axis from isocentre top position of maximum of IDD
      DoseMap = getSlice(PtsG , Plan.Beams(b).GantryAngle , Plan.Beams(b).PatientSupportAngle , Plan.CTinfo.ImagePositionPatient , Plan.Beams(b).isocenter , Plan.DoseGrid.resolution , DoseSpot);

    case 'skin'
      %Compute spot sigma at the skin level
      %Make a mask defining the skin
      ExternalROI_ID = getROIByName(misData, Plan.ExternalROI);
      body = misData(ExternalROI_ID).mask3D.value;
      SE = strel('sphere',5);
      skin = (~imerode(body,SE)) .* body; %Mask defining the position of the skin
      [Xsk,Ysk,Zsk] =  ind2sub(Plan.DoseGrid.size , find(skin)); %index of the voxels of the skin
      PtSkDC =  PXLindex2DICOM([Xsk,Ysk,Zsk] , Plan.CTinfo.Spacing , Plan.CTinfo.ImagePositionPatient );
      PtSkDC = PtSkDC'; %DICOM coordinate (mm) of the voxels
      M = matDICOM2IECgantry(Plan.Beams(b).GantryAngle , Plan.Beams(b).PatientSupportAngle , Plan.Beams(b).isocenter); %Rotate around isocentre
      PtskG = M * PtSkDC; %Coordinate (mm) of the skin voxels in IEC gantry
      [~, nearAxis] = find(sum(PtskG(1:2,:).^2 , 1) <= Plan.Beams(b).SpotSpacing.^2); %Find the voxels closest to optical axis of beamlet. Their distance to beam aixs is less than SpotSpacing
      [Zgmax , idxSkin] = max(PtskG(3,nearAxis)); %The point at the skin surface has the highest Zg (it is the further away from isocentre and towards the source) and is near the optical axis
      PtsG = [Xg(:) , Yg(:) , ones(NbPts,1) .* (Zgmax(1) - 10)]; %Shift 10mm back to be sure to be inside some tissue where dose is deposited in order to compute the sigma
      Diso2Meas = Zgmax; %Record the (un-shifted) distance to skin so that the aperture and CEF are located at the correct depth
      DoseMap = getSlice(PtsG , Plan.Beams(b).GantryAngle , Plan.Beams(b).PatientSupportAngle , Plan.CTinfo.ImagePositionPatient , Plan.Beams(b).isocenter , Plan.DoseGrid.resolution , DoseSpot);

  end


  %Lateral dose profile
  if Plan.showGraph
    figure(200)
    size(DoseMap)
    size(Xg)
    contour(Xg , Yg , reshape(DoseMap,size(Xg)),'ShowText','on')
    xlabel('X gantry (mm)')
    ylabel('Y gantry (mm)')
    title(['Lateral dose profile at ',Location,': MCsquare'])
    grid on
  end

  %Fit a bi-normal gaussian to the dose distribution
  [p,yf] = fitDoubleGauss(DoseMap , PtsG);

  if Plan.showGraph
      figure(201)
      contour(Xg , Yg , reshape(yf,size(Xg)),'ShowText','on')
      xlabel('X gantry (mm)')
      ylabel('Y gantry (mm)')
      title(['Lateral dose profile at ',Location,': Gaussian fit'])
      grid on

      figure(202)
      plot(100.*(DoseMap-yf)./Dmax,'o')
      xlabel('Point')
      ylabel('Error (%)')
      title('Residual')
  end

  %Covariance matrix
  [spotSigmaSc , sx , sy , r ] = covMatrix2sigma(p);

  switch Location
    case 'PLANE'
      spotSigma.sigma = spotSigmaSc; %| -_SCALAR_-  lateral sigma (mm) of the PBS spot along the main axis of the ellipse
      spotSigma.Sx = sx; %| -_SCALAR_- Sigma (mm) of the lateral Gaussian dose distribution in a PBS spot at the base of the CEF along the X axis
      spotSigma.Sy = sy; %| -_SCALAR_- Sigma (mm) of the lateral Gaussian dose distribution in a PBS spot at the base of the CEF along the Y axis
      spotSigma.r = r; %| -_SCALAR_- Correlation between X and Y

      case 'NOZZLE'
        spotSigma.sigma = spotSigmaSc; %| -_SCALAR_-  lateral sigma (mm) of the PBS spot along the main axis of the ellipse
        spotSigma.Sx = sx; %| -_SCALAR_- Sigma (mm) of the lateral Gaussian dose distribution in a PBS spot at the base of the CEF along the X axis
        spotSigma.Sy = sy; %| -_SCALAR_- Sigma (mm) of the lateral Gaussian dose distribution in a PBS spot at the base of the CEF along the Y axis
        spotSigma.r = r; %| -_SCALAR_- Correlation between X and Y
     otherwise
      %Just return spotSigma
      spotSigma = spotSigmaSc;
  end


end
