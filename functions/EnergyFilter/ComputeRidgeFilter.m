%% ComputeRidgeFilter
%Compute the Conformal Energy Filter (CEF) shape from the depth dose distribution along the different rays paralell to beam axis
%
%% Syntax
% |Plan = ComputeRidgeFilter(Plan, ROI, handles)|
%
%
%% Description
% |Plan = ComputeRidgeFilter(Plan, ROI, handles)| Description
%
%
%% Input arguments
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are
% stored. The following data must be present in the structure:
%   * |Plan.Beams(b).CEFbaseWET| -_SCALAR_- WET (mm) of the base of the CET acting as a range shifter
%   * |Plan.Beams(b).CEFbaseThickness| -_SCALAR_- Thickness (mm) of the base of the CEF
%   * |Plan.Beams(b).RSinfo| -_STRUCT_- Information about the range shifter
%       * |Plan.Beams(b).RSinfo.NbSlabs| -_SCALAR_- Number of slabs of range shifter
%       * |Plan.Beams(b).RSinfo.RSslabThickness|  -_SCALAR VECTOR_- |RSslabThickness(s)| Thickness (mm) of the s-th slab of the range shifter
%       * |Plan.Beams(b).RSinfo.RangeShifterMaterial| -_STRING_- Material of the range shifter
%       * |Plan.Beams(b).RSinfo.RangeShifterWET| -_SCALAR_- water equivalent thickness (mm) of the range shifter
%   * |Plan.Spike.WET| -_SCALAR_- Relative water equivalent thickness of the Conformal Energy Filter material
%   * |Plan.Spike.min_thickness| -_SCALAR_- Thickness (mm) in of the base on which the spikes are built. The base has the same |R_WET| as the spikes
%   * |Plan.Beams(b).sigmaAtCEF| -_STRUCTURE_- Sigma of the lateral Gaussian dose distribution in a PBS spot at the base of the CEF
%       * |Beams.sigmaAtCEF.Sx| -_SCALAR_- Sigma (mm) of the lateral Gaussian dose distribution in a PBS spot at the base of the CEF along the X axis
%       * |Beams.sigmaAtCEF.Sy| -_SCALAR_- Sigma (mm) of the lateral Gaussian dose distribution in a PBS spot at the base of the CEF along the Y axis
%       * |Beams.sigmaAtCEF.r| -_SCALAR_- Correlation between X and Y
%       * |Beams.sigmaAtCEF.sigma| -_SCALAR_- Sigma (mm) of the lateral Gaussian dose distribution in a PBS spot at the base of the CEF (assuming circular spot)
%   * |Plan.Beams(b).RSinfo.R_max| -_SCALAR_- Range (cm) in water to reach the distal surface of the PTV
%   * |Plan.Beams(b).BlockData|  -_SCALAR MATRIX_- |BlockData(i,:)=[x,y]|  Coordinates (mm) of the i-th point defining the contour of the aperture block projected onto the machine isocentric plane in the IEC BEAM LIMITING DEVICE coordinate system for the b-th beam.
%                                     Present only if |Plan.Beams(b).ApertureBlock==1|
%   * |Plan.Beams(b).Layers(L).Energy| -_SCALAR_- Energy (MeV) of the L-th layer
%     * |Plan.Beams(b).Layers(L).SpotPositions(k)| -_SCALAR VECTOR_- =[x,y] position (mm) of the k-th spot in the layer
%     * |Plan.Beams(b).Layers(L).SpotWeights(k)| -_SCALAR VECTOR_-  Weight (per fraction) of the k-th spot in the layer
%
% |ROI| - _struct_ - MIROpt structure containing information about all
%           volumes in the RTSTRUCT file. The following data must be present in the structure:
%     * |ROI(i).mask3D.value|- _array_ - |ROI(i).mask3D.value(x,y,z)=1| if the voxel at (x,y,z) is located inside the RT struct
%
% |handles| -_STRUCTURE_- REggui data handle. The CT scan is stored in |handles| in the image with name |Plan.CTname|.
%
%
%% Output arguments
%
% |Plan| - _struct_ - MIROpt structure with updated information:
%   * |Plan.Spike.WET| -_SCALAR_- Relative water equivalent thickness of the printed material
%   * |Plan.Spike.min_thickness| -_SCALAR_- Thickness (mm) in of the base on which the spikes are built. The base has the same |R_WET| as the spikes
%   * |Plan.Beams(b).Range| -_SCALAR_- Range (cm) in water of the spots to deliver on the CEF + range shifter
%   * |Plan.Beams(b).RangeShifterWET| -_SCALAR_- Water equivalent thickness (mm) of the range shifter to add to the CEF
%   * |Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance| -_SCALAR_- Distance (mm) from isocentre to the base of the CEF.
%   * |Plan.Beams(b).RangeModulator.CEMThicknessData| -_SCALAR MATRIX_- |CompensatorThicknessData(x,y)| Thickness (mm) of the CEF pixel at position (x;y) in the IEC beam Limiting device CS
%   * |Plan.Beams(b).RangeModulator.CEM3Dmask| -_SCALAR MATRIX_- 3D mask of the CEF. |CEM3Dmask(x,y,z)=1| if the voxel at location (x,y,z)  in the plane of the CEF for beam b belongs to the CEF.
%                                     The z index is 1 at the base of the CEM and the z index increses going up the spikes. x and y indices are aligned with the IEC gantry CS
%   * |Plan.Beams(b).RangeModulator.Modulator3DPixelSpacing| -_SCALAR VECTOR_- |CompensatorPixelSpacing = [x,y,z]| Pixel size (mm) in the plane of the CEF for the |CompensatorThicknessData| matrix in the plane of the CEM
%   * |Plan.Beams(b).RangeModulator.ModulatorMountingPosition| -_STRING_- DICOM tag Modulator Mounting Position (300D,1017). Defines the face of the CEM used as reference for |IsocenterToRangeModulatorDistance|
%                         SOURCE_SIDE is using the downstream face of the block as a reference position for expressing the isocenter to block tray distance.
%                         PATIENT_SIDE is using the upstream face
%   * |Plan.Beams(b).RangeModulator.ModulatorOrigin| -_SCALAR VECTOR_- Physical coordinate [x,y,z] the voxel |CompensatorThicknessData(1,1)| and |hedgehog3D(1,1,1)| for beam b  in the plane of the CEF.
%   * |Plan.Beams(b).RangeModulator.FrameOfReferenceUID| -_STRING_- DICOM UID of the frame of reference for the CEM
%   * |Plan.Beams(b).RangeModulator.ReferencedRangeModulator| -_STRING_- DICOM UID of the binary file describing the CEM
%   * |Plan.Beams(b).RidgeFilter| -_VECTOR of STRUCT_- Structure describing the shape of the Conformal Energy Filter filer
%      * |Plan.Beams(b).RidgeFilter(k).x_centre| -_SCALAR_- X coordinate (mm in IEC Gantry CS) of central axis of the k-th spike in the plane of the CEF
%      * |Plan.Beams(b).RidgeFilter(k).y_centre| -_SCALAR_- Y coordinate (mm in IEC Gantry CS) of central axis of the k-th spike in the plane of the CEF
%      * |Plan.Beams(b).RidgeFilter(k).apothem| -_SCALAR_- Apothem (mm) of the base of the spike
%      * |Plan.Beams(b).RidgeFilter(k).w| -_SCALAR_- Weight of the single PBS spot of energy |max(Layers(:).Energy)|  to deliver at the base of spike = weight of the SOBP
%      * |Plan.Beams(b).RidgeFilter(k).h_outside| -_SCALAR_- Height (mm) of the base of the k-th spike
%      * |Plan.Beams(b).RidgeFilter(k).a_max| -_SCALAR VECTOR_- |RF(k).a_max(L)| Max apothem (mm) of the L-th polygon (apothem of the L-th hexagon/square)
%      * |Plan.Beams(b).RidgeFilter(k).a_min| -_SCALAR VECTOR_- |RF(k).a_min(L)| Min apothem (mm) of the L-th polygon (apothem of the (L+1)-th hexagon/square)
%      * |Plan.Beams(b).RidgeFilter(k).h_step| -_SCALAR VECTOR_- |RF(k).h_step(L)| Height (mm) of the L-th polygon from the base
%      * |Plan.Beams(b).RidgeFilter(k).w_step = w_RF(kspt,k_layer)| Weight of the step k_layer at location k_RF, expressed in unit proportional to number of protons
%     * |Plan.Beams(b).RidgeFilter(k).Energy| -_SCALAR VECTOR_- |RF(k).Energy(L)| Energy (MeV) of the protons coming out of the L-th ring from the base
%
%% REFERENCES
% [1] https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html
% [2] https://en.wikipedia.org/wiki/Brass
% [3] Zhang, R., & Newhauser, W. D. (2009). Comments on “Calculation of water equivalent thickness of materials of arbitrary density, elemental composition and thickness in proton beam irradiation.��? Physics in Medicine and Biology, 54(6), 1383–1395. https://doi.org/10.1088/0031-9155/55/9/L01
%
%% Contributors
% Authors : R. Labarbe, Lucian Hotoiu (open.reggui@gmail.com)

function Plan = ComputeRidgeFilter(Plan, ROI, handles)

if nargin < 3
  CT = [];
end

CT = Get_reggui_data(handles,Plan.CTname,'images'); %Update the CT scan with the aperture block in handles

for b = 1: size(Plan.Beams,2) %Loop for each beam
    fprintf('Beam %d \n',b)

    T_max = max([Plan.Beams(b).Layers(:).Energy]);
    T_min = min([Plan.Beams(b).Layers(:).Energy]);
    rangeMin = CSDArange('Water' , T_min); %Minimum range (cm) in water
    rangeMax  = CSDArange('Water' , T_max); %Maximum range (cm) in water
    Modul = rangeMax - rangeMin; %Modulation depth (cm) in water equivalent thickness
    [~, ~, SPR] =  getMaterialPropCT(Plan.Spike.MaterialID, Plan.ScannerDirectory); %CEM relative stopping power
    CEMThickness = 10 .* Modul ./ SPR; %Physical thickness (mm) of CEM
    param = getMachineParam(Plan.BDL);
    if (CEMThickness > param.snout.CEMmaxHeight)
      fprintf('Modulation depth WET : %f mm \n', Modul .* 10)
      fprintf('CEF thickness        : %f mm \n', CEMThickness)
      fprintf('Clearance in snout   : %f mm \n', param.snout.CEMmaxHeight)
      error('Not enough space to fit the CEM')
    end

    if ~isfield(Plan,'PlanExistentFile') || (~exist(fullfile(Plan.PlanExistentFile,['Plan_ridge_tmp',num2str(b),'.mat']),'file'))

        [spotPosition , ~ , ~ , Plan] =  collectSpotsinBeamlets(Plan);

        if (isfield(Plan.Beams(b),'GridLayout'))
          GridLayout =  Plan.Beams(b).GridLayout;
        else
          GridLayout =  'HEXAGONAL';
        end

        %Compute the magnification from isocentre plane to CEF plane
        IsocenterToRangeModulatorDistance = getIsocenterToRangeModulatorDistance(Plan.Beams(b), Plan.BDL);
        fprintf('IsocenterToRangeModulatorDistance : %f mm\n', IsocenterToRangeModulatorDistance)
        mag = min(magnification(IsocenterToRangeModulatorDistance , Plan.BDL));
        fprintf('Magnification when projecting in CEF plane : %f \n', mag)
        Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance = rounding(IsocenterToRangeModulatorDistance , Plan.Spike.intrpCTpxlSize); %This is the first guess at the position of the CEF

        %Get a first estimate of spike shape based on the ratio of surface of each step
        %Ignore the lateral scattering of protons
        [Plan.Beams(b).sigmaAtCEF , Diso2Meas ] = getSpotSigma(Plan , b , [0,0] , 'PLANE' , IsocenterToRangeModulatorDistance); %REcompute the spot sigma at the base of CEF
                %We now know the distance from isocentre to base of CEF

        if (IsocenterToRangeModulatorDistance > Plan.Beams(b).sigmaAtCEF.BDL.iso2Nozzle)
          %The base of the CEM must be retracted to a distance larger than the travel distance of the snout holder
          %There is not enough room to fit the FLASH snout between nozzle and isocentre
          fprintf('IsocenterToRangeModulatorDistance : %f mm\n', IsocenterToRangeModulatorDistance)
          fprintf('Maximum travel distance of snout : %f mm\n', Plan.Beams(b).sigmaAtCEF.BDL.iso2Nozzle)
          fprintf('There is not enough space between isocentre and nozzle to fit the FLASH snout')
          fprintf('Manually define the parameter |isocenter| in the YAML to increase gap between skin surface and nozzle')
          error('There is not enough space between isocentre and nozzle to fit the FLASH snout')
        end

        fprintf('Spot spacing : %f mm \n',Plan.Beams(b).SpotSpacing);
        fprintf('Sx at base of CEF = %f mm \n',Plan.Beams(b).sigmaAtCEF.Sx);
        fprintf('Sy at base of CEF = %f mm \n',Plan.Beams(b).sigmaAtCEF.Sy);
        fprintf('r at base of CEF = %f mm \n',Plan.Beams(b).sigmaAtCEF.r);
        fprintf('sigma at base of CEF = %f mm \n',Plan.Beams(b).sigmaAtCEF.sigma);

        if (abs(Plan.Beams(b).spotSigma .* 1.5 - Plan.Beams(b).SpotSpacing) ./ Plan.Beams(b).spotSigma > 0.1 )
          %The spot spacing is 10% off the optimum spot spacing at 1.5 sigma at iso
          Plan.Beams(b).spotSigma
          Plan.Beams(b).SpotSpacing
          warning('Spot spacing is far from the optimum 1.5 * sigma iso')
        end

        if (Plan.Beams(b).sigmaAtCEF.sigma .* 2 <= mag .* Plan.Beams(b).SpotSpacing ./2)
          apothem = Plan.Beams(b).sigmaAtCEF.sigma .* 2; %radius of the base of the spike should be 2 sigma because 95% of the protons go thourhg that area
          fprintf('2 Spot sigma <= spot spacing /2 \n')
        else
          Plan.Beams(b).SpotSpacing;
          apothem = double(mag .* Plan.Beams(b).SpotSpacing ./2);
          fprintf('2 Spot sigma > spot spacing /2 \n')
          warning('Spot fluence will be clipped by neighbouring spike')
        end
        fprintf('Spike base apothem : %f mm \n',apothem);

        %Compute the optimum shape of the CEM
        nrSides = getNbSidesOfPolygon(GridLayout , Plan.Spike.SpikeType);

        [RF , R_max , ~ , ~ , mag] = calculateCEF(Plan.Beams(b) , apothem , Plan.Spike  , spotPosition{b}, GridLayout , Plan.BDL , Plan.ScannerDirectory , nrSides);

        %Define CEM properties
        Plan.Beams(b).RidgeFilter = RF; %Shape of the Conformal Energy Filter + range compensator
        Plan.Beams(b).RangeModulator.CEMThicknessData = max([RF(:).h_step]); %Compute the maximum height of the CEF. This will be used to estimate the isocentre to Range modulator distance
        IsocenterToRangeModulatorDistance  = getIsocenterToRangeModulatorDistance(Plan.Beams(b) , Plan.BDL); %First guess at CEf thickness data
        Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance = IsocenterToRangeModulatorDistance;
        Plan.Beams(b).RangeModulator.RangeModulatorID  = [Plan.CTinfo.OriginalHeader.PatientID ' : ' Plan.Beams(b).name];
        Plan.Beams(b).RangeModulator.AccessoryCode = num2str(now); %GEnerate a unique number based on the present instant

        iso2Nozzle = getNozzle2IsoDistanceFromBDL(Plan.BDL);
        if (Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance > iso2Nozzle)
          fprintf('Isocenter to Modulator Tray Distance must be %f mm \n',Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance)
          fprintf('Maximum distance between isocentre and nozzle frame (in BDL) is %f mm \n',iso2Nozzle)
          fprintf('Select a BDL (= a treatment room) with more space between isocentre and nozzle frame \n')
          fprintf('Alternatively, clear some space by moving the isocentre using the parameter isocenter in YAML \n')
          error('CEM does not fit on the nozzle')
        end


        %Check again spot spacing, knowing the magnificaiton factor for projection to CEF plane
        if (Plan.Beams(b).sigmaAtCEF.sigma .* 2 > min(mag) .* Plan.Beams(b).SpotSpacing ./2)
          Plan.Beams(b).sigmaAtCEF.sigma
          min(mag)
          Plan.Beams(b).SpotSpacing
          fprintf('2 Spot sigma > min(mag) .* spot spacing /2 \n')
          warning('Spot fluence will be clipped by neighbouring spike')
        end

        %Improve the shape design by correcting for the lateral scattering of proton in the spikes
        switch Plan.Spike.Optimization
          case 'convolution'
              %Optimisation on the fluence, using convolution to account for lateral scattering and range straggling
              for b = 1:numel(Plan.Beams)
                Plan = getOptimumSpike(Plan , b); %Update the shape of the spikes for beam b
              end
          case 'none'
              %No optimisation
        end %switch case

        fprintf('Saving temporary files \n')
        save (fullfile(Plan.output_path,['Plan_ridge_tmp',num2str(b),'.mat']),'Plan', '-v7.3')

    else %if isfield
      fprintf('Loading SOBP depth doses and Conformal Energy Filter shape  \n')
      fullfile(Plan.PlanExistentFile,['Plan_ridge_tmp',num2str(b),'.mat'])
      Plan = reloadPartialResults(fullfile(Plan.PlanExistentFile,['Plan_ridge_tmp',num2str(b),'.mat']) , Plan.YAML);
    end %if isfield
end %for b

%Compute the shape of the CEF to be exported in a DICOM plan later
if (Plan.makeSTL)
  %Additionally, save a STL file with the ridge filter
  [CompensatorThicknessData , CompensatorPixelSpacing , ModulatorOrigin , CEM3Dmask, ModulatorMountingPosition] = makeCEF(Plan , true);
else
  %Do not save the STL file
  [CompensatorThicknessData , CompensatorPixelSpacing , ModulatorOrigin , CEM3Dmask, ModulatorMountingPosition] = makeCEF(Plan , false);
end


for b = 1: size(Plan.Beams,2) %Loop for each beam

    %Store CEF data into the structure that will be exported in DICOM plan
    Plan.Beams(b).RangeModulator.CEM3Dmask = CEM3Dmask{b};
    Plan.Beams(b).RangeModulator.CEMThicknessData = CompensatorThicknessData{b};
    Plan.Beams(b).RangeModulator.Modulator3DPixelSpacing = CompensatorPixelSpacing{b};
    Plan.Beams(b).RangeModulator.ModulatorOrigin = ModulatorOrigin{b};
    Plan.Beams(b).RangeModulator.ModulatorMountingPosition = ModulatorMountingPosition;
    Plan.Beams(b).RangeModulator.FrameOfReferenceUID = dicomuid; %UID of the binary file with the CEM mask
    Plan.Beams(b).RangeModulator.ReferencedRangeModulator = dicomuid; %UID of the coordinate system of the CEM
end

end %function


%=========================
%Compute the number of sides of the polygon
%=========================
function nrSides = getNbSidesOfPolygon(GridLayout , SpikeType)

  switch GridLayout
    case 'HEXAGONAL'
      nrSides = 6;
    case 'SQUARE'
      nrSides = 4;
  end
  switch SpikeType
  case {'ellipse','smooth'}
      nrSides = 1;
  end
end
