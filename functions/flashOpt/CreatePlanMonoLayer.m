%% CreatePlanMonoLayer
% Convert a multi energy layer PBS plan into a single energy layer plan (to be used with a hedgehog)
%  * Collect together the weight of all the spots located at the same X,Y positions% and allocate the weight to the deepest energy layer.
%  * Compute the maximum energy of the proton beam
%  * Compute the WET of the range shifter required to shift the max energy down to the energy at entrance of CEF
%  * Describe the hedgehog shape in the Range Modulator Sequence of a DICOM plan. Some private DICOM tag are used in order to define the information not storable using public DICOM tags
%
% |PlanMono| can be saved as a DICOM file using the function |createDICOMPlan|. The DICOM dictionary |plugins\openMIROpt\functions\io\dicom-dict.txt|
% must be used because it constins the deifnition of the priovate DICOM tags for the export ofthe CEF
%
%
%% Syntax
% |PlanMono = CreatePlanMonoLayer(Plan , filename , protonsFullDose)|
%
%
%% Description
% |PlanMono = CreatePlanMonoLayer(Plan , filename , protonsFullDose)| Description
%
%
%% Input arguments
% |Plan| -_STRUCTURE_- Structure defining a multi energy layer PBS plan
%  * |Plan.name| -_STRING_- plan name tag inside the plan file
%  * |Plan.TargetROI_ID| -_SCALAR_- Index of the target volume, as it %   appears in the RTSTRUCT file.
%  * |Plan.fractions| -_SCALAR_- Number of fractions for the treatment.
%  * |Plan.Machine| -_STRING_- Name of the treatment machine. Must be defined in Beam Data Library
%  * |Plan.CTinfo| -_STRUCTURE_- DICOM header of the CT scan
%  * |Plan.output_path| -_STRING_- Folder in which the results are saved
%  * |Plan.TargetROI| -_STRING_- Name of the PTV structure in handles.images
%  * |Plan.CTname| -_STRING_- Name of the CT image in handles.images
%  * |Plan.BeamletsBy| -_STRING_- Algorithm used to compute the beamlets: Monte Carlo ('MCsquare') or Pencil Beam ('FoCa')
%  * |Plan.BDL = Plan.BDL| -_STRING_- Beam data library. Name of the folder in REGGUI\plugins\openMCsquare\lib\BDL
%  * |Plan.ScannerDirectory| -_STRING_- Ct scanner HU to SPR conversion parameters. Name of the folder in REGGUI\plugins\openMCsquare\lib\Scanners
%  * |Plan.Beams| -_VECTOR of STRUCTURES_- Information about the different beams in the plan
%      * |Plan.Beams(i).GantryAngle| -_SCALAR_- Angle (deg) og the i-th beam in the plan
%      * |Plan.Beams(i).PatientSupportAngle| -_SCALAR_- Yaw Angle (deg) of the treatment couch
%      * |Plan.Beams(i).isocenter| -_SCALAR VECTOR_- [x,y,z] Coordiantes (mm) of the isocentre in the planning CT scan for i-th beam
%      * |Plan.Beams(b).Layers(L).Energy| -_SCALAR_- Energy (MeV) of the L-th layer
%      * |Plan.Beams(b).Layers(L).SpotPositions(k)| -_SCALAR VECTOR_- =[x,y] position (mm) of the k-th spot in the layer
%      * |Plan.Beams(b).Layers(L).SpotWeights(k)| -_SCALAR VECTOR_-  Weight per fraction of the k-th spot in the layer
%      * |Plan.Beams(i)RSinfo| -_STRUCTURE_- Definition of the material properties of the range shifter
%         * |Plan.Beams(i).RSinfo.RangeShifterID| -_STRING_- ID of the range shifter as defined in the beam data library
%         * |Plan.Beams(i).RSinfo.RangeShifterType| -_STRING_- At the moment only BINARY is supported.
%         * |Plan.Beams(i).RSinfo.IsocenterToRangeShifterDistance| -_SCALAR_- Distance (in mm) between the downstream side of the range shifter and the isocentre
%      * |Plan.Beams(b).RangeModulator.RangeModulatorID| -_STRING_- (300A,0346) Identifier for Range Modulator. Read from the YAML configuration file
%      * |Plan.Beams(b).RangeModulator.AccessoryCode| -_STRING_- (300A,00F9) An accessory identifier to be read by a device such as a bar code reader. Read from the YAML configuration file
%      * |Plan.Beams(b).RangeModulator.CEMThicknessData| -_SCALAR MATRIX_- |CompensatorThicknessData(x,y)| Thickness (mm) of the CEF pixel at position (x;y) in the IEC beam Limiting device CS
%      * |Plan.Beams(b).RangeModulator.Modulator3DPixelSpacing| -_SCALAR VECTOR_- |CompensatorPixelSpacing = [x,y,z]| Pixel size (mm) in the plane of the CEF for the |CompensatorThicknessData| matrix  in the plane of the CEM
%       * |Plan.Beams(b).RangeModulator.ModulatorOrigin| -_SCALAR VECTOR_- Physical coordinate the voxel |CompensatorThicknessData(1,1)| and |hedgehog3D(1,1,1)| for beam b  in the plane of the CEF.
%      * |Plan.Beams(b).BlockData|  -_SCALAR MATRIX_- |BlockData{h}(i,:)=[x,y]|  Coordinates (mm) of the i-th point defining the contour of the aperture block projected onto the machine isocentric plane in the IEC BEAM LIMITING DEVICE coordinate system for the b-th beam.
%                                     If there are several holes in the block, then |BlockData{h}| descibes the controu of the h-th hole
%   *|Plan.SpotTrajectoryInfo| -_STRUCTURE_- [OPTIONAL: if provided, the spot trajectory is not optimised] Pre-computed spot trajectory. Used to accelerate computations.
%                         |SpotTrajectoryInfo{b}| Parameter for beam |b|
%     * |Plan.SpotTrajectoryInfo.beam{b}.sobpSequence| -_SCALAR VECTOR_-  Order of the indices of |spot| to sort the spots. |OrderedSOBP = spot(sobpSequence,:)|
%     * |Plan.SpotTrajectoryInfo.beam{b}.Nmaps| -_STRUCTURE_-  Topological information on the initial spot sequence
%       * |Nmaps.NeighbourghMap| -_SCALAR VECTOR_- |NeighbourghMap(d,i)| d=# of delivered spot; i=# of impacted spot; |NeighbourghMap(d,i)|  = 0 if there is an impact
%       * |Nmaps.NeighbourghWeightMap| -_SCALAR VECTOR_- |NeighbourghWeightMap(d,i)| fraction of the dose of spot d (# of delivered spot) that is also delivered at spot i (i=# of impacted spot)
%       * |Nmaps.NeighbourghTimeMap| -_SCALAR MATRIX_- |NeighbourghTimeMap(d,i)| Time (ms) required to go from spot d to spot i
%     * |Plan.SpotTrajectoryInfo.sobpPosition{b}| - CELL VECTOR_ - beamletPosition{b}(i,:) = [x,y] The i-th spot of the b-th beam is located at position [x,y] in the BEV coordinates
%     * |Plan.SpotTrajectoryInfo.weight2spot| - SCALAR MATRIX_ - |weight2spot(weightIdx,:) = [b,BeamLetNb,l]| The spot at index |weightIdx| is related to the b-th beam and the SOBP # |BeamLetNb| in layer # l
%
%
% |filename| -_STRING_- File name of the plan to be saved by createDICOMPlan.m
%
% |protonsFullDose| -_SCALAR_- Number of protons to compute the dose in Monte Carlo
%
%
%% Output arguments
%
% |PlanMono| -_STRUCTURE_- Structure defining a single energy layer PBS plan. The weigh of all Bragg peaks
%                 at the same (x,y) position are summed to generate the weight of the SOBP.
%      * |Plan.Beams(b).NumberOfRangeModulators| -_SCALAR_- = 1 (300A,0340) Number of range modulators associated with current beam.
%      * |Plan.Beams(b).RangeModulator.RangeModulatorID| -_STRING_- (300A,0346) Identifier for Range Modulator. Read from the YAML configuration file
%      * |Plan.Beams(b).RangeModulator.AccessoryCode| -_STRING_- (300A,00F9) An accessory identifier to be read by a device such as a bar code reader. Read from the YAML configuration file
%      * |Plan.Beams(b).RangeModulator.RangeModulatorType| -_STRING_- (300A,0348) Type of Range Modulator = 'FIXED'
%      * |Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance| -_SCALAR_- Private DICOM tag (300D,1012). VR type FL. Isocenter to UPSTERAM edge of CEF (mm).
%      * |Plan.Beams(b).RangeModulator.Modulator3DPixelSpacing| -_SCALAR VECTOR_- Private DICOM tag (300D,1015). VR type FL. Physical distance (in mm) between the center of each pixel of |RangeModulator.ModulatorThicknessData| in the plane of the CEM.
%                                         Specified by a numeric pair - adjacent row spacing followed by adjacent column spacing.
%      * |Plan.Beams(b).RangeModulator.ModulatorPosition| -_SCALAR VECTOR_- Private DICOM tag  (300D,1016). VR type FL. The x and y coordinates of the upper left hand corner (first pixel transmitted) of the CEF,
%                                        in the plane of the CEM in the IEC GANTRY coordinate system (mm).
%      * |Plan.Beams(b).RangeModulator.ModulatorThicknessData| -_SCALAR MATRIX_- Private DICOM tag (300D,1010). VR type FL. A data stream of the pixel samples that comprise the range compensator,
%                                         expressed as physical thickness (in mm), parallel to radiation beam axis.
%                                         The order of pixels sent is left to right, top to bottom (upper left pixel, followed by the remainder of row 1, followed by the remainder of the rows).
%       When the plan is exported by |createDICOMPlan|, the following two private tag will be added:
%                 (300D,1013) UL  Number of rows of |RangeModulator.ModulatorThicknessData|
%                 (300D,1014) UL  Number of columns of |RangeModulator.ModulatorThicknessData|
%      * |Plan.Beams(b).NumberOfRangeShifters| -_INTEGER_- = 1 Number of range shifter
%      * |Plan.Beams(b).RSinfo| -_STRUCTURE_- Information about the range shifter
%           * |Plan.Beams(b).RSinfo.RangeShifterID| -_STRING_- ID of the range shifter as defined in the beam data library
%           * |Plan.Beams(b).RSinfo.RangeShifterType| -_STRING_- At the moment only BINARY is supported. Possible types are: BINARY: Device is composed of different thickness materials that can be moved in or out of the beam in various stepped combinations. ANALOG: Device is variable thickness and is composed of opposing sliding wedges, water column or similar mechanism.
%           * |Plan.Beams(b).RSinfo.RangeShifterWET|  -_SCALAR_- Range shifter water equivalent length in mm
%           * |Plan.Beams(b).RSinfo.IsocenterToRangeShifterDistance|  -_SCALAR_- Distance (mm) from isocentre to base of range shifter
%           * |Plan.Beams(b).RSinfo.RangeShifterSetting| -_STRNG_- = 'IN' To specify that the range shifter is IN the beam
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com).


function PlanMono = CreatePlanMonoLayer(Plan , filename , protonsFullDose)

  PlanMono = Plan; %Copy all fields from Plan to PlanMono
  PlanMono = rmfield(PlanMono, 'Beams'); %REmove the beam field. This will be recreated below

  if isfield(Plan , 'scanAlgoGW')
    %The plan contains information about scanAlgo. Let's get hash for scanAlgo as well
    PlanMono.SoftwareVersion = getSWversion(Plan.scanAlgoGW); %REcord the software version in the plan
  else
    %The plan does not contain any information about scanalgo. Record hash for MIROPT and FLASH MIROPT only
    PlanMono.SoftwareVersion = getSWversion(); %REcord the software version in the plan
  end

  %Create a beam with a single energy layer
  for b = 1: size(Plan.Beams,2)

      snout = getParamSnout(Plan.Beams(b).SnoutID);

      %Copy everything
      PlanMono.Beams(b) = Plan.Beams(b);

      %Remove what is irrelevant
      PlanMono.Beams(b).Layers = struct; %

      l=1; %single energy layer

      PlanMono.Beams(b).Layers(l).Energy = max([Plan.Beams(1).Layers(:).Energy]); %The monolayer plan has the max energy of the original plan, which is also the max energy of the machine
      PlanMono.Beams(b).Layers(l).NumberOfPaintings = 1; %Force a single painting of the layer to have a high dose rate

      physicsConstants;
      ChargePerMU = MU_to_NumProtons(1, PlanMono.Beams(b).Layers(l).Energy) .* eV; %Cb per MU
      PlanMono.Beams(b).MetersetRate = 60 .* Plan.Inozzle .* 1e-9 ./ ChargePerMU; % 60 .* Cb/s / (Cb / MU) = 60 .* MU/s = MU/min

      idx = 1;
      Plan = weight2charge(Plan, Plan.BDL , PlanMono.Beams(b).Layers(l).Energy); %convert the weight from the energy of the l-th layer to the **maximum** energy. This account for the MU to charge variation as a function of eproton energy
      w = plan2weight(Plan); %rebuild the vector of optimized weights. These are the weights PER FRACTION
      weight2spot = Plan.SpotTrajectoryInfo.weight2spot; %Get the matrix linking the Bragg peak weights to the SOBP weight

      %Loop sequentially through each spot of the trajectory
      for spt = 1:numel(Plan.SpotTrajectoryInfo.beam{b}.sobpSequence)
          beamletIndex = Plan.SpotTrajectoryInfo.beam{b}.sobpSequence(spt); %Index of the SOBP beamlet in the trajectory sequence
          spotIndices = find((weight2spot(:,1) == b) .* (weight2spot(:,2) == beamletIndex)); %These are all the Bragg peaks contributing to beamlet |beamletIndex| of beam |b|
          weights = w(spotIndices); %REtrieve the weight per fraction of all these Bragg peaks

          if sum(weights) > 0
            %Add the spot only if it has a weight
            PlanMono.Beams(b).Layers(l).SpotPositions(idx,:) = Plan.SpotTrajectoryInfo.sobpPosition{b}(beamletIndex,:); %Spot position (different from spike position because they are not in the same plane)
            PlanMono.Beams(b).Layers(l).SpotWeights(idx) = sum(weights); %Weight of the SOBP = sum of weight of each BP. This is the weight PER FRACTION
            idx = idx + 1;
          end %if sum(
      end %for spt

        %Conformal energy filter
        %-----------------------
        %Define the shape of the CEF and store it in DICOM as a range modulator in a Range Modulator Sequence (300A,0342)
        %The CEF is stored in the Range Compensator sequence of the DICOM plan
        % From a philosophical point of view it is debatable whether it should be a compensaotr or a modulator.
        % Indeed, the CEF is doing both compensator and modulator.
        PlanMono.Beams(b).NumberOfRangeModulators = 1;
        PlanMono.Beams(b).RangeModulator.IBA_ConformalFLASH_energy_modulator = 'IBA ConformalFLASH energy modulator'; %Private Creator identifier
        PlanMono.Beams(b).RangeModulator.RangeModulatorType = '3D_PRINTED';
                                                        %NB: The term '3D_PRINTED' is not in the list of Defined Terms of the DICOM standard.
                                                        %We define a new term as none of these are reflecting the situation of a “variable modulation according to position��?
                                                        %If there is a term commonly or most used for this in the scientific literature, it may be worth to use it. I have seen 3D or 3D-printed for example.

        %Create private DICOM tags in order to store the description of the spikes of the hedhgehog CEF
        PlanMono.Beams(b).RangeModulator.Modulator3DPixelSpacing = Plan.Beams(b).RangeModulator.Modulator3DPixelSpacing;
        PlanMono.Beams(b).RangeModulator.ModulatorOrigin = Plan.Beams(b).RangeModulator.ModulatorOrigin;
        PlanMono.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance = Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance;
        PlanMono.Beams(b).RangeModulator.RangeModulatorID = Plan.Beams(b).RangeModulator.RangeModulatorID; %GEnerate a unique number based on the present instant
        PlanMono.Beams(b).RangeModulator.AccessoryCode = Plan.Beams(b).RangeModulator.AccessoryCode;
        PlanMono.Beams(b).RangeModulator.ModulatorThicknessData = Plan.Beams(b).RangeModulator.CEMThicknessData;
        PlanMono.Beams(b).RangeModulator.ModulatorMountingPosition = Plan.Beams(b).RangeModulator.ModulatorMountingPosition;
        PlanMono.Beams(b).RangeModulator.ModulatorMaterialID = Plan.Spike.MaterialID;

        %Range shifter
        %-----------------
        %Add the range shifter to degrade the maximum energy down to the energy at entrance of CEF
        if isfield(Plan.Beams(b),'RSinfo')

              PlanMono.Beams(b).NumberOfRangeShifters = 1;
              param = getMachineParam(Plan.BDL);
              PlanMono.Beams(b).RSinfo.RangeShifterNumber = 1;
              PlanMono.Beams(b).RSinfo.RangeShifterID = snout.AccessoryCode(sum(Plan.Beams(b).RSinfo.RSslabThickness) ./ snout.RSslabThickness(end));
              PlanMono.Beams(b).RSinfo.RSslabThickness = Plan.Beams(b).RSinfo.RSslabThickness; %mm

              %PlanMono.Beams(b).RSinfo.RangeShifterMaterial = Plan.Beams(b).RSinfo.RangeShifterMaterial;
              descrip = ['RangeShifterMaterial = ''' convertStringsToChars(Plan.Beams(b).RSinfo.RangeShifterMaterial) ''' ;  RSslabThickness = [' num2str(Plan.Beams(b).RSinfo.RSslabThickness') , '] ; snoutType = ''' param.snout.snoutType '''; '];
              PlanMono.Beams(b).RSinfo.RangeShifterDescription = descrip; %User defined description of Range Shifter.
        else
          PlanMono.Beams(b).NumberOfRangeShifters = 0; %There is no rnage shifter
        end

        %Aperture block
        %--------------
        if (Plan.Beams(b).ApertureBlock)
            %If there are several blocks, this corresponds to several holes in the same aperture block
            %Each controu of hole is represented in a different 'IonBlockSequence' but they are all located at the
            %same distance from the isocentre
            for BlckNb = 1:numel(Plan.Beams(b).BlockData)
                PlanMono.Beams(b).IonBlockSequence{BlckNb}.BlockMountingPosition = Plan.Beams(b).BlockMountingPosition;
                PlanMono.Beams(b).IonBlockSequence{BlckNb}.IsocenterToBlockTrayDistance = Plan.Beams(b).IsocenterToBlockTrayDistance; %The aperture is placed downstream to the ridge filter and the range compensator
                PlanMono.Beams(b).IonBlockSequence{BlckNb}.BlockType = 'APERTURE';
                PlanMono.Beams(b).IonBlockSequence{BlckNb}.BlockDivergence = 'ABSENT';
                PlanMono.Beams(b).IonBlockSequence{BlckNb}.BlockNumber = BlckNb;
                PlanMono.Beams(b).IonBlockSequence{BlckNb}.BlockName = Plan.Beams(b).name; %from yaml Name of the beam
                PlanMono.Beams(b).IonBlockSequence{BlckNb}.BlockTrayID = Plan.Beams(b).name; %from yaml  Name of the beam
                PlanMono.Beams(b).IonBlockSequence{BlckNb}.MaterialID = Plan.Beams(b).BlockMaterialID;
                PlanMono.Beams(b).IonBlockSequence{BlckNb}.BlockThickness = Plan.Beams(b).BlockThickness;
                PlanMono.Beams(b).IonBlockSequence{BlckNb}.BlockNumberOfPoints = size(Plan.Beams(b).BlockData{BlckNb},1);
                PlanMono.Beams(b).IonBlockSequence{BlckNb}.BlockData = Plan.Beams(b).BlockData{BlckNb};
              end
        end
  end

  %Update the trajectory information with the monolayer spot ordering
  %The spot ordering has been sorted out by the CreateMonoLayerPLan function.
  [sobpPosition , weight2spot ] =  collectSpotsinBeamlets(PlanMono); %This is a monolayer plan but call collectSpotsinBeamlets to nicely package the data in proper structures
  PlanMono.SpotTrajectoryInfo.sobpPosition =  sobpPosition;
  PlanMono.SpotTrajectoryInfo.weight2spot = weight2spot;
  PlanMono.SpotTrajectoryInfo.beam{1}.sobpSequence  = weight2spot(:,2);

end
