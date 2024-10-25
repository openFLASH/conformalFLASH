%% parseFLASHplan
% Load the DICOM PT plan for FLASH treatment
% This is a Monolayer FLASH plan with a CEM
%
%% Syntax
% |[handles, Plan] = parseFLASHplanMono(planFileName , CEMheaderFile , CEMdataFile, Plan, handles)|
%
%
%% Description
% |[handles, Plan] = parseFLASHplanMono(planFileName , []  , [] , Plan, handles)| Read the plan. The CEM is in a private tag in plan
%
% |[handles, Plan] = parseFLASHplanMono(planFileName , [] , CEMdataFile, Plan, handles)| Read the plan. The CEM is is at Erik format in a text file
%
% |[handles, Plan] = parseFLASHplanMono(planFileName , CEMheaderFile , CEMdataFile, Plan, handles)| Read the plan. The CEM is is at Claes format in two text files
%
%
%% Input arguments
% |planFileName| -_STRING_- Full path and file name to the RT DICOM plan
%
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are stored.
%
% |handles| - _STRUCT_ - REGGUI data structure.
%
%% Output arguments
%
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are stored.
%
% |handles| - _STRUCT_ - REGGUI data structure.
%
%
%% Contributors
% Authors : R. Labarbe, L. Hotoiu (open.reggui@gmail.com)

function [handles, Plan] = parseFLASHplan(planFileName , Plan, handles)

    if nargin < 2
        %Create the Plan structure
        Plan = struct;
        [~,Plan.output_path] = get_reggui_path();
        Plan.showGraph = false;
    end

    if nargin < 3
        %Create the handles structure
        handles = struct;
        handles.path = Plan.output_path;
        handles = Initialize_reggui_handles(handles);
    end

    [planDir,planFile,EXT] = fileparts(planFileName);
    handles.path = planDir;
    handles.dataPath = planDir;
    name = 'FLASHplan';

    %Load the plan from disk. Keep all the DICOM tags
    setFlashDICOMdict(); %If not already defined, load the DICOM dictionary with private FLASH tags
    monoPlan = dicominfo(planFileName);

    %Check scanAlgo software version
    SWvPlan = split(monoPlan.SoftwareVersion,'\'); %Get the list of softare version
    foundFlag = false;
    if isfield(Plan , 'scanAlgoGW')
      tmp = getSWversion(Plan.scanAlgoGW);
      tmp = split(tmp{3},':');
      SWvHere = strtrim(tmp{2});
    else
      SWvHere = []; %We did not recieve information about scanAlgo
    end

    foundFlag = false;

    for idx = 1:numel(SWvPlan)
       if ~isempty(strfind(SWvPlan{idx} , 'ScanAlgo'))
          foundFlag = true;
          data = split(SWvPlan{idx},':');
          if strcmp(strtrim(data{2}) , SWvHere)
            fprintf('Checked scanning controller version : match')
          else
            fprintf('ScanAlgo in plan : %s \n', strtrim(data{2}) )
            fprintf('ScanAlgo here    : %s \n', SWvHere )
            warning('Mismatch in scanAlgo version between plan and current installtion')
          end
       end
    end

    if ~foundFlag
      warning('Plan did not contain scanAlgo version')
    end


    %Load the plan from disk. Use a REGGUI function that will do some processing ofthe data to put it at the REGGUI format
    handles = Import_plan(planDir, [planFile EXT], 1, name, handles); %Import_plan only import TreatmentDeliveryType = TREATMENT and ignores SETUP beams.
    data = Get_reggui_data(handles,name,'plans');

    %Analyse the content of the FLASH plan and convert into a MIROPT plan
    NbBeams = numel(data);

    for b = 1:NbBeams
      %Create a list with the beam name as read by Reggui importer
      BeamNameListReggui{b} = data{b}.name;
    end

    for b = 1:length(fieldnames(monoPlan.IonBeamSequence))
      %Create list of beam name and beam number as read by DICOM importer
      itemBeam = sprintf('Item_%i',b);
      BeamNumberList(b) = monoPlan.IonBeamSequence.(itemBeam).BeamNumber;
      BeamNameListPlan{b} =  monoPlan.IonBeamSequence.(itemBeam).BeamName;
    end

    %Determine the number of fractions in which each beam is delivered
    %-----------------------------------------------------------------------
    %The ConformalFLASH beams are delivered as one beam per fraction.
    %Plan.fractions(b) : Number of fractions for the b-th beam. The spot weights vector for the beam define the weight per fraction
    %
    %This is encoded in the plan in one of the following methods:
    % METHOD 1:
    % .........
    % Each beam is duplicated in IonBeamSequence the the number of times it will be delivered. The spot weight of each duplicated beam is divided (at the TPS)
    % by the number of duplications so that the dose adds up to the correct value
    % "FractionSequenceGroup" contains one element which references all the beams to be delivered. "NumberOfFractionPLanned" for this sequence is equal to 1
    % We just need to compute the dose for each beam, with a fraction of 1 each time.
    Plan.fractions = zeros(1,NbBeams);
    for FracG = 1:length(fieldnames(monoPlan.FractionGroupSequence))
      %Loop for every fraction defined in the plan
      itemFrac = sprintf('Item_%i',FracG);
      NbFraction = monoPlan.FractionGroupSequence.(itemFrac).NumberOfFractionsPlanned; %Number of fraction for this item
      BeamSequence = monoPlan.FractionGroupSequence.(itemFrac).ReferencedBeamSequence;
      %Loop for every beam reference in the fraction
      for rb = 1:length(fieldnames(BeamSequence))
        itemRefBeam = sprintf('Item_%i',rb);
        RefBeamNb = monoPlan.FractionGroupSequence.(itemFrac).ReferencedBeamSequence.(itemRefBeam).ReferencedBeamNumber; %Beam number reference in DICOM plan
        b = find(BeamNumberList == RefBeamNb); %Find the index in the list that has that beam number
        if ~isempty(b)
          %Add the number of fraction only if this is a treatment beam that was loaded by Reggui
          bReggui = strmatch(BeamNameListPlan{b} , BeamNameListReggui ); %Find the beam that has the same name and was loaded by REGGUI
          Plan.fractions(bReggui) = Plan.fractions(bReggui) + NbFraction; %Plan.fractions(b) Number of fractions for the b-th beam. The final spot weights vector for that beam will be divided by this number
        end %if ~isempty(b)
      end %for rb
    end %for FracG

    Plan.name = monoPlan.RTPlanLabel;
    Plan.FileName = 'Plan'; %default value for the file name of the plan

    %Construct the beam structure
    %------------------------------
    for b = 1:NbBeams

        itemBeam = sprintf('Item_%i',b);
        NbLayers = numel(data{b}.spots);
        Layers = data{b}.spots;
        spot = data{b}.spots.xy;

        MachineName = monoPlan.IonBeamSequence.(itemBeam).TreatmentMachineName;

        fprintf('BDL name      : %s \n' , Plan.BDL);
        [Plan.MachineType , Plan.Machine.name] = getMachineFromBDL(Plan.BDL); %Read machine type and machine name from BDL
        fprintf('Machine name  : %s \n' , Plan.Machine.name)
        fprintf('Machine type  : %s \n' , Plan.MachineType)

        Plan.Beams(b).GantryAngle = data{b}.gantry_angle;
        Plan.Beams(b).PatientSupportAngle = data{b}.table_angle;
        Plan.Beams(b).isocenter = data{b}.isocenter;
        Plan.Beams(b).name = data{b}.name;
        Plan.Beams(b).VDSA = monoPlan.IonBeamSequence.(itemBeam).VirtualSourceAxisDistances;

        Plan.Beams(b).CEFbaseWET = 0 ; %There is no residual range shift in CEM base

        if ( ~isempty(Plan) && isfield(Plan, 'Extras') )
            if isfield(Plan.Extras, 'NbScarves')
                Plan.Beams(b).NbScarves = Plan.Extras.NbScarves;
            end
        end

        %Create the beam structure in the Plan
        %--------------------------------------
        for e = 1:NbLayers
            Plan.Beams(b).Layers(e).Energy = Layers(e).energy;
            Plan.Beams(b).Layers(e).nominalSpotPosition = Layers(e).xy;
            Plan.Beams(b).Layers(e).SpotPositions = Layers(e).xy;
            Plan.Beams(b).Layers(e).SpotWeights = (Layers(e).weight)' ; %This is understood as the weight PER fraction by MIROPT
                          % If the BDL is different in MIROPT and RayStation some rescaling of the MU definition will be required using the doseMeterSet tag
            minW = min(Layers(e).weight);
            maxW = max(Layers(e).weight);
            if Plan.showGraph
                figure(100+b)
                scatter(Layers(e).xy(:,1),Layers(e).xy(:,2) , 50 , round(255.*(Layers(e).weight-minW) ./ (maxW-minW)) , 'filled')
                hold on
            end
        end

        if Plan.showGraph
            hcb = colorbar;
            set(get(hcb,'Title'),'String','Spot charge (AU)')

            figure(100+b)
            grid on
            title(['Spot grid for beam ' num2str(b) '@ isocenter'])
            xlabel('X (mm)')
            ylabel('Y (mm)')
            hold off
            drawnow
        end

        physicsConstants;
        maxE =max([Plan.Beams(b).Layers(:).Energy]);
        ChargePerMU = MU_to_NumProtons(1, maxE) .* eV; %Cb per MU
        if ~isfield(Plan , 'Inozzle')
            if isfield(monoPlan.IonBeamSequence.(itemBeam).IonControlPointSequence.Item_1 , 'MetersetRate')
              fprintf('Loading Inozzle from DICOM plan \n')
              Plan.Inozzle = monoPlan.IonBeamSequence.(itemBeam).IonControlPointSequence.Item_1.MetersetRate .* ChargePerMU .* 1e9 ./ 60; %Proton beam current (nA)
            else
              warning('Plan is missing MetersetRate. Using default value.')
              Plan.Inozzle = 500; %nA
            end
          else
            fprintf('Inozzle already defined in variable Plan \n')
          end
        fprintf('Proton beam current : %f nA \n', Plan.Inozzle)


        %Read the aperture data
        %-------------------------------
        if isfield(monoPlan.IonBeamSequence.(itemBeam) , 'IonBlockSequence')
            %There is an aperture block
            BlockSeq = fieldnames(monoPlan.IonBeamSequence.(itemBeam).IonBlockSequence);
            Plan.Beams(b).ApertureBlock = 1; %set the flag for the aperture block

            for BlckNb = 1:numel(BlockSeq)
              block = getfield(monoPlan.IonBeamSequence.(itemBeam).IonBlockSequence, ['Item_' num2str(BlckNb)]);
              Plan.Beams(b).BlockMountingPosition = block.BlockMountingPosition;
              Plan.Beams(b).BlockMaterialID = block.MaterialID;
              Plan.Beams(b).BlockThickness = block.BlockThickness;

              switch block.BlockMountingPosition
                case 'SOURCE_SIDE'
                    %The distance is defined on the downstream side of the aperture.
                    %In MIROPT, IsocenterToBlockTrayDistance is defined PATIENT_SIDE (i.e. the upstream side)
                    %Need to adjust data from plan: increase the distance
                    Plan.Beams(b).IsocenterToBlockTrayDistance = block.IsocenterToBlockTrayDistance + Plan.Beams(b).BlockThickness;
                case 'PATIENT_SIDE'
                    warning('Plan defines aperture on PATIENT_SIDE')
                    Plan.Beams(b).IsocenterToBlockTrayDistance = block.IsocenterToBlockTrayDistance;
              end

              dataBlock = reshape(block.BlockData,2,block.BlockNumberOfPoints);
              Plan.Beams(b).BlockData{BlckNb} = dataBlock';

              if Plan.showGraph
                  figure(100+b)
                  hold on
                  plot(dataBlock(1,:) , dataBlock(2,:) , '-r')
                  hold off
                  drawnow
              end
            end
        end

        %Get snout information
        %---------------------
        Plan.Beams(b).SnoutID = monoPlan.IonBeamSequence.(itemBeam).SnoutSequence.Item_1.SnoutID;
        % if ~strcmp(Plan.Beams(b).SnoutID , 'FLASH_Snout_S')
        %   fprintf('SnoutID in the plan : %s \n',Plan.Beams(b).SnoutID)
        %   warning('This is not a FLASH snout. Overwriting snout ID')
        %   Plan.Beams(b).SnoutID = 'FLASH_Snout_S';
        % end
        %The plan defines the snout position on the UPSTREAM side of the aperture block
        Plan.Beams(b).SnoutPosition = monoPlan.IonBeamSequence.(itemBeam).IonControlPointSequence.Item_1.SnoutPosition;

        %Define range shifter properties
        %-------------------------------
        Plan.Beams(b).NumberOfRangeShifters = monoPlan.IonBeamSequence.(itemBeam).NumberOfRangeShifters;
        if Plan.Beams(b).NumberOfRangeShifters
              %There is a range shifter
              Plan.Beams(b).RSinfo = monoPlan.IonBeamSequence.(itemBeam).RangeShifterSequence.Item_1;
              snout = getParamSnout(Plan.Beams(b).SnoutID);
              Plan.Beams(b).RSinfo.RSslabThickness = snout.RSslabThickness(snout.RangeShifterSlabs(Plan.Beams(b).RSinfo.RangeShifterID));
              Plan.Beams(b).RSinfo.NbSlabs = numel(find(Plan.Beams(b).RSinfo.RSslabThickness));
              Plan.Beams(b).RSinfo.SlabOffset = snout.RangeShifterOffset(1:Plan.Beams(b).RSinfo.NbSlabs) - snout.RangeShifterOffset(1) + Plan.Beams(b).RSinfo.RSslabThickness(1) ; %Offset from |IsocenterToRangeShifterDistance| and the upstream side of the i-th slab
              fprintf('Range shifter thickness : %f mm \n', Plan.Beams(b).RSinfo.RSslabThickness)
              fprintf('Number of slabs : %d \n', Plan.Beams(b).RSinfo.NbSlabs)

              Plan.Beams(b).RSinfo.IsocenterToRangeShifterDistance = monoPlan.IonBeamSequence.(itemBeam).IonControlPointSequence.Item_1.RangeShifterSettingsSequence.Item_1.IsocenterToRangeShifterDistance; %Distance from isocenter to downstream surface of range shifter
              Plan.Beams(b).RSinfo.RangeShifterSetting = monoPlan.IonBeamSequence.(itemBeam).IonControlPointSequence.Item_1.RangeShifterSettingsSequence.Item_1.RangeShifterSetting;
              Plan.Beams(b).RSinfo.RangeShifterMaterial = snout.RangeShifterMaterial;
              fprintf('Range shifter material : %s \n', Plan.Beams(b).RSinfo.RangeShifterMaterial)


              if ~isfield(monoPlan.IonBeamSequence.(itemBeam).IonControlPointSequence.Item_1.RangeShifterSettingsSequence.Item_1, 'RangeShifterWaterEquivalentThickness')
                %The range shifter WET is not defined in the plan. Compute it
                water = materialDescription('water');
                [~, ~, SPRrs] =  getMaterialPropCT(Plan.Beams(b).RSinfo.RangeShifterMaterial, Plan.ScannerDirectory); %CEM relative stopping power
                WET_RS = Plan.Beams(b).RSinfo.RSslabThickness .* SPRrs;  %Water equivalent thickness (mm) of the range shifter
                Plan.Beams(b).RSinfo.RangeShifterWET = sum(WET_RS); %Range shifter WET in mm
              else
                %the range shifter WET is in the plan. Just copy it
                Plan.Beams(b).RSinfo.RangeShifterWET = double(monoPlan.IonBeamSequence.(itemBeam).IonControlPointSequence.Item_1.RangeShifterSettingsSequence.Item_1.RangeShifterWaterEquivalentThickness); %Range shifter WET in mm
              end
              if isfield(monoPlan.IonBeamSequence.(itemBeam).IonControlPointSequence.Item_1.RangeShifterSettingsSequence.Item_1 , 'IsocenterToRangeShifterDistance')
                  if round(double(Plan.Beams(b).RSinfo.IsocenterToRangeShifterDistance),1) ~= round(double(monoPlan.IonBeamSequence.(itemBeam).IonControlPointSequence.Item_1.RangeShifterSettingsSequence.Item_1.IsocenterToRangeShifterDistance),1)
                    fprintf('Isocenter To RangeShifter Distance from snout position        : %f mm \n', round(double(Plan.Beams(b).RSinfo.IsocenterToRangeShifterDistance),1))
                    fprintf('Isocenter To RangeShifter Distance from DICOM tag (300A,0364) : %f mm \n', round(double(monoPlan.IonBeamSequence.(itemBeam).IonControlPointSequence.Item_1.RangeShifterSettingsSequence.Item_1.IsocenterToRangeShifterDistance),1))
                    warning('Isocenter To RangeShifter Distance inconsistent with FLASH snout')
                  end
              else
                warning('DICOM tag (300A,0364) Isocenter To RangeShifter Distance is missing')
              end
        end

        %Define parameters for CEM
        %--------------------------

        %The hedgehog is defined in a private tag in the plan
        itemCEM = sprintf('Item_%i',1);
        Plan.Beams(b).NumberOfRidgeFilters = monoPlan.IonBeamSequence.(itemBeam).NumberOfRangeModulators;
        if Plan.Beams(b).NumberOfRidgeFilters
            %There is a range modulator
            Plan.RidgeFilter = true;
            Plan.Beams(b).RangeModulator.AccessoryCode = monoPlan.IonBeamSequence.(itemBeam).RangeModulatorSequence.(itemCEM).AccessoryCode;
            Plan.Beams(b).RangeModulator.RangeModulatorID = monoPlan.IonBeamSequence.(itemBeam).RangeModulatorSequence.(itemCEM).RangeModulatorID;
            Plan.Beams(b).RangeModulator.RangeModulatorType = monoPlan.IonBeamSequence.(itemBeam).RangeModulatorSequence.(itemCEM).RangeModulatorType;
            if ~strcmp(Plan.Beams(b).RangeModulator.RangeModulatorType , '3D_PRINTED') & ~strcmp(Plan.Beams(b).RangeModulator.RangeModulatorType , 'FIXED')
              %In the first software reelase, the ConformalFLASH plan plans were exported with '3D_PRINTED'
              %In the new software release, the plans are exported with 'FIXED' because '3D_PRINTED' is not in the list of Defined Terms of the DICOM standard and Aria reject it
              Plan.Beams(b).RangeModulator.RangeModulatorType
              error('Wrong type of ConformalFLASH energy modulator')
            end

            Plan.Spike.MaterialID = remove_bad_chars(getPrivateTag('300D' , '0018' , 'IBA'  , monoPlan.IonBeamSequence.(itemBeam).RangeModulatorSequence.(itemCEM) , 'ModulatorMaterialID'));
            Plan.Spike.MaterialID = strip(Plan.Spike.MaterialID);
            if (strcmp(Plan.Spike.MaterialID(end),"_"))
                Plan.Spike.MaterialID = Plan.Spike.MaterialID(1:end-1);
            end
            Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance = Plan.Beams(b).SnoutPosition + snout.CEMOffset ; %| -_SCALAR_- Distance (mm) from isocentre to the base of the CEF.
            fprintf('Distance isocenter to base of CEM : %3.1f mm \n',Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance)

            nrPixelsY = double(getPrivateTag('300D' , '0013' , 'IBA'  ,monoPlan.IonBeamSequence.(itemBeam).RangeModulatorSequence.(itemCEM) , 'ModulatorRows'));
            nrPixelsX = double(getPrivateTag('300D' , '0014' , 'IBA'  ,monoPlan.IonBeamSequence.(itemBeam).RangeModulatorSequence.(itemCEM) , 'ModulatorColumns'));

            ModulatorPixelSpacing =  getPrivateTag('300D' , '0015' , 'IBA'  ,monoPlan.IonBeamSequence.(itemBeam).RangeModulatorSequence.(itemCEM) , 'ModulatorPixelSpacing'); %pixel spacing in the plane of the CEM
            ModulatorPixelSpacing = flip(ModulatorPixelSpacing,2); %Stored in DICOM as [Y,X]

            %Define Z resolution: this is the smallest dZ step between two terraces of the tower
            z_unif_val = 1; % we consider 1 mm as uniform step value in tower
            Modulator3DPixelSpacing = round(double([ModulatorPixelSpacing', z_unif_val]),1); %| -_SCALAR VECTOR_- |CompensatorPixelSpacing = [x,y,z]| Pixel size (mm) in the plane of the CEF for the |CompensatorThicknessData| matrix in the plane of the CEM
            fprintf('CEM pixel size from Dicom plan : (%3.1f, %3.1f, %3.1f) mm \n',Modulator3DPixelSpacing(1),Modulator3DPixelSpacing(2),Modulator3DPixelSpacing(3))

            if(Modulator3DPixelSpacing(1) ~= Modulator3DPixelSpacing(2))
              error('Pixels of the elevation map are not square')
            end
            Plan.Beams(b).RangeModulator.Modulator3DPixelSpacing = Modulator3DPixelSpacing;

            %Convert the 3D elevation map from DICOM file into a 3D mask.
            %elvMap2mask takes care of the flip of the Y axis
            [CEM3Dmask1 , CEMThicknessData] = elvMap2mask(double( getPrivateTag('300D' , '0020' , 'IBA', monoPlan.IonBeamSequence.(itemBeam).RangeModulatorSequence.(itemCEM) , 'ModulatorThicknessData') ) , nrPixelsX , nrPixelsY , Modulator3DPixelSpacing);

            ModulatorPosition = getPrivateTag('300D' , '0016' , 'IBA'  ,monoPlan.IonBeamSequence.(itemBeam).RangeModulatorSequence.(itemCEM) , 'ModulatorPosition');
            Plan.Beams(b).RangeModulator.ModulatorOrigin = [ModulatorPosition' , 0]; %| -_SCALAR VECTOR_- Physical coordinate [x,y,z] the voxel |CompensatorThicknessData(1,1)| and |hedgehog3D(1,1,1)| for beam b  in the plane of the CEF.
            Plan.Beams(b).RangeModulator.ModulatorOrigin(2) = Plan.Beams(b).RangeModulator.ModulatorOrigin(2) - (nrPixelsY-1) .* Modulator3DPixelSpacing(2); %We will flip the Y axis, so change sign of origin

            if isfield(monoPlan.IonBeamSequence.(itemBeam).RangeModulatorSequence.(itemCEM), 'ModulatorMountingPosition')
              Plan.Beams(b).RangeModulator.ModulatorMountingPosition = getPrivateTag('300D' , '0030' , 'IBA'  ,monoPlan.IonBeamSequence.(itemBeam).RangeModulatorSequence.(itemCEM), 'ModulatorMountingPosition');
              if ~strcmp(getPrivateTag('300D' , '0030' , 'IBA'  ,monoPlan.IonBeamSequence.(itemBeam).RangeModulatorSequence.(itemCEM), 'ModulatorMountingPosition') , 'PATIENT_SIDE')
                error('ModulatorMountingPosition must be PATIENT_SIDE')
              end
            else
              warning('ModulatorMountingPosition not in plan. Using PATIENT_SIDE as default')
              Plan.Beams(b).RangeModulator.ModulatorMountingPosition = 'PATIENT_SIDE';
            end

            Plan.Beams(b).RangeModulator.CEMThicknessData = CEMThicknessData;  %| -_SCALAR MATRIX_- |CompensatorThicknessData(x,y)| Thickness (mm) of the CEF pixel at position (x;y) in the IEC beam Limiting device CS
            Plan.Beams(b).RangeModulator.CEM3Dmask = CEM3Dmask1; % | -_SCALAR MATRIX_- 3D mask of the CEF. |CEM3Dmask(x,y,z)=1| if the voxel at location (x,y,z)  in the plane of the CEF for beam b belongs to the CEF.
                                                    % Z=0 at the base of CEF. Z increase in the smae way as Zg if the spike point toward the proton source

            %Display the elevation map of the CEM
            X = 1:nrPixelsX;
            Y = 1:nrPixelsY;
            ModulatorOrigin(1) = Plan.Beams(b).RangeModulator.ModulatorOrigin(1);
            ModulatorOrigin(2) = Plan.Beams(b).RangeModulator.ModulatorOrigin(2);
            X = (X-1) .* ModulatorPixelSpacing(1) + ModulatorOrigin(1);
            Y = (Y-1) .* ModulatorPixelSpacing(2) + ModulatorOrigin(2);
            CEMcontourPlot(48+b*2 , X , Y, CEMThicknessData , Plan.Beams(b).BlockData, Plan.Beams(b).VDSA , Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance);


            if isfield(monoPlan.IonBeamSequence.(itemBeam).IonControlPointSequence.Item_1, 'RangeModulatorSettingsSequence')
                if isfield(monoPlan.IonBeamSequence.(itemBeam).IonControlPointSequence.Item_1.RangeModulatorSettingsSequence.Item_1 , 'IsocenterToRangeModulatorDistance')
                    if round(double(Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance),1) ~= round(double(monoPlan.IonBeamSequence.(itemBeam).IonControlPointSequence.Item_1.RangeModulatorSettingsSequence.Item_1.IsocenterToRangeModulatorDistance),1)
                      fprintf('Isocenter to Modulator Tray Distance from snout position        : %f mm \n', round(double(Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance),1))
                      fprintf('Isocenter to Modulator Tray Distance from DICOM tag (300D,1012) : %f mm \n', round(double(monoPlan.IonBeamSequence.(itemBeam).IonControlPointSequence.Item_1.RangeModulatorSettingsSequence.Item_1.IsocenterToRangeModulatorDistance),1))
                      warning('CEM position inconsistent with FLASH snout')
                    end
                else
                    warning('DICOM tag (300D,1012) Isocenter to Modulator TrayDistance is missing')
                end
            end
      end

      Plan.Beams(b).spotSigma = 10; %mm It is only used to determine neighbourgh spots. The exact value is not too critical

    end %for b
end
