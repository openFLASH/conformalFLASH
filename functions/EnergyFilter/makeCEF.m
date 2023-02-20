%% makeCEF
% Compute a 2D elevation map and a 3D mask describing the CEM. The masks have the spatial resolution of |Plan.Spike.intrpCTpxlSize|
%
% In addition, export a STL file for 3D printer with the shape of the Conformal Energy Filter.
% There will be one file per beam in the plan.
% The STL file describes the shape of the spikes, the vertical sides of the block and the bottom surface of the CEF. This defines a solid object for a 3D printer.
% The facet coordinates are in MILLIMETERS. % NB: The STL format does not define a length unit [3].
% In the STL file, the X,Y dimension of the CEF are demagnified so that the CEF is projected at the position of the accessory drawer
%
% The file name will be |Plan.Beams().RangeModulator.AccessoryCode| with the '.stl' extension
% If empty, the function return the description of the CEF height in |CEMelevation|
%
%% Syntax
% |[CEM2DElvMap , pixelSize , origin , CEM3Dmask , ModulatorMountingPosition] = makeCEF(Plan , SaveSTL)|
%
%
%% Description
% |[CEM2DElvMap , pixelSize , origin , CEM3Dmask , ModulatorMountingPosition] = makeCEF(Plan , SaveSTL)| Description
%
%
%% Input arguments
% |Plan| - _struct_ - MIROpt structure with updated information:
%   * |Plan.bevPTV| - _CELL VECTOR OF SCALAR MATRIX_ - |proj(x,y)| Value of the paralell projection of the |target| on a a surface perpendicular to the |beam| axis. The surface is centered on |beam| axis. It is a square with size equal to |bev_size|*|bev_size|
%   * |Plan.bev_x| -_CELL VECTOR OF SCLAR VECTOR_- |Plan.bev_x{b}(i)| X coordinate of Plan.bevPTV{b}(i,:) or Y coordinate of Plan.bevPTV{b}(:,i) for b-th beam
%   * |Plan.Spike.intrpCTpxlSize| -_SCALAR_- Size (mm) of the voxels in the interpolated CT used to compute dose through CEF. The CEF mask will be generated with this resolution
%   * |Plan.Spike.MinThickness| -_SCALAR_- Thickness (mm) in of the base on which the spikes are built. The base has the same |R_WET| as the spikes
%   * |Plan.Spike.SpikeType| -_STRING_- Type of spike to be designed. The centre of the spike corresponds to the BP with smaller range ('up') or the largest range ('down') or randomise pixel column ('random'), apply Gaussian filter to "smear"('smooth'), draw elliptical spike ('ellipse')
%   * |Plan.Beams(b).GridLayout| -_STRING_- Layout of the PBS spots on the grid. Options: HEXAGONAL (default), SQUARE
%   * |Plan.Beams(b).RidgeFilter| -_VECTOR of STRUCT_- Structure describing the shape of the Conformal Energy Filter filer
%     * |Plan.Beams(b).RidgeFilter(k).x_centre| -_SCALAR_- X coordinate (mm in IEC Gantry CS) of central axis of the k-th spike projected in the plane of the CEF
%     * |Plan.Beams(b).RidgeFilter(k).y_centre| -_SCALAR_- Y coordinate (mm in IEC Gantry CS) of central axis of the k-th spike projected in the plane of the CEF
%     * |Plan.Beams(b).RidgeFilter(k).h_outside| -_SCALAR_- Height (mm) of the base of the k-th spike
%     * |Plan.Beams(b).RidgeFilter(k).a_max| -_SCALAR VECTOR_- |RF(k).a_max(L)| Outer diameter (mm) of the L-th ring (= external diameter of the L-th polygon/square)
%     * |Plan.Beams(b).RidgeFilter(k).a_min| -_SCALAR VECTOR_- |RF(k).a_min(L)| Inner diameter (mm) of the L-th ring (= external diameter of the (L+1)-th polygon/square)
%     * |Plan.Beams(b).RidgeFilter(k).h_step| -_SCALAR VECTOR_- |RF(k).h_step(L)| Height (mm) of the L-th ring
%
% |SaveSTL| -_BOOL_- TRUE : Save the STL file to disk
%
%% Output arguments
%
% |CEM2DElvMap| -_SCALAR MATRIX_- |CEM2DElvMap{b}(x,y)| Physical height (mm) of the CEF at the position (x,y) (Gantry CS)  in the plane of the CEF for beam b. Voxel ordered in increasing coordinates
%
% |pixelSize| -_SCALAR MATRIX_- |pixelSize{b} = [x,y,z]| Pixel size of the CEF in |CEMelevation| for beam b. The CEF is specified in the plane of the CEF
%
% |origin|  -_SCALAR MATRIX_- |origin{b}= [x,y,z]| Coordinate (mm) in CEM CS of the the voxel |CEM2DElvMap(1,1)| and |CEM3Dmask(1,1,1)| for beam b (in the plane of the CEF).
%
% |CEM3Dmask|  -_SCALAR MATRIX_- 3D mask of the CEF. |CEM3Dmask(x,y,z)=1| if the voxel at location (x,y,z)  in the plane of the CEF for beam b belongs to the CEF.
%                                     The z index is 1 at the base of the CEM and the z index increses going up the spikes. x and y indices are aligned with the IEC gantry CS
%
% |ModulatorMountingPosition| -_STRING_- DICOM tag Modulator Mounting Position (300D,1017). Defines the face of the CEM used as reference for |IsocenterToRangeModulatorDistance|
%                         SOURCE_SIDE is using the downstream face of the block as a reference position for expressing the isocenter to block tray distance.
%                         PATIENT_SIDE is using the upstream face
%
%% Contributors
% Authors : R. Labarbe, Lucian Hotoiu (open.reggui@gmail.com)
%
% REFERENCES
% [1] https://www.wikihow.com/Calculate-the-Area-of-a-Hexagon
% [2] https://en.wikipedia.org/wiki/STL_(file_format)
% [3] https://wps3dprinter.wordpress.com/designing/stl-units/
% [4] https://www.alibreforum.com/forum/index.php?threads/stl-file-units-and-size-problems.19469/

function [CEM2DElvMap , pixelSize , origin , CEM3Dmask , ModulatorMountingPosition] = makeCEF(Plan , SaveSTL)

%If undefined, define the spike orientation
if ~isfield(Plan.Spike , 'SpikeOrientation')
  fprintf('Default spike orientation towards isocentre \n')
  Plan.Spike.SpikeOrientation = 'PATIENT_SIDE';
end

param = getMachineParam(Plan.BDL); %Get the snout parameter for the design of the CEM

%Loop for each beam
for b = 1:numel(Plan.Beams)
    fprintf('Creating Conformal Energy Filter for beam %d \n',b)

    %Find the X,Y voxel size and number of pixels
    % coordinate (in IC Gantry CS) of the centre of the Spikes are already projected onto the plane of the CEF, taking into account the anamorphic projection
    SpikeCentres = [[Plan.Beams(b).RidgeFilter(:).x_centre]' , [Plan.Beams(b).RidgeFilter(:).y_centre]' ];

    NbSpikes = size(SpikeCentres,1);
    minCentre = min(SpikeCentres); %Lower bound for centre of Spikes
    maxCentre = max(SpikeCentres); %Higher bound for centre of Spikes
    maxZCEF = max([Plan.Beams(b).RidgeFilter(:).h_step],[],'all'); %Height of the tallest spike

    if isfield(Plan.Beams(b), 'RangeCompensator')
      %If a range compensator is defined, then add it to the CEM
      RangeCompensator = Plan.Beams(b).RangeCompensator;
      RangeCompensator.IsocentertoCompensatorTrayDistance = Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance; % The range compensator is projected into the plane of the CEM
      RangeCompensator.BDL = Plan.BDL;

      maxZCEF = maxZCEF + max(RangeCompensator.CompensatorThicknessData , [] , 'all'); %Make the grid higher to accomodate the range shifter

    else
      RangeCompensator = [];
    end

    stepZ = Plan.Spike.intrpCTpxlSize; %Height resolution of CEF mask.
    Z_cem = -2:stepZ:maxZCEF+2; %The elevation (mm) of the CEF slices. Go to negative numbers to define a border for STL
      %Allow for space on both side of the CEF so that the isosurface will generate a closed STL object
      %|Z_cem=0| at the base of CEM. |Z_cem| increases when moving up in the CEM

    % Try first guess step
    stepXY = Plan.Spike.intrpCTpxlSize; %Lateral resolution of CEF mask.
      %Define a mask with the same resolution as the hig resolution Ct scan to avoid aliasing problems

    [maxR0 , ~ , nrSides , GridLayout ] = getConvGridParam(Plan , b);
    maxR = maxR0 .* 1.5; %Make the spike grid sufficiently large to avoid clipping
    border = 10; %mm  Add a 10mm border around the spikes in order for the base to be less likely to bend
    wallThick = 1; %mm Thickness of the surrounding wall
    maxBlock = maxCentre + [maxR ,  maxR] +  border .* 2 + wallThick; %mm
    minBlock = minCentre - [maxR ,  maxR] - (border .* 2 + wallThick); %mm

    % adjust number of steps to be an integer
    nrSteps = round(maxBlock ./ (2 .* stepXY)) .*2 + 1; %Make sure this is an odd number so that the matrix copying is symetrical
    stepXY = round(maxBlock ./ nrSteps , 1); %voxel size at 0.1mm resolution

    %Define the dimensions of the Conformal Energy modulator
    Xvec = minBlock(1) : stepXY(1) : maxBlock(1);
    Yvec = minBlock(2) : stepXY(2) : maxBlock(2);

    [Ycem, Xcem] = meshgrid( Yvec , Xvec);
    CEMSize = size(Xcem); %CEMSize(dimX, dimY)

    %Define the dimension of each spike
    Xspike = -maxR:stepXY(1):maxR;
    Yspike = -maxR:stepXY(2):maxR;

    [Ypk,Xpk] = meshgrid(Yspike,Xspike); %meshgrid outputs the data in the format Xpk(y,x). The transpose make sure that the first index of CEMelevation will the the X axis and the second index the Y axis
    pts_Spike = [Xpk(:)';Ypk(:)'];
    SpikeSize = size(Xpk); %SpikeSize(dimX, dimY)


    %The Z height for each X,Y
    CEMelevation = zeros(CEMSize); %CEMelevation(x,y)
    CEF3dMask = zeros(CEMSize(1),CEMSize(2),numel(Z_cem)); %3D mask defining shape of CEF

    %Get the coordinate of the spikes where there is dose in the SOBP
    SptWithDoseX = [Plan.Beams(b).RidgeFilter.x_centre] ; %The position are already projected into the plane of the CEF
    SptWithDoseX = SptWithDoseX([Plan.Beams(b).RidgeFilter.w] > 0); %X Coordinate of spot with dose
    SptWithDoseY = [Plan.Beams(b).RidgeFilter.y_centre] ;
    SptWithDoseY = SptWithDoseY([Plan.Beams(b).RidgeFilter.w] > 0); %Y coordiante of spot with dose
    SptWithDose = [SptWithDoseX' , SptWithDoseY']; %List of the spikes with dose in SOBP


    %Loop for each spike
    %Draw each spike separately and add it to the CEM
    for SpikeIdx = 1:NbSpikes
        fprintf('Spike %d of %d \n',SpikeIdx,NbSpikes)

        [ElvMapSpk , MaskSpk ] = riseOneSpike(Plan , b , SpikeIdx , pts_Spike , SpikeSize , Z_cem  , RangeCompensator);

        if Plan.showGraph
          figure(100)
          surf(Xpk , Ypk , ElvMapSpk)
          title(['Conformal Energy Filter beam ',num2str(b), '-- Spike ', num2str(SpikeIdx) ,' of ', num2str(NbSpikes)])
          xlabel('X (mm)')
          ylabel('Y (mm)')
          zlabel('Height (mm)')

          figure(112)
          image(squeeze(MaskSpk(:,round(size(MaskSpk,2)/2),:)).*255)
          title(['Conformal Energy Filter beam ',num2str(b), '-- Spike ', num2str(SpikeIdx) ,' of ', num2str(NbSpikes)])
          ylabel('X (pxl)')
          xlabel('Height (pxl)')
          drawnow
        end

        %Copy the elevation map of the spike into the CEM map
        wNzero = find(ElvMapSpk ~= 0);
        [IwNzero, JwNzero] = find(ElvMapSpk ~= 0);
        [idx , Xidx , Yidx] = getPxlIndex(SpikeCentres(SpikeIdx,:) , pts_Spike(:,wNzero) , minBlock , stepXY , CEMSize);
        CEMelevation(idx) = ElvMapSpk(wNzero); %Copy only the non zero pixels of the polygon

        %Copy the 3D spike into the 3D mask
        %Apply spike rotation to follow beamlet angle
        Spacing = [stepXY(1) , stepXY(2) , stepZ];
        ImageOrigin = [min(Xspike) , min(Yspike) , min(Z_cem)]; %Origin of the spike coordinate system
        P = getDICOMcoord(MaskSpk, Spacing./2, ImageOrigin , [0,0,maxZCEF]);
        %Get the CEM coordiantes (mm) of all pixel == 1. The rotation is done at the base of the CEF
                                               %MaskSpk has a resolution double of the CEF3dMask

       % align the Z coordinate to the IEC Z gantry
       % so that the divergence of the tower matches the angle of the proton beam
       %-----------------------------------------------------
       %Make sure the Z_cem points in the same direction as Zg
       switch Plan.Spike.SpikeOrientation
         case 'SOURCE_SIDE'
           %The Z_cem is pointing in the same direction as the Zg
         case  'PATIENT_SIDE'
           %The Z_cem is pointing opposite to the Zg
           P(:,3) = -P(:,3); %flip Z_cem. Keep the origin at the base of the hedgehog
       end

        %rotate the spike to mae it aligned with proton beam
        switch Plan.Spike.SpikeType
          case 'fractal'
              %Do not rotate the spikes for the fractal hedgehog

          otherwise
            %In all other cases, align the spike with the beam direction
            M = getRotMatrix(SpikeCentres(SpikeIdx,1) , SpikeCentres(SpikeIdx,2) , Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance , Plan.BDL); %Rotation matrix to move spike in the direction of the scanned proton beam
            P = M * P;
        end

        %Bring back Z_cem in its original direction
        switch Plan.Spike.SpikeOrientation
          case 'SOURCE_SIDE'
            %The Z_cem is pointing in the same direction as the Zg
          case  'PATIENT_SIDE'
            %The Z_cem is pointing opposite to the Zg
            P(:,3) = - P(:,3); %flip Z_cem. Keep the origin at the base of the hedgehog
        end


        P(1,:) = P(1,:) + SpikeCentres(SpikeIdx,1); %Move spike to its position in IEC gantry
        P(2,:) = P(2,:) + SpikeCentres(SpikeIdx,2);
        P(3,:) = P(3,:) + maxZCEF; %Move back to the base of the CEF

        ImageOrigin = [min(Xvec) , min(Yvec) , min(Z_cem)]; %Origin of the mask for full CEF
        Axyz = DICOM2PXLindex(P' , Spacing , ImageOrigin , true); %Convert back into pixel coordinates for the full CEF 3D mask
        idx = sub2ind(size(CEF3dMask) , Axyz(:,1) , Axyz(:,2) , Axyz(:,3)); %Convert [X,Y,Z] coordinates into a linear index
        CEF3dMask(idx) = 1; %Set the spike into the 3D mask

    end %for SpikeIdx


    if Plan.showGraph
      figure(112)
      image(squeeze(CEF3dMask(:,round(size(CEF3dMask,2)./2),:)).*255)
      title('Slice through the CEF')
    end

    % Interpolate and extrapolate the height of the pixels that are still =0
    % Their height will be equal to the height of the closest non-zero pixel
    fprintf('Dilating the block .... \n')
    NbIter = 2.* round(Plan.Beams(b).sigmaAtCEF.sigma ./ min(3 .* stepXY)); %Dilate the CEF by spot sigma at CEF base beyond the spike contours
    SE = strel('disk',3,0); %Convolution kernel is a disk with 3 pixel radius
    dilCEMelevation = ~~CEMelevation; %Convert CEMelevation into a binary mask

    %The dilation is applying iteratively a convolution filter of 3 pixels
    %in order to reduce the computation time. It is faster than a single convolution with a big kernel
    for it = 1:NbIter
        fprintf('iteration %d of %d \n',it,NbIter)
        dilCEMelevation = imdilate(dilCEMelevation,SE);
    end

    idxExtrPol = find(~CEMelevation .* dilCEMelevation); %linear index of the voxels to be interpolated or extrapolated: dilate CEMelevation~=0 by 2 sigma
    [~,idxClosestPxl] = bwdist(CEMelevation); %|idxClosestPxl(x,y)| is the linear index (in CEMelevation) of the non zero pixel which is closest to (x,y)
    CEMelevation(idxExtrPol) = CEMelevation(idxClosestPxl(ind2sub(size(CEMelevation),idxExtrPol))); %Replace the pixel to be interpolated by the value of the closest non zero pixel
    CEF3dMask = insertObjectIn3dMask(idxExtrPol , CEMelevation , CEF3dMask , Z_cem);
    fprintf('done \n')


    % Set base height to min thickness
    wzero = find(CEMelevation == 0);
    CEMelevation(wzero) = Plan.Spike.MinThickness;
    idx = find((Z_cem <= Plan.Spike.MinThickness) .* (Z_cem >= 0));
    CEF3dMask(:,:,idx) = 1; %Fill the mask of the spike in the 'CT' of the CEF

    %Remove the resin below the base of the CEF
    idx = find(Z_cem < 0);
    CEF3dMask(:,:,idx) = 0;

    %Remove the pixels which are beyond the maximum radius of the snout
    mask = (Ycem.^2 + Xcem.^2) <= param.snout.CEMmaxRadius.^2; %Identify the voxels inside the ring of the holder
    CEMelevation = CEMelevation .* mask; %Clip the mask which is located outside of the allowed circle
    mask3D = repmat(mask,1,1,size(CEF3dMask,3));
    CEF3dMask = CEF3dMask .* mask3D; %REmove the voxels out of the allowed circle in the 3D mask

    %Remove voxels all around the base in order to close the STL object
    CEF3dMask(1,:,:) = 0; %Remove the resin below the base of the CEF
    CEF3dMask(end,:,:) = 0; %Remove the resin below the base of the CEF
    CEF3dMask(:,1,:) = 0; %Remove the resin below the base of the CEF
    CEF3dMask(:,end,:) = 0; %Remove the resin below the base of the CEF

    %Draw 2 orthogonal sections going through the isocentre
    if Plan.showGraph
        figure(59+b)
        DistToIso = sum(SpikeCentres.^2,2);
        [~ , NearIso] = min(DistToIso); %Find the index of the spot closest to isocentre
        centralSpike = SpikeCentres(NearIso,:); %Coordinate of the spike closest to optic axis
        [~ , Xc] = min((Xvec - centralSpike(1)).^2);
        [~ , Yc] = min((Yvec - centralSpike(2)).^2);

        plot(Xvec,CEMelevation(:,Yc))
        xlabel('X (mm)')
        ylabel('Height (mm)')
        title(['Thickness profile at Y=',num2str(Yvec(Yc)),'mm'])
        grid on

        figure(60+b)
        plot(Yvec,CEMelevation(Xc,:))
        xlabel('Y (mm)')
        ylabel('Height (mm)')
        title(['Thickness profile at X=',num2str(Xvec(Xc)),'mm'])
        grid on

        figure(112)
        image(squeeze(CEF3dMask(:,round(size(CEF3dMask,2)./2),:)).*255)
        title('Slice through the CEF')
      end

    %Export the STL file to disk
    if (SaveSTL)
      pathName = getOutputDir(Plan.output_path , b);
      filename = fullfile(pathName,[matlab.lang.makeValidName(Plan.Beams(b).RangeModulator.AccessoryCode),'.stl'])
      exportCEM2STL(CEF3dMask , round([stepXY , stepZ] , 2) , ...
                                [rounding(minBlock(1),stepXY(1)) , rounding(minBlock(2),stepXY(2)) , 0] ,...
                                Plan.Beams(b).RangeModulator.AccessoryCode , filename)

end

    %Draw a contour plot of the CEF
    if Plan.showGraph

        CEMcontourPlot(19+b , Xvec, Yvec, CEMelevation , Plan.Beams(b).BlockData , Plan.Beams(b).VDSA , Plan.Beams(b).RangeModulator.IsocenterToRangeModulatorDistance)
        plot(SptWithDoseX , SptWithDoseY,'+k'); % '+' where dose is delivered
        plot(SpikeCentres(:,1) , SpikeCentres(:,2) , 'ok'); % 'o' where a spot was planned
      end

    % Define the orientation of the CEM spikes
    %The modulator (CEM) has a flat base and spike on the other direction
    %This is the flat base that should be used to define the distance from isocenter.
    switch Plan.Spike.SpikeOrientation
      case 'SOURCE_SIDE'
        %The CEM points towards the source. The flat base is fixed on the 'SOURCE_SIDE' of the modulator tray
        ModulatorMountingPosition = 'SOURCE_SIDE'; %The modulator is pointing towards the patient.
      case  'PATIENT_SIDE'
        %the CEM poitns towards the patient. The flat base is fixed on the 'PATIENT_SIDE' of the modulator tray
        ModulatorMountingPosition = 'PATIENT_SIDE'; %The modulator is pointing towards the patient.
      otherwise
        Plan.Spike.SpikeOrientation
        error('Unknown Plan.Spike.SpikeOrientation')

    end

    %Check that the total height of the CEF will fit in the FLASH snout
    if max(CEMelevation,[],'all') > param.snout.CEMmaxHeight
      fprintf('Maximum height of the CEM            : %3.2f mm \n',max(CEMelevation,[],'all'))
      fprintf('MAximum clearance in the FLASH snout : %3.2f mm \n',param.snout.CEMmaxHeight)
      error('No enough space to fit the CEM into the FLASH snout')
    end

    CEM2DElvMap{b} = CEMelevation;
    CEM3Dmask{b} = CEF3dMask; %The index z increases in the direction of increasing Zg. x and y indices are aligned with the IEC gantry CS
    pixelSize{b} = round([stepXY , stepZ] , 2); %Round to 2 significative digit below mm

    origin{b} = [rounding(minBlock(1),pixelSize{b}(1)) ,...
                  rounding(minBlock(2),pixelSize{b}(2)) , ...
                  rounding(min(Z_cem) , pixelSize{b}(3))];  %Coordinate of the first pixel of the CEF

  end %for b
end %end of function


%========================================================
function [idx , Xidx , Yidx] = getPxlIndex(SpikeCentres , pts_Spike , minBlock , stepXY, CEMSize)
 NbPxl = size(pts_Spike,2);

 PxlPos = pts_Spike + repmat(SpikeCentres',1,NbPxl); %Position [x,y] (mm) in the block of each pixel of the Spike grid map
 PxlIdx = (PxlPos - repmat(minBlock',1,NbPxl));

 Xidx = 1 + round(PxlIdx(1,:)./ stepXY(1));
 Yidx = 1 + round(PxlIdx(2,:)./ stepXY(2));

 idx = sub2ind(CEMSize  , Xidx , Yidx);

end


%===========================
% Insert an object defined by a 2D elevation map into a 3D mask
%
% INPUT
% |idxPxl2d| -_SCALAR VECTOR_- |idxPxl2d(i)| linear index of the i-th pixel in 2D map |CEMelevation| that must be copied into the 3D map
% |CEMelevation|  -_SCALAR VECTOR_- |CEMelevation(i)| Elevation of the i-th pixel. The voxels of |CEF3dMask| at location xi,yi and from Z_cem=0 to Z_cem = |CEMelevation(i)| will be set =1
% |CEF3dMask| -_SCALAR MATRIX_- |CEF3dMask(x,y,z)| The voxel at (x,y,z) =1 if they are included inside the CEF
% |Z_cem| -_SCALAR VECTOR_- |Z_cem(k)| is the elevation of the voxels in slice |CEF3dMask(:,:,k)|
%===========================
function CEF3dMask = insertObjectIn3dMask(idxPxl2d  , CEMelevation , CEF3dMask , Z_cem)

  for idx = 1:numel(idxPxl2d)
   [xi,yi] = ind2sub(size(CEMelevation),idxPxl2d(idx));
   idxZ = find( (Z_cem <= CEMelevation(idxPxl2d(idx)) ) .* (Z_cem >= 0) );
   CEF3dMask(xi,yi,idxZ) = 1; %Fill the mask of the spike in the 'CT' of the CEF
  end

end

%==============================
% Add a label in the 2D mask and the 3D mask
%
%INPUT
% |stepXY| -_SCALAR_- pixel size (mm)
% |border| -_SCALAR_- size (mm) of the border added around the spikes
% |LabelText| -_STRING_- Text to be engraved i nthe CEF
% |BaseThickness| -_SCALAR_- Thickness of the base on which the label is added
% |CEMelevation|  -_SCALAR VECTOR_- |CEMelevation(i)| Elevation of the i-th pixel. The voxels of |CEF3dMask| at location xi,yi and from Z_cem=0 to Z_cem = |CEMelevation(i)| will be set =1
% |CEF3dMask| -_SCALAR MATRIX_- |CEF3dMask(x,y,z)| The voxel at (x,y,z) =1 if they are included inside the CEF
% |Z_cem| -_SCALAR VECTOR_- |Z_cem(k)| is the elevation of the voxels in slice |CEF3dMask(:,:,k)|
%==============================
function [CEMelevation , CEF3dMask] = addLabelToCEM(stepXY ,border , LabelText , direc, BaseThickness , CEMelevation , CEF3dMask , Z_cem)

  NbPxlLabel = round(border ./ (2 .* min(stepXY))); %The width of the label is half the number of pixels of the border
  FontSize = round(NbPxlLabel ./2); %Font size in nb pixels
  %Xoffset = 10; %pxl Offset of the origin of the label wrt to CEF edge
  Label = BaseThickness .* ones(NbPxlLabel , round(FontSize .* numel(LabelText))); %Prepare a blank label
  LabelThickness = BaseThickness .* 2; %The height of the "engraving" will be twice the thickness of the base

  %Label = AddTextToImage(Label , String       ,Position,  Color          , Font    , FontSize)
  Label  = AddTextToImage(Label , LabelText , [1,1] , LabelThickness , 'Arial' , FontSize);
  sizeCEMelevation = size(CEMelevation);

  %Remove any trailing spaces
  Ymax = max(find(sum(Label-BaseThickness,1))); %This is the last pixel used in Y
  Label = Label(:,1:Ymax); %Clip the label to remove all trailing space. This will make centering easier

  switch direc
    case 'X'
      %Rotate the label by 90° to place it along the X axis
      Label = Label';
    case 'Y'
      %Label already correctly aligned
    end

  sizeLabel = size(Label);
  if sizeLabel(1) > sizeCEMelevation(1) % if label size in X direction is larger than the available space in X, then cut of the label to fit
    Label = Label(1:sizeCEMelevation(1),:);
  end

  if sizeLabel(2) > sizeCEMelevation(2) % if label size in Y direction is larger than the available space in Y, then cut of the label to fit
    Label = Label(:,1:sizeCEMelevation(2));
  end

  Xpos = 1:size(Label,1);
  Ypos = 1:size(Label,2);

  %CEnter the label onto its axis
  switch direc
    case 'X'
      XposCent = Xpos + round(sizeCEMelevation(1) ./2) - round(size(Label,1) ./ 2); %Shift to middle of the edge
      YposCent = Ypos;
    case 'Y'
      XposCent = Xpos;
      YposCent = Ypos + round(sizeCEMelevation(2) ./2) - round(size(Label,2) ./ 2);
    end


  CEMelevation(XposCent,YposCent) = Label; %Copy the label on the height map of the CEF along the Y axis
  %The label will also be "imprinted" in the range modulator map exported in the DICOM RT-ion plan
  [yi,xi] = meshgrid(YposCent,XposCent);
  CEF3dMask = insertObjectIn3dMask(sub2ind(size(CEMelevation),xi(find(Label==LabelThickness)),yi(find(Label==LabelThickness)))  , CEMelevation , CEF3dMask , Z_cem);


end

%====================================
% Get the rotation matrix to convert
% from the spike coordinate system (daughter)
% to the IEC gantry CS (mother)
% so that the main axis of the spike is aligned wiht the proton beam axis
%
% INPUT
% |X| -_SCALAR_- X coordinate (IEC gantry, mm) of the spike in the plane of isocentre
% |Y| -_SCALAR_- Y coordinate (IEC gantry, mm) of the spike in the plane of isocentre
% |IsocenterToRangeModulatorDistance| -_SCALAR_- Distance from isocentre to the base of the CEF
% |BDL_file| -_STRING_- path to the BDL file
%====================================
function M = getRotMatrix(X , Y , IsocenterToRangeModulatorDistance , BDL_file)

  az = IsocenterToRangeModulatorDistance;
  mag = magnification(az , BDL_file); %Magnification facotr from isocentre plane to CEF

  ax = X .* (1-mag(1));
  alpha = - rad2deg(atan(ax ./ az)); %Rotation around the Y IEC gantry axis

  bz = sqrt(az.^2 + ax.^2);
  by = Y .* (1-mag(2));
  bet = rad2deg(atan(by ./ bz)); %Rotation around the X IEC gantry axis

  M = roll(alpha,[0,0,0]) * pitch(bet,[0,0,0]); %Rotation matrix IECgantry = M * spike
        % M = Ry * Rx

end

%===============================================
% Add 4 fiducial markers in the elevation map and in the 3D mask
%
% INPUT
% |CEMelevation|  -_SCALAR MATRIX_- |CEM2DElvMap{b}(x,y)| Physical height (mm) of the CEF at the position (x,y) (Gantry CS)  in the plane of the CEF for beam b. Voxel ordered in increasing coordinates
% |CEF3dMask| -_SCALAR MATRIX_- |CEM3Dmask{b}(x,y,z)| 3D mask of the CEF.  |CEM3Dmask{b}(x,y,z)=1| if the voxel at location (x,y,z)  in the plane of the CEF for beam b belongs to the CEF.
% |FiduHeight| -_SCALAR_- Height (pixel) of the fiducial markers
% |border| -_SCALAR_- Width (pixel) of the border around the CEF
% |Z_cem| -_SCALAR VECTOR_- |Z_cem(k)| is the elevation of the voxels in slice |CEF3dMask(:,:,k)|
%
% OUTPUT
% |CEMelevation|  -_SCALAR MATRIX_- |CEM2DElvMap{b}(x,y)| Physical height (mm) of the CEF at the position (x,y) (Gantry CS)  in the plane of the CEF for beam b. Voxel ordered in increasing coordinates
% |CEF3dMask| -_SCALAR MATRIX_- |CEM3Dmask{b}(x,y,z)| 3D mask of the CEF.  |CEM3Dmask{b}(x,y,z)=1| if the voxel at location (x,y,z)  in the plane of the CEF for beam b belongs to the CEF.
%
%===============================================
function [CEMelevation , CEF3dMask] = addFiducials(CEMelevation , CEF3dMask , FiduHeight , border , Z_cem)

    sCEF = size(CEMelevation);
    loc = round(border ./2);
    mask = ones(5,5) .* FiduHeight; %Mask defining the shape of the fiducal marker
    pxl = -2:2;

    mark(1,:) = [loc         , loc];
    mark(2,:) = [sCEF(1)-loc , loc];
    mark(3,:) = [loc         , sCEF(2)-loc];
    mark(4,:) = [sCEF(1)-loc , sCEF(2)-loc];

    for m = 1:size(mark,1)
      CEMelevation(mark(m,1)+pxl , mark(m,2)+pxl) = mask;
      [Y,X] = meshgrid(mark(m,2)+pxl , mark(m,1)+pxl);
      idx = sub2ind(size(CEMelevation) , X , Y);
      CEF3dMask = insertObjectIn3dMask(idx , CEMelevation , CEF3dMask , Z_cem);
    end

end


%=================================
% Draw a rise rectangle on the CEF
%=================================
function maskRect = drawRectangle(X , Y , border)
  maxX = max(X,[],'all');
  maxY = max(Y,[],'all');
  minX = min(X,[],'all');
  minY = min(Y,[],'all');
  maskRect = (X > (minX + border)) .* (X < (maxX - border)) .* (Y > (minY+ border)) .* (Y < (maxY - border));
end


%=================================
% Draw a rise face on the CEF for support
%=================================
function maskFaceWall = drawFaceWall(X, border)

      minX = min(X,[],'all');
      maxX = max(X,[],'all');
      maskFaceWall = (X > (minX + border));
end


%=================================
% Draw a rise face on the CEF for support
%=================================
function [CEMelevation , CEF3dMask] = drawFeetSupport(CEMelevation , CEF3dMask , FiduHeight, Z_cem)

    sCEF = size(CEMelevation);
    loc = 16;
    mask = ones(31,31) .* FiduHeight; %Mask defining the shape of the fiducal marker
    pxl = -15:15;

    mark(1,:) = [loc, loc];
    mark(2,:) = [loc, sCEF(2)-loc];

    for m = 1:size(mark,1)
      CEMelevation(mark(m,1)+pxl , mark(m,2)+pxl) = mask;
      [Y,X] = meshgrid(mark(m,2)+pxl , mark(m,1)+pxl);
      idx = sub2ind(size(CEMelevation) , X , Y);
      CEF3dMask = insertObjectIn3dMask(idx , CEMelevation , CEF3dMask , Z_cem);
    end
end


%=================================
% Draw a positioning markers on the CEF base
%=================================
function [CEMelevation , CEF3dMask] = addMarkers(CEMelevation , CEF3dMask , FiduHeight, Z_cem)

    sCEF = size(CEMelevation);
    loc = 7;
    mask = ones(13,13) .* FiduHeight; %Mask defining the shape of the fiducal marker
    pxl = -6:6;

    %mark(1,:) = [loc, loc];
    mark(1,:) = [round((sCEF(1)-loc)/2), loc];
    mark(2,:) = [sCEF(1)-loc, loc];
    mark(3,:) = [sCEF(1)-loc, round((sCEF(2)-loc)/2)];
    mark(4,:) = [sCEF(1)-loc, sCEF(2)-loc];
    mark(5,:) = [round((sCEF(1)-loc)/2), sCEF(2)-loc];
    %mark(7,:) = [loc, sCEF(2)-loc];
    mark(6,:) = [loc, round((sCEF(2)-loc)/2)];

    for m = 1:size(mark,1)
        CEMelevation(mark(m,1)+pxl , mark(m,2)+pxl) = mask;
        [Y,X] = meshgrid(mark(m,2)+pxl , mark(m,1)+pxl);
        idx = sub2ind(size(CEMelevation) , X , Y);
        CEF3dMask = insertObjectIn3dMask(idx , CEMelevation , CEF3dMask , Z_cem);
    end

end
