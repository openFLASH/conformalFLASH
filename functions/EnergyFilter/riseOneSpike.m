%% riseOneSpike
% Compute the 2d elevation map and the 3D mask for one single spike of the CEM.
% different spike geomeries are available:
%     * Aztec pyramid with the centre of the spike corresponds to the BP with smaller range (|Plan.Spike.SpikeType ='up'|)
%     * Aztec pyramid with the centre of the spike corresponds to the BP with the largest range (|Plan.Spike.SpikeType ='down'|)
%     * Cenote geometry |Plan.Spike.SpikeType ='ellipse'|
%     * Fractal (Manathan) geoemtry |Plan.Spike.SpikeType ='fractal'|
%
%% Syntax
% |[CEM2DElvMap , CEM3Dmask , RcThickness , ZgHigh] = riseOneSpike(Plan , b , SpikeIdx , pts_Spike , SpikeSize , Z_cem , RangeCompensator)|
%
%
%% Description
% |[CEM2DElvMap , CEM3Dmask , RcThickness , ZgHigh] = riseOneSpike(Plan , b , SpikeIdx , pts_Spike , SpikeSize , Z_cem , RangeCompensator)| Description
%
%
%% Input arguments
% |Plan| - _struct_ - MIROpt structure with updated information:
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
% |b| -_SCALAR_- Number of the beam to process
%
% |SpikeIdx| -_SCALAR_- Number of the spike |RidgeFilter(k)| to build
%
% |pts_Spike|  -_SCALAR MATRIX_- |pts_Spike(:,i)= [x,y]| Coordinate (mm) of the i-th pixel in |CEM2DElvMap(sub2ind(SpikeSize,i))|
%
% |Z_cem|  -_SCALAR VECTOR_- Z coordinate (mm) in CEM object CS. |Z_cem=0| at the base of CEM. |Z_cem| increases when moving up in the CEM
%
% |RangeCompensator| -_STRUCTURE_- [OPTIONAL. Only needed if a range compensator is to be added to the CEM]
%        * |RangeCompensator.IsocentertoCompensatorTrayDistance| -_SCALAR_- Distance (mm) from isocentre to upstream side of the range compensaot
%        * |RangeCompensator.CompensatorColumns| -_SCALAR_- Number of pixels along the X dimension
%        * |RangeCompensator.CompensatorRows| -_SCALAR_- Number of pixels along the Y dimension
%        * |RangeCompensator.CompensatorPixelSpacing| -_SCALAR VECTOR_- [x,y] size (mm) of the pixels projected into the isocentre plane
%        * |RangeCompensator.CompensatorPosition| -_SCALAR VECTOR_- [x,y] Coordinate (mm) of the voxels in the IEC gantry CS
%        * |RangeCompensator.CompensatorThicknessData| -_SCALAR MATRIX_- |CompensatorThicknessData(x,y)| Height (mm) of the range compensator at pixel (x,y)

%
%% Output arguments
%
% |CEM2DElvMap| - _SCALAR MATRIX_ - |CEM2DElvMap(i,j)| Height (mm) of the column at position (i,j). The height is measured from |Z_cem=0|
%
% |CEM3Dmask|  -_SCALAR MATRIX_- Binary mask |CEM3Dmask(i,j,k) = 1| if the voxel is contained inside the CEM defined by the elevation map
%               |CEM3Dmask(i,j,k)| is inside the CEM if |CEM2DElvMap(i,j) > Z_cem(k)|
%               |CEM3Dmask| has the double resolution of the 2D grid |CEM2DElvMap| and of the |Z_cem| vector
%
% |RcThickness| - _SCALAR MATRIX_ - |RcThickness(i,j)| Height (mm) of the range compensator at position (i,j).
%
% |ZgHigh|  -_SCALAR VECTOR_- Z coordinate (mm) in CEM3Dmask object CS. |ZgHigh=0| at the base of CEM. |ZgHigh| increases when moving up in the CEM
%                         |ZgHigh| has double the resolution of the |Z_cem| vector
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [CEM2DElvMap , CEM3Dmask , RcThickness , ZgHigh] = riseOneSpike(Plan , b , SpikeIdx , pts_Spike , SpikeSize , Z_cem , RangeCompensator)

%Get the spike properties from the plan
[~ , ~ , nrSides , GridLayout ] = getConvGridParam(Plan , b);
SpikeType = Plan.Spike.SpikeType;
RidgeFilter = Plan.Beams(b).RidgeFilter(SpikeIdx);
param.Sx_Sy = Plan.Beams(b).sigmaAtCEF.BDL.SpotSize1x ./ Plan.Beams(b).sigmaAtCEF.BDL.SpotSize1y;
param.ang = Plan.Beams(b).sigmaAtCEF.BDL.SpotTilt; %Rotation angle (radian) of the ellipse; %Rotation angle (radian) of the ellipse
param.radius = max(Plan.Beams(b).RidgeFilter(SpikeIdx).a_max);

  if nargin < 7
     RangeCompensator = [] ;
  end

  if min(diff(Z_cem)) <0
    error('Z_cem must be a continuously increasing vector')
  end

  %Create a Z vector (in CEM coordiante) with a resolution 2 x resolution of Z_cem
  stLow = (max(Z_cem) - min(Z_cem) ) ./ (numel(Z_cem)-1); %Spatial resolution of Z_cem, assuming uniform spacing of the points
  stHigh = stLow ./ 2; %Step size of the high resolution Z vector.
  ZgHigh = min(Z_cem):stHigh:max(Z_cem); %Height coordinate of the high resolution mask. Make the spacing uniform in case it was not uniform in |Z_cem|


  % Compute the diagonal of one pixel
  %The wall of a chimney will be made wider by a couple of pixels to avoid the pixelisation issues
  %that lead to top of chimney being columns of 1 pixels
  dX = min(diff(unique(pts_Spike(1,:)))); %Smallest X step
  dY = min(diff(unique(pts_Spike(2,:)))); %Smallest Y step
  dR = sqrt(dX.^2 + dY.^2) .* 0;

  %Create an empty elevation map
  CEM2DElvMap = zeros(SpikeSize);

  %Rise floor of whole spike to halp modulation depth
  switch GridLayout
    case 'HEXAGONAL'
      nrSidesLattice = 6;
    case 'SQUARE'
      nrSidesLattice = 4;
  end

  if (RidgeFilter.w ~= 0)
    %There is a SOBP delivered at this spike position
    %Do not add spots that have zero weight

    %Draw the zone out side of the chimney / spike. This has a shape defined by the lattice. It has a height defined by |RidgeFilter.h_outside|
    [CEM2DElvMap2, ~] = drawPolygon(nrSidesLattice , pts_Spike, SpikeSize, RidgeFilter.latticePeriod./2 , RidgeFilter.h_outside , param);
    idx = find(CEM2DElvMap2);
    CEM2DElvMap(idx) = CEM2DElvMap2(idx);

    %The following layers receive a hexagon/square
    switch SpikeType %All spike have the same shape. We cann switch based on type of 1st spike

      case 'fractal'

          %Multiple Manathan towers within one cell
          colWidth = 2 .* max(RidgeFilter.a_max) ./ RidgeFilter.NbColumns; %mm width of individual spikes
          L=1;
          for Xi = 0:RidgeFilter.NbColumns-1
            for Yi = 0:RidgeFilter.NbColumns-1
                X = -max(RidgeFilter.a_max) + colWidth./2 + Xi .* colWidth;
                Y = -max(RidgeFilter.a_max) + colWidth./2 + Yi .* colWidth;
                pts_Spike_offst(1,:) = pts_Spike(1,:) - X;
                pts_Spike_offst(2,:) = pts_Spike(2,:) - Y;
                CEM2DElvMap = risePolygon(colWidth ./2 , nrSides, CEM2DElvMap, pts_Spike_offst , SpikeSize , RidgeFilter.h_grid(L) , param);
                L = L+1;
            end
          end

      otherwise
        %The spike have a stepped shape
        %Collect small steps together till their width reaches the map resolution
        [a_min , a_max , h_step ] = clip2SpikeResolution(RidgeFilter.a_max , RidgeFilter.a_min , RidgeFilter.h_step , min([dX,dY]));
        if(prod(diff(a_max)<0))
          %Radii are in increasing order. Flip these vectors. They must be sorted in decreasing order
          a_min = flipdim(a_min,1);
          a_max  = flipdim(a_max,1);
          h_step  = flipdim(h_step,1);
        end

        % Sort the steps from maximum radius to minimum radius
        % We draw the largest polygon and then draw polygon with more and more reduced diameter
        [a_max , Isort] = sort(a_max,'descend');
        h_step = h_step(Isort);

        %The spike is made of a sequence of steps
        %The radii must be in decreasing order: from largest to smallest
        flagFisrtStep = 1;
        for L = 1:numel(a_max)
            %Loop for each hexagon/square of spike |SpikeIdx|
              if flagFisrtStep
                  %This is the most external step. Make it 2 pixels wider to avoid having a too narrow border
                  %This will avoid the top of the chimney to look like a column of pixels.
                  %The most external step is anyway larger than 2 sigma, so there is supposed to be no proton here
                  CEM2DElvMap = risePolygon(a_max(L) + dR , nrSides, CEM2DElvMap, pts_Spike, SpikeSize , h_step(L) , param);
                  flagFisrtStep = 0;
              else
                  CEM2DElvMap = risePolygon(a_max(L) , nrSides, CEM2DElvMap, pts_Spike, SpikeSize , h_step(L) , param);
              end
        end %for L

    end %switch / case

  end %if (RidgeFilter.w ~= 0)

  %Add the shape of the range compensator
  if ~isempty(RangeCompensator)
    % Coordinate in the isocentre plane
    Xrc = (1:RangeCompensator.CompensatorColumns) .*  RangeCompensator.CompensatorPixelSpacing(1) + RangeCompensator.CompensatorPosition(1); %Coordinate (mm) of the voxels in the IEC gantry CS
    Yrc = (1:RangeCompensator.CompensatorRows) .*  RangeCompensator.CompensatorPixelSpacing(2) + RangeCompensator.CompensatorPosition(2);

    %Project into the plane of the range modulator
    mag = magnification(RangeCompensator.IsocentertoCompensatorTrayDistance , RangeCompensator.BDL); %The magnification is to the base of hte CEM
    Xrc = Xrc .* mag(1); %Project range compensator shape to the plane of the CEM
    Yrc = Yrc .* mag(2);

    Xspk = reshape(pts_Spike(1,:),SpikeSize) + RidgeFilter(1).x_centre;
    Yspk = reshape(pts_Spike(2,:),SpikeSize) + RidgeFilter(1).y_centre;

    RcThickness = interpn(Xrc,Yrc,RangeCompensator.CompensatorThicknessData,Xspk,Yspk,'linear',0); %Interpolate the compensator thickness on the grid of the CEM
    RcThickness(CEM2DElvMap ==0) = 0; %Add the range compensator only to pixels where there is a CEM added

    CEM2DElvMap = CEM2DElvMap + RcThickness;
  else
    %If there is no range compensator, then the thickness is zero
    RcThickness = 0;
  end


  %Add the CEF in the mask
  CEM2DElvMapHigh = imresize(CEM2DElvMap,2,'nearest' ); %Double the resolution of the 2D map
  CEM3Dmask = elevationMap2mask3D(CEM2DElvMapHigh , ZgHigh);

end

%===========================
% Draw the polygon inside the provided |CEM2DElvMap|
% Contrarily to |drawPolygon.m| which create a new |CEM2DElvMap| at each call,
% the function |risePolygon| add a polygon to the provided |CEM2DElvMap|
%
%INPUT
% apothem_spike -_SCALAR_- Apothem of the polygon
% nrSides -_SCALAR_- Number of sides of the polygon
% CEM2DElvMap -_SCALAR MATRIX_- |CEM2DElvMap(x,y)| Height of the pixel (x,y) of the spike
% pts_Spike  -_SCALAR MATRIX_- |pts_Spike(:,i)= [x,y]| Coordinate (mm) of the i-th pixel in |CEM2DElvMap(sub2ind(SpikeSize,i))|
% SpikeSize -_SCALAR_- SpikeSize = [dimX, dimY] number of pixel along each dimension of the spike
% h_step -_SCALAR_- height of the centre of the step of the spike
%
% OUTPUT
% CEM2DElvMap -_SCALAR MATRIX_- Updated |CEM2DElvMap(x,y)| Height of the pixel (x,y) of the spike
%===========================
function CEM2DElvMap = risePolygon(apothem_spike , nrSides, CEM2DElvMap, pts_Spike, SpikeSize, h_step , r)

  [~, GridCEM3Dmask] = drawPolygon(nrSides, pts_Spike, SpikeSize, apothem_spike, h_step , r);
  CEM2DElvMap(logical(GridCEM3Dmask)) =  h_step; %Set the height of the ring

end

%===========================================================
%Determine the number of pixel per step of the stair case
% Make sur to get integer number of pixels per step
%
% INPUT
% |weight| -_SCALAR VECTOR_- |weight(i)| weight of PBs spot corresponding to the i-th step
% |CellWidth| -_SCALAR_- Width (mm) of a cell
% |pxlSize| -SCALAR- Size (mm) of one pixel
% |SugNbColumns| -_INTEGER_- Suggested number of columns
%
% OUTPUT
% |NbPxlInSteps| -_INTEGER VECTOR_- |NbPxlInStaircase(i)|  Number of pixelk to allocate to the i-th step
% |NbPxlInStaircase| -_INTEGER_- Total number of pixel in the staircase
% |NbPxlInCell| -_SCALAR_- Number of pixels in a cell
%===========================================================
function [NbPxlInSteps , NbPxlInStaircase , NbPxlInCell] = getNbPixelPerStaircase(weight , CellWidth , pxlSize , SugNbColumns)

  Threshold = 0.1; %Do not care about steps which are < 10% of maximum step width
  BinStepsMore = true; % Flag indicating whether we still need to reduce the number of StepSize
        %We need an integer number of pixels in each step and at the same time keep the proper area ratio between steps
        %In addition, the number of staircases in a cell should be roughly hat was defined in |SugNbColumns|
        %To achieve that, the smallest steps will be binned together into larger steps untill a reasonable number of step is obtained

  while BinStepsMore
      resolution = max(weight) .* Threshold;
      a_max = cumsum(weight);
      a_min = [0 , a_max(1:end-1)];
      [r_min , r_max , ~ , hidx ] = clip2SpikeResolution(a_max , a_min , 1:numel(a_max) , resolution);

      StepWeights = zeros(1,numel(weight));
      StepWeights(hidx) = r_max - r_min; %Step width with the smallest one removed

      NbPxlInSteps = round(StepWeights ./ min( r_max - r_min)); %Nb pixel for each step of staircase so that there is an integer nb of pixel in each step
      NbPxlInStaircase = sum(NbPxlInSteps); %Nb of columns in a staircase to have integer Nb pixel per stair

      NbPxlInCell = floor(CellWidth ./ pxlSize); %Nb of pixels in the cell
      SuggNbPxlInStaircase = NbPxlInCell ./ SugNbColumns; %Suggested Nb pixel in a staircase
      mul = SuggNbPxlInStaircase ./ NbPxlInStaircase;
      if mul >= 0.8
        %The small steps are sufficiently binned
        %We have reached a reasonable number of steps in one staircase
        BinStepsMore = 0;
      else
        %If mul < 1, this indicates that there are too many small steps
        %Try to bin the steps again with a coarser resolution
        Threshold = Threshold + 0.1; %Increase threshold by 10%
        if (Threshold > 1)
          Threshold
          NbPxlInCell
          error('Could not create the staircase ')
        end
      end
  end %while BinStepsMore
  NbPxlInStaircase = round(mul) .* NbPxlInStaircase; %Find a multiple to get close to the suggested Nb of staircase
  NbPxlInSteps = round(mul) .* NbPxlInSteps;

end

%======================================
% Rise one staircase
%======================================
function CEM2DElvMap = riseStaircase(StairDir , h_step , StepSize , coord , CEM2DElvMap , SpikeSize , r)

  for h = 1:numel(h_step)
    %Loop for every step in one staircase
    if prod(StepSize(h,:)) > 0

      %Shift the origin to the centre of the step
      switch StairDir
        case 0
            coord(1,:) = coord(1,:) - StepSize(h,1);
        case 1
            coord(2,:) = coord(2,:) - StepSize(h,2);
        end
        %Draw a rectangle only is a weight was given to that step
        CEM2DElvMap = risePolygon(StepSize(h,:) , 41 , CEM2DElvMap, coord , SpikeSize , h_step(h) , r);
              %                             41 defines a rectangle

        %Shift the origin to the edge of the step
        switch StairDir
          case 0
              coord(1,:) = coord(1,:) - StepSize(h,1);
          case 1
              coord(2,:) = coord(2,:) - StepSize(h,2);
          end
    end
  end

end

%=================================
% GEt the size of the steps of the staircase
%=================================
function [StepSize , pts_Spike_offst ] = getStepSize(StairDir , pts_Spike_offst , dX , dY , NbPxlInStepsX , NbPxlInStepsY , colWidthX , colWidthY)

  switch StairDir
    case 0
        StepSize(:,1) = NbPxlInStepsX .* dX ./ 2; %Size (mm) of each step in X direction
        StepSize(:,2) = colWidthY ./ 2; %Size (mm) of each step in Y direction
        pts_Spike_offst(2,:) = pts_Spike_offst(2,:) - colWidthY ./2;
    case 1
        StepSize(:,2) = NbPxlInStepsY .* dY ./ 2; %Size (mm) of each step in X direction
        StepSize(:,1) = colWidthX ./ 2; %Size (mm) of each step in Y direction
        pts_Spike_offst(1,:) = pts_Spike_offst(1,:) - colWidthX ./2;
    end

  end
