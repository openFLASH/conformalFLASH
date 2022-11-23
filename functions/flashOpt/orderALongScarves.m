%% orderALongScarves
% Define the scanning order of spots by placing paralell scarves
% * Identify the major axis of the scanned area
% * Set the main axis of the scarves paralell to the major axis
% * Place multiple scarves next to each other, aligned with the major axis
% * The number of paralell scarves is qual to the minor axis length divided by the width of a scarf
%
%% Syntax
% |OrderedSpot  = orderALongScarves(spot, ScarfWidth , plotNb)|
%
%
%% Description
% |OrderedSpot  = orderALongScarves(spot, ScarfWidth , plotNb)| Description
%
%
%% Input arguments
% |spot| - _SCLAR MATRIX_ - The i-th spot to deliver is spot(i,:) = [x,y]
%
% |ScarfWidth| -_SCALAR_- Maximum width (in |spot| length units) of a scarf to obtain high dose rate inside a scarf
%
% |NbScarves| -_NTEGER_- Number of scarve to lay paralell to the main ellipse axis. Each scarf will be scanned sequentially
%
% |plotNb| -_SCALAR_- [OTPIONAL] Figure number where to plot the trajectory
%
%
%% Output arguments
%
% |OrderedSpot| - _SCLAR MATRIX_ - The ordered sequence. i-th spot to deliver is OrderedSpot(i,:) = [x,y]
%
% |OrderedCoord| -_INTEGER VECTOR_- Ordered sequence of spot indices |OrderedSpot = spot(OrderedCoord,:)|
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [OrderedSpot , OrderedCoord] = orderALongScarves(spot, V , NbScarves ,  plotNb)

  if (nargin < 3)
    plotNb = [];
  end

  [Lat , Ts ] = getLatticeInfo(spot(:,1:2)); %Get lattice info in the space domain

  CoordType = 'space'; %The trajectory is optimised in physical distance
  %CoordType = 'time'; %The trajectory is optimised in time distance
  switch(CoordType)
    case 'space'
      spotT = spot;
    case 'time'
      Vap = getV(Lat , V);
      spotT = spot ./ repmat(Vap,size(spot,1),1); %Convert distances into approximate delivery times
  end

  stats = fitEllipse(spotT(:,1:2)); %Find the properties of the ellipse enclosing all PBS spots

  [latticeSpacing , Ts , GridSize] = getLatticeInfo(spotT(:,1:2)); %Get lattice info in the time domain

  if isempty(NbScarves)
    %No number of scarve is defined. Let's compute it
    ScarfWidth = 100; %ms
    NbScarves = round(stats.MinorAxisLength ./ ScarfWidth); %Number of paralell scarves. Round nb of scarves to closest integer
    if (NbScarves==0)
      NbScarves =1;
    end
  else
    ScarfWidth = 1.3 .* stats.MinorAxisLength ./ NbScarves; %Cover a width larger than the minior axis length with the scarves
    if (ScarfWidth==0)
      ScarfWidth =100; %ms
    end
  end
  fprintf('Scarf width = %d ms \n',ScarfWidth)

  %Loop for each scarf and identify the spine axis of the scarf
  alph = deg2rad(stats.Orientation); %Angle between ellipse main axis axis and X axis
  Wvec = [cos(alph+90) , sin(alph+90)];
  Wvec = Wvec ./ norm(Wvec); %Unit Vector pointing along the width of the scarf

  Lvec = [cos(alph) , sin(alph)];
  Lvec = Lvec ./ norm(Lvec); %Unit Vector pointing along the spine axis of the scarf

  %Create the branches = the main axis of the scarf
  direct = 1;
  for scf = 1:NbScarves
    %Loop for each paralell scarf
    FirstPt = stats.Centroid + Wvec .* ScarfWidth .* ((scf-1) + 0.5 - NbScarves./2);  %point on the scarf axis and aligned with centroid
    spts = walkScarfAxis(FirstPt , Lvec , stats.MajorAxisLength , GridSize , latticeSpacing , spotT(:,1:2) , Ts , direct);
    branches.Path{scf} = spts;
    branches.BranchLength(scf) = size(branches.Path{scf},1);
    direct = direct .* (-1); %Walk sucessive scarves so that they are head to tail. This will increase dose rate in the part of the overlap region between scarves
  end

  %There are no branching points
  branches.branchingPoint = zeros(numel(branches.BranchLength),1);
  branches.BranchingCoord = [];
  BranchScanOrder = 1:NbScarves;

  %Plot the red line representing the main axis of the scraf
  if (~isempty(plotNb))
    figure(plotNb)
    for i = 1:numel(branches.Path)
      hold on
      pt = branches.Path{i};
      plot(pt(:,1),pt(:,2),'-^r')
    end
  end

  %Build the ribs perpendicular to the branches (= scarf axes)
  OrderedCoord = walkBranches(spotT , branches , ScarfWidth./2, BranchScanOrder);
  OrderedSpot = spot(OrderedCoord,:);

  % hold on
  % plot(spotT(OrderedCoord,1),spotT(OrderedCoord,2),'-k')
  % pause

end

%==========================================
% Walk by one unit at a time along the identified main axis of the scarf
% At each step, find the closest point in the mapImage
% Check that the pixel indices are >0 and smaller than mapImage size
%
% INPUT
% FirstPt
% |Lvec|  -_SCALAR VECTOR_- UInit vector aligned with the major axis of the ellipse
% |MajorAxisLength| -_SCALAR_- 	Length (in mm) of the major axis of the ellipse that has the same normalized second central moments as the region
% |GridSize| -_SCALAR MATRIX_- Dimension of the grid on which the spot are placed |GridSize(:,1) = [minX , minY]| |GridSize(:,2) = [maxX , maxY]|
% |latticeSpacing| _SCALAR VECTOR_ - Spatial period (mm) of the spot lattice. |lat(i)| is the spacing along vectior Ts(i,:)
% |spotCoord| - _SCLAR MATRIX_ - The i-th spot to deliver is spot(i,:) = [x,y]
% |Ts| -_SCALAR MATRIX_- |Ts(i,:)= [x,y]| Unit vector defining the i-th latice vector.
%                           |Ts(1,:)| is mostly paralell to the X axis
%                           |Ts(2,:)| is mostly paralell to the Y axis
% |direct| -_SCALAR_- +1 or -1 Define the direction of Lvec. In order to scan the spine in different directions
%
% OUTPUT
% |seq| - _SCLAR MATRIX_ - The i-th spot to deliver is seq(i,:) = [x,y]
%==========================================
function seq = walkScarfAxis(FirstPt , Lvec , MajorAxisLength , GridSize , latticeSpacing , spotCoord , Ts , direct)

  seq = [];

  idxL = getTalignedWithL(Ts , Lvec); %Latice vector that is mostly aligned with scarf axis
  Ta = Ts(idxL,:); %Lattice vector mostly paralell to Lvec
  idx = [1,2];
  idx(idxL) = [];
  Tb = Ts(idx,:); %Lattice vector mostly perpendicular to Lvec

  % By moving along the Lvec axis (Lvec is mostly paralell to Ta), we will intersect the periodic lattice lines paralell to Tb
  % The intersection point between the Lvec vector and the Tb period line is given by the vector equation:
  % w .* Lvec = Ta + b .* Tb
  % where w and b are scalars. We can rewrite the equation as:
  % w .* Lvec - b .* Tb = Ta
  %Let's solve this equation to find w which is the periodic distance along Lvec where the intersection occurs
  X = linsolve([Lvec' , Tb'] , Ta');
  w = abs(X(1));
  if (w==0)
    w=1;
  end


  [d , SpotIndex] = getDistances2(spotCoord , FirstPt); %Find the closest point to the scarf axis
  Central = spotCoord(SpotIndex,1:2); %This is a point on the grid that is the closest to the scarf axis
  stepS = latticeSpacing(idxL) .* w; %Move along spine so as to jump from one lattice column to the next

  fac = 2.5;
  NbSteps = round(fac .* MajorAxisLength ./ (2 .* stepS) ); %Make the offset a multiple of the step
  StartSpine = Central - direct .* Lvec .* NbSteps .* stepS; %Move to the begining of the spine
  L = 0;

  while (L < round(fac .* MajorAxisLength));
    Pt = StartSpine + direct .* Lvec .* L; %The spine does not need to be on the grid. Only the rib must scan on the grid
    seq = [seq ; Pt]; %Add the spot to the end of the sequence
    L = L + stepS; %Jump by one lattice step
  end

% plot(seq(:,1),seq(:,2),'-*g')
% pause

end




%===============================
% Compute apparent scanning speed
% taking intoa ccount the sweep time and the spot delivery time
%=================================
function Vap = getV(Lat , V)
  Td = 20; %ms Assumed average delivery time of a spot %TODO can we have a better estimate ?
  Vap = Lat .* V ./ (Lat + V .* Td);
  fprintf('Apparent scanning speed : X = %f m/s  Y = %f m/s \n',Vap(1),Vap(2));
end
