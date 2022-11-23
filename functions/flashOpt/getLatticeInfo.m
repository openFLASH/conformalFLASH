%% getLatticeInfo
% Compute the lattice spacing of PBS spots from the list of all PBS spots
% Compute the 2 vectors defining the lattice grid
%
%% Syntax
% |latticeSpacing = getLatticeInfo(spotCoord)|
%
%
%% Description
% |latticeSpacing = getLatticeInfo(spotCoord)| Description
%
%
%% Input arguments
% |spotCoord| - _SCLAR MATRIX_ - Coordinates (mm) of the i-th spot to deliver is spot(i,:) = [x,y]
%
%
%% Output arguments
%
% |lat| - _SCALAR VECTOR_ - Spatial period (mm) of the spot lattice. |lat(i)| is the spacing along vectior Ts(i,:)
%
% |Ts| -_SCALAR MATRIX_- |Ts(i,:)= [x,y]| Unit vector defining the i-th latice vector.
%                           |Ts(1,:)| is mostly paralell to the X axis
%                           |Ts(2,:)| is mostly paralell to the Y axis
%
% |GridSize| -_SCALAR MATRIX_- Dimension of the grid on which the spot are placed |GridSize(:,1) = [minX , minY]| |GridSize(:,2) = [maxX , maxY]|
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [lat , Ts , GridSize] = getLatticeInfo(spotCoord)

  if (size(spotCoord,1) > 1)

        %There are at least 2 spots
        DistMat = interSpotDistance(spotCoord);

        %Get the main axis of the lattice. The rib will be aligned along these main axes
        [d , idx]= getDistances2(spotCoord(:,1:2) , [0,0]); %Get the point nearest to isocentre
        Iso = spotCoord(idx(1),1:2); %If there are several spots near isocetre, arbitrarily pick one
        dis2Lattice = round(DistMat(idx(1),:),3); %Distance from isocentre to all other spots. Round to the closest um
        dis2Lattice(idx(1))= NaN;  %Remove the distance from isocentre to isocentre
        d2iso =  min(dis2Lattice); %Find the nearest points surrounding the isocentre
        [Tstmp , latticeSpacing] = getAllTs(d2iso , dis2Lattice , spotCoord, Iso);

        if (size(Tstmp , 1) <= 2)
          %The lattice is asymetrical. The lattice period in X is different from Y
          %Find a second Ts
          d = dis2Lattice;
          d(mod(d,d2iso)==0) = []; %Remove multiple of the first period
          d2iso =  min(d); %Find the nearest points surrounding the isocentre
          [TsTstmp2 , latticeSpacing2] = getAllTs(d2iso , dis2Lattice , spotCoord, Iso);
          latticeSpacing = [latticeSpacing ; latticeSpacing2];
          Tstmp = [Tstmp ; TsTstmp2];
        end

        %Get the Ts vectors that are mostly pointing along X and Y axes
        [idX , idY]= getTalignedWithL(Tstmp , [1,0]); %Find the T that is mostly aligned with X axis

        % idX = getTalignedWithL(Tstmp , [1,0]); %Find the T that is mostly aligned with X axis
        if ~isempty (idX)
          Ts1 = Tstmp(idX,:); %This T is mostly aligned with X
          lat(1) = latticeSpacing(idX);
        else
          %No TS is paralell to [1,0]
          Ts2 = Tstmp(idY,:); %This T is mostly aligned with X
          Ts1 =  [-Ts2(2), Ts2(1)]; %Define Ts1 as orthogonal to Ts2
          lat(1) = latticeSpacing(idY);
        end

        if ~isempty(idY)
          Ts2 = Tstmp(idY,:); %This T is mostly aligned with X
          lat(2) = latticeSpacing(idY);
        else
          %The points do not form a lattice but a line
          %Take a vector that is orthogonal to TS1
          Ts2 =  [-Ts1(2), Ts1(1)];
          lat(2) = latticeSpacing(idX);
        end

        Ts = [Ts1 ; Ts2];

        %Minimum and maximum coordinates of the spots
        GridSize(1,:) = min(spotCoord(:,1:2),[],1);
        GridSize(2,:) = max(spotCoord(:,1:2),[],1);
  else
      %There is only one spot. Make arbitrary lattice along the IEC axes
     lat = [1,0];
     Ts = [0,1];
     GridSize = [40,40]; %Set a lattice spacing sufficiently large so that calculateCEF uses the double integral option
  end
end

%-----------------------
% Get the list of vectoir pointing to all the neightbourgh of a specified point
%
% INPUT
% |dis2Lattice| Distance from isocentre to all other spots
% |d2iso| Distance to the nearest points surrounding the isocentre
%
% OUTPUT
% |Ts| -_SCALAR MATRIX_- |Ts(i,:)= []x,y| Unit vector pointing to one of the neighbourgh in the lattice
% |latticeSpacing| - _SCALAR_ - Distance (mm) between one spot and all its closest neighbourghs
%-----------------------
function [Ts , latticeSpacing] = getAllTs(d2iso , dis2Lattice , spotCoord , Iso)
  idxAroundIso = find(dis2Lattice == d2iso);
  LaticeApx = spotCoord(idxAroundIso , 1:2); %Coordinates of all points surrounding the isocentre
  Ts = [LaticeApx(:,1)-Iso(1) , LaticeApx(:,2)-Iso(2)]; %Vectors from isocentre to all apexes of lattice element
  nTs = sqrt(sum(Ts.^2,2));
  latticeSpacing = nTs;
  Ts = [Ts(:,1) ./ nTs , Ts(:,2) ./ nTs]; %Normalise the vector
end
