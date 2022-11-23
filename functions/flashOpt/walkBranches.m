%% walkBranches
% Order the PBS spots by walking each branch sequentially.
% When walking on branch, scan the spots along "rib" line perpendicualr to the spine of the branch
% Sucessive ribs are scanned in opposite direction
% Branches are scanned in the order defined by |BranchScanOrder|
% If a branch ends at a branching point, then the next branch is the shortest branch reaching this branching point
%
%% Syntax
% |OrderedIdx = walkBranches(spotCoord , branches , ScarfHalthWidth , BranchScanOrder)|
%
%
%% Description
% |OrderedIdx = walkBranches(spotCoord , branches , ScarfHalthWidth , BranchScanOrder)| Description
%
%
%% Input arguments
% |spotCoord| - _SCLAR MATRIX_ - The i-th spot to deliver is spot(i,:) = [x,y,k] where k is the spot index in the calling function
%
% |branches| -_STRUCTURE_- Information about the skeleton of the scanned area
%  * |branches.BranchLength| -_SCALAR VECTOR_- |branches.BranchLength(b)| Length of the b-th branch
%  * |branches.branchingPoint| -_SCALAR VECTOR_- |branches.branchingPoint(i)=p| The i-th branch is connected to the p-th branching point. p=0 means that the branch is connected to no branching point
%  * |branches.BranchingCoord| -_SCALAR MATRIX_- |branches.BranchingCoord(p)=[x,y]| The (x,y) coordinate of the p-th branching point
%  * |branches.Path| -_CELL VECTOR_- |branches.Path{b}(j) = [x,y]| The j-th point of the b-th branch has coordinates [x,y]
%                     The point are ordered so that sucessesive point on the line have indices b and b+1. b=0 is one end point. b=end is a branching point or another end point
%
% |ScarfHalthWidth| -_SCALAR_- Width (mm) along the halft width of a scarf pattern
%
% |BranchScanOrder| -_SCALAR VECTOR_- Order in which the branches are selected
%
%
%% Output arguments
%
% |OrderedIdx| -_SCALAR VECTOR_- |OrderedIdx(i)| = N indicates that the i-th spot to deliver is the spot spotCoord(N,:)
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function OrderedIdx = walkBranches(spotCoord , branches , ScarfHalthWidth , BranchScanOrder)

  DistMat = interSpotDistance(spotCoord);

  %Get information about the spot lattice size and spacing
  [Lat , Ts , GridSize] = getLatticeInfo(spotCoord(:,1:2));
  AlreadyScanned = zeros(numel(branches.BranchLength),1); %Keep track of the branches we already scanned
  BrOrder = [];
  [Bidx , BranchScanOrder , T , S , latticeSpacing] = pickNextBranch(BranchScanOrder , branches , Ts , Lat);

  % figure(100)
  % hold on
  % plot(branches.Path{Bidx}(:,1),branches.Path{Bidx}(:,2),'r-^')


  %Scan all branches sucessively
  %-----------------------------
  invOrder = false;
  direc = 1; %Define the scanning direction along a rib
  OrderedIdx = [];
  TurnPt = []; %Flag points which are at the end of ribs

  [idxA , idxB]= getTalignedWithL(Ts , [1,0]);
  Ta = Ts(idxA,:)  .* Lat(idxA); %Lattice vector mostly paralell to X IEC gantry axis. Length is equal to the lattice period
  if isempty(idxB)
    %The point do not form a lattice but a line
    %Take a vector that is orthogonal to Ta
    Tb = [-Ta(2), Ta(1)];
  else
    Tb = Ts(idxB,:) .* Lat(idxB); %Lattice vector mostly paralell to Y IEC gantry axis
  end

  [~ , isoIdx] = getDistances2(spotCoord , [0,0]); %Find the point closest to isocentre
  C = spotCoord(isoIdx(1),1:2); %This is a point on the grid that is the closest to the scarf axis


  while (~prod(AlreadyScanned)) %Loop untill all branches are scanned
    for pts = 1:branches.BranchLength(Bidx)
        switch invOrder
            case false
            %Walk from end point
            ptsIdx = pts;
            case true
            %Walk from branching point
            ptsIdx = branches.BranchLength(Bidx) - pts + 1 ;
        end
        %Loop for every point on the branch
        P1 = branches.Path{Bidx}(ptsIdx,:);

        %Jump to the closest point on the spot lattice
        X = linsolve([Ta' , Tb'] , (P1 - C)');
        X = round(X);
        P1 = C + X(1) .* Ta + X(2) .* Tb;

        %Scan the rib
        seq = OrderTherib(P1 , T , direc , GridSize , spotCoord(:,1:2) , ScarfHalthWidth , latticeSpacing); %Order the spots along one rib at position P1. Scan along vector T
        if ~isempty(seq)
          tp = zeros(1,size(seq,2));
          tp(1) = -1; %Flag -1 at the index of the begining of a rib
          tp(end) = +1; %Flag +1 at the index of the end of a rib. If the rib has a single spot, then there will only be a +1 flag for this rib
          OrderedIdx = [OrderedIdx , seq];
          TurnPt = [TurnPt , tp];
        end

        direc = direc .* (-1); %Change the scanning direction along the rib. Sucessive ribs are scanned head to tail
    end %for pts

    AlreadyScanned(Bidx) = 1; %The branch Bidx is now scanned

    %Choose a new branch. Take the next longest one
    len = branches.BranchLength;
    if (~prod(AlreadyScanned))
        %There are still unscanned branches
        [Bidx , BranchScanOrder , T , S, latticeSpacing] = pickNextBranch(BranchScanOrder , branches , Ts, Lat);
        invOrder = false;
    end
  end

  %Remove duplicate spots in the list
  %----------------------------------
  [OrderedIdx , TurnPt ]= removeDuplicates(OrderedIdx , TurnPt);

  %Check for missing spots
  %-----------------------
  unscanned = spotCoord(:,1:2);
  indx = 1:size(unscanned,1);
  unscanned = [unscanned , indx'];
  unscanned(OrderedIdx,:)=[]; %List of the spot that were missed by the scarves

  if (~isempty(unscanned))
    for idx = 1: size(unscanned,1)
      scanned = spotCoord(OrderedIdx,1:2); %Update List of the olready scanned spots
      d  = getDistances2(scanned, unscanned(idx,1:2)); %Find the closest spot in the already scanned list
      d_filer = d' .* (TurnPt ~=0);
      d_filer(d_filer==0)= NaN;
      [dmin , posInSequence] = min(d_filer, [] , 'omitnan'); %Find the closest spot that is at the begining or end of a rib
      switch (TurnPt(posInSequence))
        case -1
            OrderedIdx = [OrderedIdx(1:posInSequence-1) , round(unscanned(idx,3)) , OrderedIdx(posInSequence:end)]; %Insert spot in the scanning sequence
            TurnPt = [TurnPt(1:posInSequence-1) , 0 , TurnPt(posInSequence:end)];

        case 1
            OrderedIdx = [OrderedIdx(1:posInSequence) , round(unscanned(idx,3)) , OrderedIdx(posInSequence+1:end)]; %Insert spot in the scanning sequence
            TurnPt = [TurnPt(1:posInSequence) , 0 , TurnPt(posInSequence+1:end)];
      end

    end
  end
end

%==================================
% Find a vector T nearly orthogonal to the spine (P2-P1) but that points toward
% one of the closest neighbourgh of the P1 point surch that
%
% INPUT
% |P1| -_SCALAR VECTOR_- |P1=[x,y]| First point of the spine vector
% |P2| -_SCALAR VECTOR_- |P2=[x,y]| Second point of the spine vector
% |spotCoord| -_SCALAR MATRIX_- |spotCoord(i,:)=[]x,y| The coordinate of the i-th PBS spots
% |DistMat| - _SCALAR MATRIX_ - DistMat(i,j) distance between the i-th and j-th spot
%
% OUTPUT
% |T| -_SCALAR VECTOR_- Vector roughthly orthogonal to (P2-P1) and pointing to the closest neightbourgh to P1
% |S| -_SCALAR VECTOR_- VEctor pointing along the spine from P1 to P2
%==================================
function [T , S , indx] = findVec2Neighbor(P1, P2, Ts)

      S = P2 - P1; %Vector tangent to the branch at point P1 and pointing towards P2
      S = S ./ norm(S); %Unit vector pointing along the spine

      [~ , indx]= getTalignedWithL(Ts , S);
      T = Ts(indx , :); %Vector nearly normal to spine and pointing along a line of the lattice

end


%============================
% Order the spot along a rib
% INPUT
% |P1| -_SCALAR VECTOR_- |P1= [x,y]| Index coordinate of the point on the spine from which the rib is scanned
%
% |T| -_SCALAR VECTOR_- |T= [x,y]| Vector pointing from P1 to the closest neighbourgh
%
% |direc| -_SCALAR_- +1 or -1 defines the scanning direction along the rib
%
% |mapImage| -_SCALAR MATRIX_- mapImage(x,y) = N if the N-th spot is a spot delivered at location (x,y)
%
%
% |ScarfHalthWidth| -_INTEGER_- %Number of spots along the halft width of a scarf pattern
%
% OUTPUT
% |seq| -_SCALAR VECTOR_- |seq(i) = N| The i-th spot to be scanned is labelled N in the |mapImage|

%============================
function seq = OrderTherib(P1 , T , direc , GridSize , spotCoord  , ScarfHalthWidth , latticeSpacing)

  seq = [];
  Nstep = round(ScarfHalthWidth ./ latticeSpacing); %Number of spot in a halkf scarf width

  %Deal with the case when the scarf is one spot wide
  wSingleLine = find(diff(GridSize,1)==0); %Along this axis, the scarf is one spot wide
  if (~isempty(wSingleLine))
    GridSize(1,wSingleLine)= GridSize(1,wSingleLine) - latticeSpacing;
    GridSize(2,wSingleLine)= GridSize(2,wSingleLine) + latticeSpacing;
  end


  for rib = -Nstep:1:Nstep
    %Search for spot along one direction of the rib
    Pt = P1 + direc .* T .* rib .* latticeSpacing; %Coordinate of a point along the rib

    % figure(100)
    % hold on
    % plot(Pt(1),Pt(2),'+')


    if (isPtinGrd(Pt,GridSize .* 1.1))
        [d , SpotIndex] = getDistances2(spotCoord , Pt);
        if (sqrt(d(SpotIndex(1))) < latticeSpacing./2)
          %There is a spot delivered close to the selected spot
          seq = [seq , SpotIndex(1)]; %Add the spot to the end of the sequence
          % figure(100)
          % hold on
          % text(Pt(1),Pt(2),num2str(SpotIndex(1)))
        end

        %pause
    end
  end

  seq = unique(seq,'stable'); %Make sure that the spots are counted only once in the sequence. Do not upset the order of the spots

end

%======================
%Scan the branch in the order defined in BranchScanOrder
%======================
function [Bidx , BranchScanOrder , T , S , latticeSpacing] = pickNextBranch(BranchScanOrder , branches , Ts , Lat)
  Bidx = BranchScanOrder(1);
  BranchScanOrder(1) = [];

  P1 = branches.Path{Bidx}(1,:);
  P2 = branches.Path{Bidx}(end,:); %The branch is a straight line. S is the first and last vector
  [T , S, indx] = findVec2Neighbor(P1, P2, Ts);
  latticeSpacing = Lat(indx);

end
