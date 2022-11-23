%% function
% Get the index of the Ts vector mostly aligned with Lvec
% ALso get the index of the Ts vector mostly orthogonal to Lvec
%
%% Syntax
% |idxMax = getTalignedWithL(Ts , Lvec)|
%
%
%% Description
% |idxMax = getTalignedWithL(Ts , Lvec)| Description
%
%
%% Input arguments
% |Ts| -_SCALAR MATRIX_- |Ts(i,:)= []x,y| Unit vector pointing to the i-th neighbourgh in the lattice
%
% |Lvec| -_SCALAR VECTOR_- |Lvec = [x,y]| Unit vEctor defining the axis direction
%
%% Output arguments
%
% |idxMax| - _INTEGER_ - Index of the vector Ts that is mostly paralell to Lvec  |Ts(idxMax)|
%
% |idxMin| - INTEGER_ - Index of the vector Ts that is mostly orthogonal to Lvec  |Ts(idxMax)|. If all Ts form a line, |idxMin=[]| is empty
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)


function [idxMax , idxMin]= getTalignedWithL(Ts , Lvec)

  %ScalProd = abs(sum([Ts(:,1) .* Lvec(1) , Ts(:,2) .* Lvec(2)] , 2)); %Scalar product dot(S,Ts) for each vecotr Ts
  nTs = sqrt(sum(Ts.^2,2)); %Norm of the T vectors
  Ts = Ts ./ repmat(nTs , 1 , size(Ts,2)); %Normalise the vectors Ts
  Lvec = Lvec ./ norm(Lvec);

  ScalProd = round(scalarProd(Ts,Lvec),5); %The scalar product depens only on the angle because these are normalised vectors
  [~ , idxMax] = max(abs(ScalProd), [] , 1); %Find the Ts which is mostly paralell or anti paralell to the scarf main axis

  if (numel(idxMax) > 1)
    %There are several paralell or anti paralell vector
    %Take the one with positive sign (i.e. the positive one)
    [~ , idxS] = sort(ScalProd(idxMax) , 'descend');
    idxMax = idxMax(idxS(1));
  end

  if (ScalProd(idxMax) == 0)
    % All Ts are orthogonal to Lvec
    idxMin = idxMax;
    idxMax = []; %No vector paralell to Lvec
    return
  end



  %Remove Ts that are paralell or anti-paralell to the vector we have already identified
  Tstmp = Ts;
  idxTmp = 1:size(Ts,1);
  ScalProd = round(abs(scalarProd(Tstmp,Ts(idxMax,:))),5);
  Tstmp((ScalProd == 1),:) = [];
  idxTmp(ScalProd == 1) = [];

  %Within the remaining Ts, search for the one mostly perpendicular to Lvec
  if ~isempty(Tstmp)
    ScalProd = round(scalarProd(Tstmp,Lvec),5); %The scalar product depens only on the angle because these are normalised vectors
    [~ , idxMin] = min(abs(ScalProd), [] , 1); %Find the Ts which is mostly perpendicular to the scarf main axis. We take the absolute value of the sclar product to find 0
    idxMin = idxTmp(idxMin);
  else
    %All vectors in Ts are paralell to Lvec
    idxMin = [];
  end

end
