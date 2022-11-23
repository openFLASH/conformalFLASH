%% makespotList
% Make a list of spot with position and weight
%
%% Syntax
% |spot = makespotList(XposIn,YposIn,weightIn)|
%
%
%% Description
% |spot = makespotList(XposIn,YposIn,weightIn)| Description
%
%
%% Input arguments
% |XposIn| - _SCLAR VECTOR_ - |XposIn(i)| X-position (m) of the i-th row
%
% |YposIn| - _CELL VECTOR_ - |YposIn{i}(j)| Y position of the j-th spot in the i-th row
%
% |weightIn| - _CELL VECTOR_ - |weightIn{i}(j)| weight of the j-th spot in the i-th row
%
%
%% Output arguments
% |spot| - _SCLAR MATRIX_ - The i-th spot to deliver is spot(i,:) = [x,y,w ,i]
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function spot = makespotList(XposIn,YposIn,weightIn)
  index = 1;
  NbStepsX = length(XposIn);

  for Xscan = 1:NbStepsX %Nb of steps for spots
    NbStepsY = length(YposIn{Xscan}); %there can be a different number of spot in each line
    for Yscan =1:NbStepsY
      spot(index,:) = [XposIn(Xscan) , YposIn{Xscan}(Yscan) , weightIn{Xscan}(Yscan) , index];
      index = index +1;
    end
  end
end
