%% scalarProd
% Scalar product of one vector Lvec with each vecotr in a matrix Ts
%
%% Syntax
% |ScalProd = scalarProd(Ts,Lvec)|
%
%
%% Description
% |ScalProd = scalarProd(Ts,Lvec)| Description
%
%
%% Input arguments
% |Ts| - SCALAR MATRIX_ - |Ts(i,:)= [x,y,,...]| The i-th vector Ts
%
% |Lvec| - SCALAR VECTOR_ - |Lvec= [x,y,,...]| Vector
%
%% Output arguments
%
% |ScalProd| - _SCALAR VECTRO_ - |ScalProd(i)| Scalar product of Ts(i,:) with Lvec
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function ScalProd = scalarProd(Ts,Lvec)

  A = Ts * repmat(Lvec',1,size(Ts,1));
  ScalProd = A(:,1);

end
