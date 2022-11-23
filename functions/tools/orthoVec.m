%%orthoVec
% Find the vector b orthogonal to the vector a in 2D
%
%% Syntax
% |b = orthoVec(a)|
%
%
%% Description
% |b = orthoVec(a)| Description
%
%
%% Input arguments
% |a| - _SCALAR VECTOR_ - |a=[x,y]|  Components of the vector
%
%
%% Output arguments
%
% |b| - _SCALAR VECTOR_ - |b=[x,y]|  Components of the vector normal to |a|
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function b = orthoVec(a)

  if(sum(a.^2)==0)
    error('No orthogonal vector')
  end
  if a(2) ~= 0
    b(1)=1;
    b(2)= - b(1) .* a(1) ./a(2);
    b = b ./norm(b);
  else
    b(2)=1;
    b(1)= - b(2) .* a(2) ./a(1);
    b = b ./norm(b);
  end

end
