%% removeDuplicates
% Remove duplicate entries in vector |a| without changing the order of elements in |a|
%
%% Syntax
% |b = removeDuplicates(a)|
%
%
%% Description
% |b = removeDuplicates(a)| Description
%
%
%% Input arguments
% |a| - _SCALAR VECTOR_ - The vector from which duplicate number shall be removed
%
%
%% Output arguments
%
% |b| - _SCALAR VECTOR_ - Vector with elements in same order as |a| but only the first occurence of each number is kept
%
%%REFERENCE
% [1] https://nl.mathworks.com/matlabcentral/answers/16667-how-to-remove-repeating-elements-from-an-array
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [b , fl ]= removeDuplicates(a , flag)

[b,m1,n1] = unique(a,'first');
[c1,d1] =sort(m1);
b = b(d1);

bf = flag(m1);
fl = bf(d1);

end
