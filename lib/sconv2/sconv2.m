function C  = sconv2(A, B, shape)
% C = sconv2(A, B, shape)
%
% Like conv2 but suitable for convolution of sparse matrices
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Date creation: 15/April/2013
%
% See also: conv2

if nargin < 3
    shape = 'full';
end

[m, n] = size(A);
[p, q] = size(B);

[i, j, a] = find(A);
[k, l, b] = find(B);

[I, K] = ndgrid(i, k);
[J, L] = ndgrid(j, l);
C = a(:)*b(:).';

switch lower(shape)
    case 'full'
        C = sparse(I(:)+K(:)-1,J(:)+L(:)-1, C(:), m+p-1, n+q-1);
    case 'valid'
        mnc = max([m-max(0,p-1),n-max(0,q-1)],0);
        i = I(:)+K(:)-p;
        j = J(:)+L(:)-q;
        b = i > 0 & i <= mnc(1) & ...
            j > 0 & j <= mnc(2);
        C = sparse(i(b), j(b), C(b), mnc(1), mnc(2)); %The function sparse ADDS the elements of the vector C if there are seveal elements assigned to the same cell i(b), j(b)
     case 'same'
        i = I(:)+K(:)-ceil((p+1)/2);
        j = J(:)+L(:)-ceil((q+1)/2);
        b = i > 0 & i <= m & ...
            j > 0 & j <= n;
        C = sparse(i(b), j(b), C(b), m, n);
end

end % sconv2
