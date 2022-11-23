%% sigmoid
% The sigmoid responce function:
% y = 1 ./ (1 +  exp(-(x-X50).*gamma))
%
%
%% Syntax
% |y = sigmoid(x,X50,gamma)|
%
%% Description
% |y = sigmoid(x,X50,gamma)| Computes the sigmpid function
%
%% Input arguments
% |x|  - _SCALAR VECTOR_ - |x(i)| is the value of the i-th independent
%
% |X50| - _SCALAR_ - The abscisse of the 50% rise of the sigmoid function.
%
% |gamma| - _SCALAR_ - control the slope of the sigmoid function.
%
%% Output arguments
%
% |y|  - _SCALAR VECTOR_ - |y(i)| Value of the sigmoid at the i-th abcisse point
%
%
%% References
% [1] https://en.wikipedia.org/wiki/Logistic_regression
%
%% Contributors
% Authors : Rudi Labarbe (open.reggui@gmail.com)


function y = sigmoid(x,X50,gamma)
  y = 1 ./ (1 +  exp(-(x-X50).*gamma));
end
