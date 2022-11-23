%% assignResid2smallest
% Increase (or decrease) the value of the smallest elements in a sequence of integers
% so that the sum over the sequence is equal to a specified number.
%
%% Syntax
% |IntSeq = assignResid2smallest(IntSeq , TotNb)|
%
%
%% Description
% |IntSeq = assignResid2smallest(IntSeq , TotNb)| Description
%
%
%% Input arguments
% |IntSeq| -_INTEGER VECTOR_ - Sequence of integers
%
% |TotNb| -_INTEGER_- Value of the sum of integer to reach
%
%
%% Output arguments
%
% |IntSeq| -_INTEGER VECTOR_ - Updated sequence of integers so that the sum is equl to |TotNb|
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function IntSeq = assignResid2smallest(IntSeq , TotNb)

  RoundErr = TotNb - sum(IntSeq); %Difference between sum of sequence and the requested total number of elements
  Test = IntSeq;
  Test(Test==0) = NaN; %Remove zero weights
  [~, I] = sort(Test); % sort the elements of |Test| in ascending order
  idx = 1;
  while (RoundErr ~=0)
    CheckNeg = IntSeq(I(idx)) + RoundErr; %Check that there will be no negative number after removing elements from smallest value
    if (CheckNeg < 0)
      IntSeq(I(idx)) = 0; %The smallest element will be zero
      idx = idx + 1; %We will remove the residual number of elements from the next smallest value
    else
      IntSeq(I(idx)) = IntSeq(I(idx)) + RoundErr; %Adjust the smallest weight to get the correct total number of spikes
    end

    RoundErr = TotNb - sum(IntSeq); %Check whether there is a residaul number of elements to remove
  end

end
