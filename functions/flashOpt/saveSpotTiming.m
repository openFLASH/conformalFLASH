%% saveSpotTiming
% Save a text file with the timing of the delivery of the PBS spots
%
%% Syntax
% |saveSpotTiming(fileName , spotSequence , spot , Tstart , Tend)|
%
%
%% Description
% |saveSpotTiming(fileName , spotSequence , spot , Tstart , Tend)| Description
%
%
%% Input arguments
% |fileName| -_STRING_- Name of the text file where the spot time is to be saved
%
% |spotSequence| -_SCALAR VECTOR_- Order of the indices of |spot| to sort the spots. |OrderedSpot = spot(spotSequence,:)|
%
% |spot| - _SCLAR MATRIX_ - The i-th spot to deliver is spot(i,:) = [x,y]
%
% |Tstart| -_SCALAR VECTOR_- |Tstart(i)| Timing (ms) at the begining of the i-th spot delivery in |OrderedSpot|
%
% |Tend| -_SCALAR VECTOR_- |Tend(i)| Timing (ms) at the end of the i-th spot delivery in |OrderedSpot|
%
%% Output arguments
%
% None
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function saveSpotTiming(fileName , spotSequence , spot , Tstart , Tend)

  folder = fileparts(fileName);
  if (~exist(folder,'dir'))
    %The folder does not exist. Create it
    mkdir (folder)
  end

  fid = fopen(fileName , 'w');

  fprintf(fid, 'Spot Nb , X(mm) , Y(mm) , Tstart (ms) , Tend (ms)\n');
  for sptID = 1: numel(spotSequence )
      fprintf(fid, '%d , %3.2f , %3.2f , %3.2f , %3.2f \n', spotSequence(sptID) , spot(spotSequence(sptID),1), spot(spotSequence(sptID),2) , Tstart(spotSequence(sptID)) , Tend(spotSequence(sptID)));
  end

  fclose(fid);

end
