%% getRangeShifterWet
% Compute the WET  thickness of the range shifter
%
%% Syntax
% |[RangeShifterWET , R_beam , CEM_WET , NbSlabs , CEMbaseThickness , RSslabThickness , RangeShifterMaterial] = getRangeShifterWet(BDL , R_max , CEMMaterialID , MinThickness)|
%
%
%% Description
% |[RangeShifterWET , R_beam , CEM_WET , NbSlabs , CEMbaseThickness , RSslabThickness , RangeShifterMaterial] = getRangeShifterWet(BDL , R_max , CEMMaterialID , MinThickness)| Compute number of slabs of range shifter and allocate residual range shifting to base of CEM
%
% |[RangeShifterWET , R_beam , CEM_WET , NbSlabs , CEMbaseThickness , RSslabThickness , RangeShifterMaterial] = getRangeShifterWet(BDL , R_max )| Compute the range shifting for the range shifter only
%
%% Input arguments
%
% |BDL| -_STRING_- Beam data library. Name of the folder in REGGUI\plugins\openMCsquare\lib\BDL
%
% |R_max| -_SCALAR_- Range (cm) in water of the bema coming out of the range shifter
%
% |CEMMaterialID| -_STRING_- [OPTIONAL] ID of the material used for CEF as defined in MCsquare material list
%
% |MinThickness| -_SCALAR_- [OPTIONAL] Thickness (mm) in of the base of the CEM.
%
%% Output arguments
%
% |RangeShifterWET| -_SCALAR_-  WET (cm) of the range shifter to add to the CEF
%
% |R_beam| -_SCALAR_-  Range (cm) of the proton delivered by the machine
%
% |CEM_WET| -_SCALAR_-  WET (cm) of the base of the CEF
%
% |NbSlabs| -_SCALAR_- Number of slabs of range shifter
%
% |CEMbaseThickness| -_SCALAR_- Thickness (mm) of the base of the CEf to provide the residual range shifting
%
% |RSslabThickness|  -_SCALAR VECTOR_- |RSslabThickness(s)| Thickness (mm) of the s-th slab of the range shifter
%
% |RangeShifterMaterial| -_STRING_- Material of the range shifter
%
%
%% Contributors
% Authors : R. Labarbe, L. Hotoiu (open.reggui@gmail.com)

function [RangeShifterWET , R_beam , CEM_WET , NbSlabs , CEMbaseThickness , RSslabThickness , RangeShifterMaterial] = getRangeShifterWet(Plan, BDL, R_max, CEMMaterialID, MinThickness)

  if nargin < 5
    CEMMaterialID = [];
    MinThickness = [];
  end

  param = getMachineParam(BDL);

  Plan = getRangeShifterMaterialFromBDL(Plan);
  RangeShifterMaterial = Plan.Beams(1).RSinfo.RangeShifterMaterial;

  water = materialDescription('water');
  R_beam = energy2range(param.MAXenergy, water.alpha,water.p); %Range (cm) in water of the proton delivered by the machine
  EdownStr = range2energy(R_max, water.alpha,water.p); %Energy of the beam coming out of the range shifter

  %Compute the range in the range shifter material of the incoming proton and the outgoin proton
  RSThickness = getRSThickness( param.MAXenergy , EdownStr , RangeShifterMaterial); %Compute thickness (mm) of RS to reduce energy

  if ~isempty(CEMMaterialID)
    %There is a definition of CEM material
    %Let's make sure that we have an integral number of range shifter slabs
    %and that the residual range shifting is done in the CEM
    SlabSet = sort(param.snout.RSslabThickness); %Set of all available slabs
    NbSlabs = floor(RSThickness ./ SlabSet(end)); %Number of slabs of range shifter

    %Generate all sequence of range shifter slabs around the nominal value
    slabSeq = generateSlabsSequences(NbSlabs , SlabSet );

    %Identify the slab sequence that will give the thinnest CEM base
    [RS_WET , EoutRS ] = WETFromThickness(RangeShifterMaterial , param.MAXenergy , sum(slabSeq ,2)); %WET of the different sequences of slabs
    for i = 1:numel(EoutRS)
      CEMthick(i) = getRSThickness( EoutRS(i) , EdownStr , CEMMaterialID); %Compute the CEM thickness (mm) required for each slab combination
    end


    CEMthick(CEMthick < MinThickness) = Inf; %Remove the too thin CEM from the possible selections
    [~ , wSeq] = min(CEMthick); %Select the range slab sequence leading to the thinner CEM base

    RangeShifterWET = RS_WET(wSeq);
    RSslabThickness = slabSeq(wSeq,:);
    RSslabThickness(RSslabThickness==0) = []; %Remove the zeros from the list
    RSslabThickness = RSslabThickness';
    NbSlabs = numel(find(slabSeq(wSeq,:))); %Number of slabs (thick and thin)

    %Let's make the base of the CEF thicker to make the residual range shifting
    CEMbaseThickness = CEMthick(wSeq);
    CEM_WET = WETFromThickness(CEMMaterialID , EoutRS(wSeq) , CEMbaseThickness);

  else
    %There is no definition of CEM for FLASH
     %The range shifter will be supposed to be one single slab of range shifter
     RSslabThickness = RSThickness; %Thickness (mm) of the range shifter
     NbSlabs = 1;
     R_downRS = energy2range(EdownStr, water.alpha,water.p); %Range (cm) in water to reach downstream side of PTV
     RangeShifterWET = R_beam - R_downRS; %WET of the pile of range shifter slabs

     CEM_WET = 0;
     CEMbaseThickness = 0;

  end

end


%-----------------------------------------------
% Generate all the sequence of slab thicknesses around the nominal sequence
%
% INPUT
% |NbSlabs| -_INTEGER_- Nominal number of thick slabs
% |RSThick| -_SCALAR VECTOR_- Available thickness (mm) of the range shifter slabs. In increasing order.
%
% OUTPUT
% |slabSeq| -_SCALAR MARIX_- |slabSeq(sq, sl)| Thickness (mm) of the sl-th slab in the sq-th sequence
%-----------------------------------------------
function slabSeq = generateSlabsSequences(NbSlabs , RSThick)

  RSThick = sort(RSThick); %Make sure the thicknesses are ordered in increasing sequence
  RSThin = RSThick(1:end-1); %All the thin slabs
  NbThinSlb = numel(RSThin); %Number of thin slabs

  if NbSlabs > 0
      NbSequences = 1            + numel(RSThin) + NbThinSlb; %Total number os slab sequences to try
                  %Only thick     Thick slab +     replace last thick
                  %slabs          one thin         by ine thin
      slabSeq = zeros(NbSequences , NbSlabs + 1);
      slabSeq(:                            , 1:NbSlabs ) =  RSThick(end); %Place all the thick slabs
      slabSeq(2:2+NbThinSlb-1              , NbSlabs +1) =  RSThin; %Add one thin slab at the end
      slabSeq(2+NbThinSlb:2+2.*NbThinSlb-1 , NbSlabs   ) =  RSThin; %Replace last slab by a thin one
  else
      NbSequences =    1          +                 NbThinSlb; %Total number of slab sequences to try
      slabSeq = zeros(NbSequences , NbSlabs + 1);
      slabSeq(1               , 1   ) =  0     ; %Try with no slab
      slabSeq(2:2+NbThinSlb-1 , 1   ) =  RSThin; %Try one thin slab
  end



end
