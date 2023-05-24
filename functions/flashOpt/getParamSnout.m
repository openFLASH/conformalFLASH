%% getParamSnout
% Return a description fo the parameter of the FLASH accessory holder
%
%% Syntax
% |snout = getParamSnout(snoutType)|
%
%
%% Description
% |snout = getParamSnout(snoutType)| Description
%
%
%% Input arguments
% |snoutType| - _STRING_ - Type of the FLASH accessory holder for which the parameter are requested. The options are:
%         * 'FLASH_SNOUT' : universal nozzle with a 80x80mm field size
%
%% Output arguments
%
% |snout| - _STRUCTURE_ -  Description of the properties of the FLASH accessory holder
%   * |snout.RSslabThickness| -_SCALAR VECTOR_- Available thickness (mm) of the range shifter slabs. In increasing order.
%   * |snout.RangeShifterMaterial| -_STRING_- Material of the range shifter
%   * |snout.ApertureOffset| -_SCALAR_- mm Reference position is the upstream surface of the aperture block
%   * |snout.CEMOffset| -_SCALAR_- Distance (mm) from upstram aperture surface to the upstream surface of the i-th slab of range shifter. The numbering start from the upstream slab. 1 is close to proton source. N is close to patient
%   * |snout.CEMmaxHeight| -_SCALAR_- Maximum height (mm) of the CEF to fit in the holder
%   * |snout.snout.CEMmaxRadius| -_SCALAR_- Maximum radius (mm) of the CEM to fit in the holder
%   * |snout.AccessoryCode| -_MAP_- |snout.AccessoryCode(NbSlabs)| returns the range shifter accessory code corresponding to the specified number of slabs. |NbSlabs| is a _DOUBLE_ indicating the number of slabs with |max(RSslabThickness)|
%   * |snout.RangeShifterSlabs| -_MAP_- |snout.RangeShifterSlabs(AccessoryCode)| returns a -_SCALAR VECTOR_- with the indices of the slabs for the corresponding |AccessoryCode|
%               |snout.RSslabThickness(snout.RangeShifterSlabs(AccessoryCode))| returns a -_SCALAR VECTOR_- with the thickness (mm) of each slab of the range shifter with the specified  accessory code
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function snout = getParamSnout(snoutType)
  switch snoutType
  case 'FLASH_SNOUT'
    snout.snoutType = 'FLASH_SNOUT';
    snout.RangeShifterType = 'BINARY'; % These are slabs. So its binary: slabs in or out
    snout.RSslabThickness = sort([4 , 8 , 12 , 16]); %mm Thickness of the individual slabs of the range shifter
    snout.RangeShifterOffset = [194.5 , 177.5 , 160.5 , 143.5 , 126.5 , 109.5 , 92.5 , 75.5 , 58.5 , 41.5]; % mm The numbering start from the upstream slab. 1 is close to proton source. N is close to patient
    snout.RangeShifterOffset = flip(snout.RangeShifterOffset,2); %The slab are inserted from the aperture towards the source
                      %RangeShifterOffset(i) : Distance (mm) from upstream aperture surface (= snout position) to the upstream surface of the i-th slab of range shifter

    snout.CEMOffset = 312.5; %Distance (mm) from upstream aperture surface to upstream surface of the CEF. This includes the 7mm thickness of the CEM base.
    snout.CEMmaxHeight = 95; % Maximum height (mm) of the CEF to fit in the holder
    snout.CEMmaxRadius = 109./2; % Maximum radius (mm) of the CEM to fit in the holder

    AccessoryCode = {'1 Slab'    ,'2 Slabs'    ,'3 Slabs'    ,'4 Slabs'    ,'5 Slabs'    ,'6 Slabs'    ,'7 Slabs'    ,'8 Slabs'    ,'9 Slabs'    ,'10 Slabs' , ...
                     '1.25 Slab' ,'2.25 Slabs' ,'3.25 Slabs' ,'4.25 Slabs' ,'5.25 Slabs' ,'6.25 Slabs' ,'7.25 Slabs' ,'8.25 Slabs' ,'9.25 Slabs' , ...
                     '1.5 Slab'  ,'2.5 Slabs'  ,'3.5 Slabs'  ,'4.5 Slabs'  ,'5.5 Slabs'  ,'6.5 Slabs'  ,'7.5 Slabs'  ,'8.5 Slabs'  ,'9.5 Slabs'  , ...
                     '1.75 Slab' ,'2.75 Slabs' ,'3.75 Slabs' ,'4.75 Slabs' ,'5.75 Slabs' ,'6.75 Slabs' ,'7.75 Slabs' ,'8.75 Slabs' ,'9.75 Slabs' };

    RangeShifterThickness = [1    ,2    ,3    ,4    ,5    ,6    ,7    ,8    ,9    ,10 , ...
                     1.25 ,2.25 ,3.25 ,4.25 ,5.25 ,6.25 ,7.25 ,8.25 ,9.25 , ...
                     1.5  ,2.5  ,3.5  ,4.5  ,5.5  ,6.5  ,7.5  ,8.5  ,9.5  , ...
                     1.75 ,2.75 ,3.75 ,4.75 ,5.75 ,6.75 ,7.75 ,8.75 ,9.75 ]; %Slab composition of the range shifter

    RangeShifterSlabs = {[4]  ,[4 4  ],[4 4 4  ],[4 4 4 4  ],[4 4 4 4 4  ],[4 4 4 4 4 4  ],[4 4 4 4 4 4 4  ],[4 4 4 4 4 4 4 4  ],[4 4 4 4 4 4 4 4 4  ],[4 4 4 4 4 4 4 4 4 4],...
                           [4 1],[4 4 1],[4 4 4 1],[4 4 4 4 1],[4 4 4 4 4 1],[4 4 4 4 4 4 1],[4 4 4 4 4 4 4 1],[4 4 4 4 4 4 4 4 1],[4 4 4 4 4 4 4 4 4 1],...
                           [4 2],[4 4 2],[4 4 4 2],[4 4 4 4 2],[4 4 4 4 4 2],[4 4 4 4 4 4 2],[4 4 4 4 4 4 4 2],[4 4 4 4 4 4 4 4 2],[4 4 4 4 4 4 4 4 4 2],...
                           [4 3],[4 4 3],[4 4 4 3],[4 4 4 4 3],[4 4 4 4 4 3],[4 4 4 4 4 4 3],[4 4 4 4 4 4 4 3],[4 4 4 4 4 4 4 4 3],[4 4 4 4 4 4 4 4 4 3]};

    snout.RangeShifterSlabs = containers.Map(AccessoryCode,RangeShifterSlabs);
    snout.AccessoryCode = containers.Map(RangeShifterThickness,AccessoryCode);
    snout.RangeShifterMaterial = 'aluminium';    

  otherwise
    snoutType
    error('Unknown snout model')
  end
end
