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
%         * 'FLASH_Snout_S' : universal nozzle with a 80x80mm field size
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
  case 'FLASH_Snout_S'
    snout.snoutType = 'FLASH_Snout_S';
    snout.RangeShifterType = 'BINARY'; % These are slabs. So its binary: slabs in or out
    snout.RSslabThickness = sort([4 , 8 , 12 , 16]); %mm Thickness of the individual slabs of the range shifter
    snout.RangeShifterOffset = [194.5 , 177.5 , 160.5 , 143.5 , 126.5 , 109.5 , 92.5 , 75.5 , 58.5 , 41.5]; % mm The numbering start from the upstream slab. 1 is close to proton source. N is close to patient
    snout.RangeShifterOffset = flip(snout.RangeShifterOffset,2); %The slab are inserted from the aperture towards the source
                      %RangeShifterOffset(i) : Distance (mm) from upstream aperture surface (= snout position) to the upstream surface of the i-th slab of range shifter

    snout.CEMOffset = 312.5; %Distance (mm) from upstream aperture surface to upstream surface of the CEF. This includes the 7mm thickness of the CEM base.
    snout.CEMmaxHeight = 95; % Maximum height (mm) of the CEF to fit in the holder
    snout.CEMmaxRadius = 109./2; % Maximum radius (mm) of the CEM to fit in the holder

    AccessoryCode = {'1 Slab [Al]'    ,'2 Slabs [Al]'    ,'3 Slabs [Al]'    ,'4 Slabs [Al]'    ,'5 Slabs [Al]'    ,'6 Slabs [Al]'    ,'7 Slabs [Al]'    ,'8 Slabs [Al]'    ,'9 Slabs [Al]'    ,'10 Slabs [Al]' , ...
                     '1.25 Slab [Al]' ,'2.25 Slabs [Al]' ,'3.25 Slabs [Al]' ,'4.25 Slabs [Al]' ,'5.25 Slabs [Al]' ,'6.25 Slabs [Al]' ,'7.25 Slabs [Al]' ,'8.25 Slabs [Al]' ,'9.25 Slabs [Al]' , ...
                     '1.5 Slab [Al]'  ,'2.5 Slabs [Al]'  ,'3.5 Slabs [Al]'  ,'4.5 Slabs [Al]'  ,'5.5 Slabs [Al]'  ,'6.5 Slabs [Al]'  ,'7.5 Slabs [Al]'  ,'8.5 Slabs [Al]'  ,'9.5 Slabs [Al]'  , ...
                     '1.75 Slab [Al]' ,'2.75 Slabs [Al]' ,'3.75 Slabs [Al]' ,'4.75 Slabs [Al]' ,'5.75 Slabs [Al]' ,'6.75 Slabs [Al]' ,'7.75 Slabs [Al]' ,'8.75 Slabs [Al]' ,'9.75 Slabs [Al]' };

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

  case '40'
    snout.snoutType = '40';
    snout.RangeShifterType = 'BINARY'; % These are slabs. So its binary: slabs in or out
    snout.RSslabThickness = [0 , 50] ; %mm Thickness of the individual slabs of the range shifter
    snout.RangeShifterOffset = 0;
    snout.RangeShifterMaterial = 'PMMA';

    AccessoryCode = {'None' , 'FlashRS'};
    RangeShifterThickness = [0 , 1];
    RangeShifterSlabs = {[1], [2]};
    snout.RangeShifterSlabs = containers.Map(AccessoryCode,RangeShifterSlabs);
    snout.AccessoryCode = containers.Map(RangeShifterThickness,AccessoryCode);


  otherwise
    snoutType
    error('Unknown snout model')
  end
end
