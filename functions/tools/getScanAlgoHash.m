%% getScanAlgoHash
% Create a MD5 hash representing the configuration of the scanAlgo gateway.
% getScanAlgoHash sends to the specified scanAlgo gateway a test set of 289 spots.
% The returning data structure is MD5 hashed into an 32 characters representing an hexadecimal number.
%
% If the config of scanAlgo is changed or if the processing algorithms of scanAlgo are changed, this will result in a different output
% and therefore a differen MD5 hash.
%
%% Syntax
% |H = getScanAlgoHash(scanAlgoGW)|
%
%
%% Description
% |H = getScanAlgoHash(scanAlgoGW)| Description
%
%
%% Input arguments
% |scanAlgoGW| -_STRUCTURE_- Information about the scanAlgo gateway
%    * |scanAlgoGW.scanalgoGateway_IP| -_STRING_- IP address, including port, to the scanAlkgo gatewat
%    * |scanAlgoGW.room_id| -_STRING_- Room ID as defined  inthe gateway
%    * |scanAlgoGW.snout_id|  -_STRING_- snout ID as defined in the gateway
%    * |scanAlgoGW.spot_id| -_STRING_-  beam supply point as defined in the site config jar site.properties
%
%
%% Output arguments
%
% |H| - _STRING_ - 32 characters representing the MD5 hash as an hexadecimal number
%
%
%%REFERENCE
% [1] https://nl.mathworks.com/matlabcentral/answers/3314-hash-function-for-matlab-struct
% [2] https://nl.mathworks.com/matlabcentral/answers/45323-how-to-calculate-hash-sum-of-a-string-using-java
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function H = getScanAlgoHash(scanAlgoGW)

  %scanalgo gateway parameters
  nb_paintings = 1;
  Energy = 226; %MeV

  %Create a test spot map
  Plan{1}.spots.nb_paintings = nb_paintings;
  Plan{1}.spots.energy = Energy;

  xv = -80:10:80;
  yv = -80:10:80;
  [X,Y] =  ndgrid(xv,yv);
  Plan{1}.spots.xy = [X(:) , Y(:)]; %pbs_convert_ScanAlgo takes IEC gantry coordinates
  Plan{1}.spots.weight = (X(:).^2 + Y(:).^2) ./100;

  %Get a data structure back from scanAlgo with this test spot map
  delivery = pbs_convert_ScanAlgo(Plan,{'nb_paintings',nb_paintings,'gateway_IP',scanAlgoGW.scanalgoGateway_IP,'room_id',scanAlgoGW.room_id,'spot_id',scanAlgoGW.spot_id,'snout_id',scanAlgoGW.snout_id,'sortSpot','false'} , false);

  %Make a hash of this spot map
  H = DataHash(delivery);
  H = sprintf('%.2x', H);   % To hex string

end

%----------------------------------------------
% Compute a hash of the structure
% Code form [1]
%
% INPUT
% |Data| -_STRUCT_- Structure to be hashed
%
% OUTPUT
% |H| -_HEX STRING_- The hexadeciaml Hash of |Data|
%----------------------------------------------
function H = DataHash(Data)
  Engine = java.security.MessageDigest.getInstance('MD5');
  %Engine = java.security.MessageDigest.getInstance('SHA-256');
  H = CoreHash(Data, Engine);

end

function H = CoreHash(Data, Engine)

  % Consider the type of empty arrays:
  S = [class(Data), sprintf('%d ', size(Data))];
  Engine.update(typecast(uint16(S(:)), 'uint8'));
  H = double(typecast(Engine.digest, 'uint8'));
  if isa(Data, 'struct')
     n = numel(Data);
     if n == 1  % Scalar struct:
        F = sort(fieldnames(Data));  % ignore order of fields
        for iField = 1:length(F)
           H = bitxor(H, CoreHash(Data.(F{iField}), Engine));
        end
     else  % Struct array:
        for iS = 1:n
           H = bitxor(H, CoreHash(Data(iS), Engine));
        end
     end
  elseif isempty(Data)
     % No further actions needed
  elseif isnumeric(Data)
     Engine.update(typecast(Data(:), 'uint8'));
     H = bitxor(H, double(typecast(Engine.digest, 'uint8')));
  elseif ischar(Data)  % Silly TYPECAST cannot handle CHAR
     Engine.update(typecast(uint16(Data(:)), 'uint8'));
     H = bitxor(H, double(typecast(Engine.digest, 'uint8')));
  elseif iscell(Data)
     for iS = 1:numel(Data)
        H = bitxor(H, CoreHash(Data{iS}, Engine));
     end
  elseif islogical(Data)
     Engine.update(typecast(uint8(Data(:)), 'uint8'));
     H = bitxor(H, double(typecast(Engine.digest, 'uint8')));
  elseif isa(Data, 'function_handle')
     H = bitxor(H, CoreHash(functions(Data), Engine));
  else
     warning(['Type of variable not considered: ', class(Data)]);
  end

end
