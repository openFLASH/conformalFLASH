%% getSWversion
% Return a string with the current sofware version of FLASH MIROPT
% The version is the SHA1 hashes of the current GIT repository
%
%   * For the FLASH MIROPT repository: use the folder of getSWversion
%   * For MIROPT repositiory: use the folder of REGGUI
%   * If |scanAlgoGW| is defined, a MD5 hash obtained with the function getScanAlgoHash.m
%
%% Syntax
% |SWver = getSWversion()|
%
% |SWver = getSWversion(scanAlgoGW)|
%
%
%% Description
% |SWver = getSWversion(scanAlgoGW)| Description
%
%
%% Input arguments
% |scanAlgoGW| -_STRUCTURE_- [OPTIONAL. Default = []] Information about the scanAlgo gateway
%    * |scanAlgoGW.scanalgoGateway_IP| -_STRING_- IP address, including port, to the scanAlkgo gatewat
%    * |scanAlgoGW.room_id| -_STRING_- Room ID as defined  inthe gateway
%    * |scanAlgoGW.snout_id|  -_STRING_- snout ID as defined in the gateway
%    * |scanAlgoGW.spot_id| -_STRING_-  beam supply point as defined in the site config jar site.properties
%
%
%% Output arguments
%
% |SWver| - _CELL VECTOR of STRING_ - SHA1 of git repository of FLASH-MIROPT and MIROPT. MD5 has of the scanAlgo gateway
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function SWver = getSWversion(scanAlgoGW)

if nargin < 1
  scanAlgoGW = [];
end

%Get GIT version for FLASH-MIROPT code
folder = which ('getSWversion');
folder = fileparts(folder);
folder = fullfile(folder , '..' , '..' , '.git');
SWver{1} = readGitSHA(folder , 'FLASH-MIROPT');

%Get GIT version for MIROPT code
folder = which ('reggui');
folder = fileparts(folder);
folder = fullfile(folder , '.git');
SWver{2} = readGitSHA(folder , 'MIROPT');

if ~isempty(scanAlgoGW)
  H = getScanAlgoHash(scanAlgoGW);
  SWver{3} = ['ScanAlgo : ' H];
end


end


%------------------------------------
% REad the SHA1 code of the specified GIT repository
%------------------------------------
function SWver = readGitSHA(folder , moduleName)

  [status,result] = system (['git --git-dir ' folder ' rev-parse HEAD']);
  if status == 0
    SWver = [moduleName ' : ' strtrim(result)];
  else
    SWver = [moduleName ' : ' Unknown' ];
  end
end
