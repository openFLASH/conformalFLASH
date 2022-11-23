%% getNozzle2IsoDistanceFromBDL
% Return the distance from isocentre to the nozzle exit of the nozzle defined in the beam data model.
%
%% Syntax
% |iso2Nozzle = getNozzle2IsoDistanceFromBDL(BDL_file)|
%
%
%% Description
% |iso2Nozzle = getNozzle2IsoDistanceFromBDL(BDL_file)| Return the nozzle exit to isocentre distance for the specified beam data model
%
%
%% Input arguments
% |BDL_file| - _STRING_ -  File name of the file containing the beam data model.  The file name must contain the path to the file or be located in a directory contained in 'path'.
%
%% Output arguments
%
% |iso2Nozzle| - _SCALAR_ - Distance (mm) from the nozzle exit to the isocentre
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function iso2Nozzle = getNozzle2IsoDistanceFromBDL(BDL_file)

fid = fopen(BDL_file);

  while 1
      tline = fgets(fid);
      if ~ischar(tline), break, end
      if (not(isempty(strfind(tline,'Nozzle exit to Isocenter distance'))))
          tline = fgetl(fid);
          eval(['iso2Nozzle = ',tline,';']);
      end
  end

end
