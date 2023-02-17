%% getPrivateTag
% Return the value of a private DICOM tag in a DICOM structure.
%
% Vendors do not reserve a private group (i.e. the first number in the DICOM tag, gggg) but a block in this group (i.e. a sequence of numbers in the second number of the tag, dddd)
% and this reservation is dynamic, that is to say that each time a DICOM application wants to encode a private attribute,
% it must first check if this block is free (are there already attributes encoded in the object to be modified or in the object being created)
% and if this is not the case, reserve the next block.
% The reservation is made by filling in a Private Creator element of an odd group gggg in the range (gggg,0010-0FF).
% For example (300D,0010) reserves the block range (300D,1000-10FF).
%
% Therefore, it is not possible to rely on the DICOM dictionry of Matlab to convert a DICOM file into a structure.
% If the private block has been move to another range than the one defined in the DICOM dictionary, then the DICOM importer of MAtlab
% add the field 'Private_gggg_xxdd' to the structure.
%
% The function |getPrivateTag| search for the tag in the structure at the field name expected by the DICOM dictionary.
% If the tag cannot be found, then it search for the provided Private Creator element.
% If it finds it, then it reads the field 'Private_gggg_xxdd' and return its value.
%
%% Syntax
% |value = getPrivateTag(gggg , PrivCreatDataElTg , PrivCreatDataElVal  , dcminfo , FieldName)|
%
%
%% Description
% |value = getPrivateTag(gggg , PrivCreatDataElTg , PrivCreatDataElVal  , dcminfo , FieldName)| The value of the private tag with name |FieldName| as defined in DICOM dictionary
%                                     This field is defined in the group |gggg|. It is contained in a private block identified by a Private Creator element with tag (|gggg| , |PrivCreatDataElTg|) in the DICOM dictionary.
%                                     In the DICOM file, the  private block identified contains the value |PrivCreatDataElVal|
%
%
%% Input arguments
% |gggg| -_STRING_- Hexadecimal number of the private DICOM group (gggg)
% |PrivCreatDataElTg| -_STRING_- Hexadecimal number of the second number of the tag (00dd) that we are searching
% |PrivCreatDataElVal| -_STRING_- Value stored in the Private Creator element to identify the private block we are tracking
% |dcminfo| -_STRUCTURE_- Structure created by the DICOM importer with the content of the DICOM file.
% |FieldName| -_STRING_- Name of the private tag, as defined in the DICOM dictionary, for which we want the value
%
%
%% Output arguments
%
% |value| - _????_ - Value of the private tag that is searched. The type of the data is based on the value representation defined in the DICOM dictionary.
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function value = getPrivateTag(gggg , PrivCreatDataElTg , PrivCreatDataElVal  , dcminfo , FieldName)

  % nameOut = dicomlookup(gggg,PrivCreaIDTag) %Get the name of the Private Creator identifier in the DICOM dictionary
  % isfield(dcminfo,nameOut)

  value = [];

  if isfield(dcminfo,FieldName)
    %The private block has the same block range as defined in the MIROPT DICOM dictionary
    %The standard Matlab DICOM importer has already properly converted the fields
    %Just read the value of the field from the structure
    value = getfield(dcminfo , FieldName);
    return
  end

  %The private block is moved to another range
  %Let's find which range
  ResvBlckNb = 0; %We do not know which block range has been reserved for our private tags
  for blck = 16:255
    PriCreatorFieldNameUp = ['Private_' upper(gggg) '_' upper(dec2hex(blck)) 'xx_Creator'];
    PriCreatorFieldNameLw = ['Private_' lower(gggg) '_' lower(dec2hex(blck)) 'xx_Creator'];

    PrvBlckName = [];
    if isfield(dcminfo , PriCreatorFieldNameUp)
      %Is there a private block present with upper case hexadecimal symbols ?
      PrvBlckName = getfield(dcminfo , PriCreatorFieldName);
      gggg = upper(gggg);
    elseif isfield(dcminfo , PriCreatorFieldNameLw)
      %Is the private block present with lower case hexadecimal symbols ?
      PrvBlckName = getfield(dcminfo , PriCreatorFieldNameLw);
      gggg = lower(gggg);
    end

    if ~isempty(PrvBlckName)
      %The private block is present
      if strcmp(PrvBlckName , PrivCreatDataElVal)
        %We found the number of our private block
        ResvBlckNb = dec2hex(blck);
        break; %stop the loop
      end
    end
  end

  if ~ResvBlckNb
    %This private block is not present in this DICOM data set
    warning (['Private block ' PrivCreatDataElVal ' is not present'])
    return
  end

  [tag , repr] = getFieldTags(upper(gggg) ,  PrivCreatDataElVal , FieldName ); %Get the tag of the field with name |FieldName|
                      %In the dictionary, the hexadecimal number use upper case

  PrivFieldName = ['Private_' gggg '_' ResvBlckNb tag(3:end)];
  if isfield(dcminfo , PrivFieldName)
    %The field is present. Let's check whether is is the one we want
    value = getfield(dcminfo , PrivFieldName);

    %Cast into the right variable type
    switch repr
      case {'SH','LO', 'LT','CS'}
        value = char(value');
      case {'FD'}
        value = typecast(value , 'double');
      case {'FL'}
        value = typecast(value , 'float');
      case {'UL'}
        value = typecast(value , 'uint32');

    end
  end

end


%--------------------------------------
% Search the DICOM dictionary for the block number of a private tag specified by
% the value stored in the Private Creator element |PrivCreatDataElVal|
%
% INPUT
% |gggg| -_STRING_- Hexadecimal number of the private DICOM group (gggg)
% |PrivCreatDataElVal| -_STRING_- Value stored in the Private Creator element to identify the private block we are tracking
% |FieldName| -_STRING_- Name of the private tag, as defined in the DICOM dictionary, for which we want the value
%
% OUTPUT
% |tag| -_STRING_- Block number |dd??| of the field |FieldName|
% |repr| -_STRING_- Value representation of the field
%--------------------------------------
function [tag , repr] = getFieldTags(gggg , PrivCreatDataElVal  , FieldName )

  fn=cell(1,4);
  [fn{:}]=textread(dicomdict('get'),repmat('%s',1,4),'delimiter','\t','commentstyle','shell');

  %Find the value reserved for our private block in the DICOM dictionary
  string =remove_bad_chars(PrivCreatDataElVal);
  IndexTag = find(cellfun(@(s) ~isempty(strfind(s,string)), fn{3})); %Index of the tag with name |PrivCreatDataElVal| in the DICOM dictionary
  privateBlock = fn{1}{IndexTag};
  privateBlock = privateBlock(9:10); %This is the value |00dd| reserved for the private block |PrivCreatDataElVal| in the DICOM dictionary

  string = [gggg , ',' , privateBlock]; %In the DICOM file, we now search for tag with this string in their number

  ListInPrvBlck = cellfun(@(s) ~isempty(strfind(s,string)), fn{1}); %List of field names in the private block
  ListWithName = cellfun(@(s) ~isempty(strfind(s,FieldName)), fn{3}); %List of field names matching the request
  IndexTag = find(ListInPrvBlck .* ListWithName); %Index of the field with a matching name and belonging to the private block

  names = split(fn{1}(IndexTag),',');
  tag = names{2}(1:end-1);

  %Identify the type of data
  repr = char(fn{2}(IndexTag));

end
