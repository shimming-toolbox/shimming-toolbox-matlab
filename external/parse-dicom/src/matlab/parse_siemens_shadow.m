function [img, ser, mrprot] = parse_siemens_shadow(varargin)
% [img, ser, mrprot] = parse_siemens_shadow(dcm, debugOutput=false)
% function to parse siemens numaris 4 shadow data
% returns three structs with image, series header, mrprot info
% does not work with arrayed dcm()
%    dependencies: parse_mrprot.m
%                  c_str.m
%                  mread.m
%   E. Auerbach, CMRR, Univ. of Minnesota, 2013

debugOutput = false;

if ((nargin < 1) || (nargin > 2))
    error('invalid input arguments');
else
    dcm = varargin{1};
    if (nargin == 2), debugOutput = varargin{2}; end
end

if (size(dcm,2) > 1)
    error('parse_siemens_shadow does not work on arrayed dicominfo data!')
end

ver_string = private_field_str_fix(dcm.Private_0029_1008);
csa_string = private_field_str_fix(dcm.Private_0029_10xx_Creator);

if (strcmp(ver_string,'IMAGE NUM 4'))
    if (strcmp(csa_string,'SIEMENS CSA HEADER'))
        img = parse_shadow_func(dcm.Private_0029_1010, debugOutput);
        ser = parse_shadow_func(dcm.Private_0029_1020, debugOutput);
    else
        error('shadow: Invalid CSA HEADER identifier: %s',csa_string);
    end
elseif (strcmp(ver_string,'SPEC NUM 4'))
    if (strcmp(csa_string,'SIEMENS CSA NON-IMAGE'))
        if isfield(dcm,'Private_0029_1210')
            img = parse_shadow_func(dcm.Private_0029_1210, debugOutput);
            ser = parse_shadow_func(dcm.Private_0029_1220, debugOutput);
        else %VB13
            img = parse_shadow_func(dcm.Private_0029_1110, debugOutput);
            ser = parse_shadow_func(dcm.Private_0029_1120, debugOutput);
        end
    else
        error('shadow: Invalid CSA HEADER identifier: %s',csa_string);
    end
else
    error('shadow: Unknown/invalid NUMARIS version: %s',ver_string);
end

% now parse the mrprotocol
if isfield(ser, 'MrPhoenixProtocol') % VB13
    MrProtocol = char(ser.MrPhoenixProtocol);
else
    MrProtocol = char(ser.MrProtocol);
end
spos = strfind(MrProtocol,'### ASCCONV BEGIN ###');
slen = 22;
if (isempty(spos)) % new VD11 format ==> ### ASCCONV BEGIN <arbitrary text> ###
    spos = strfind(MrProtocol,'### ASCCONV BEGIN');
    slen = strfind(MrProtocol(spos+3:end),'###')+6;
end
epos = strfind(MrProtocol,'### ASCCONV END ###');
if ((isempty(spos)) || (isempty(epos))), error('parse_siemens_shadow error: can''t find MrProtocol!'); end
MrProtocol = MrProtocol(spos+slen:epos-2);
mrprot = parse_mrprot(MrProtocol);

%--------------------------------------------------------------------------

function hdr = parse_shadow_func(varargin)
% internal function to parse shadow header

% input (dcm) is uint8 using little endian ordering; since this could be
% run on a little or big endian machine, we need to interpret

debugOutput = false;

if ((nargin < 1) || (nargin > 2))
    error('invalid input arguments');
else
    dcm = varargin{1};
    if (nargin == 2), debugOutput = varargin{2}; end
end

% scan through the data byte by byte

fp = 1;

[hdr_ver, fp] = mread(dcm, fp, 4, 'char');                              % version string? 4 chars 'SV10'
[unknown_str, fp] = mread(dcm, fp, 4, 'char');                          % subversion string? 4 chars '\4\3\2\1'

if (~strcmp(hdr_ver,'SV10') || ~strcmp(unknown_str,char([4 3 2 1])))
    error('this is not a recognized SV10 format header');
end

[nelem, fp] = mread(dcm, fp, 1, 'uint32-le');                           % # of elements uint32
[sig1, fp] = mread(dcm, fp, 1, 'uint32-le');                            % unknown uint32 (signature? always 77?)

if (sig1 ~= 77)
    error('unrecognized format (%d != 77 following nelem', sig1);
end

%fprintf('Found %d elements\n', nelem);

for y=1:nelem
    %data_start_pos = fp;
    [tag, fp] = mread(dcm, fp, 64, 'c_str');                            % element name tag c_str[64]
    tag = strrep(tag,'-','_'); % remove invalid chars from field name
    [vm, fp] = mread(dcm, fp, 1, 'uint32-le');                          % VM (value multiplier) uint32
    [vr, fp] = mread(dcm, fp, 4, 'c_str');                              % VR (value representation) c_str[4]
    [SyngoDT, fp] = mread(dcm, fp, 1, 'uint32-le');                     % SyngoDT uint32 (seems to just map to VR)
    [NoOfItems, fp] = mread(dcm, fp, 1, 'uint32-le');                   % NoOfItems uint32
    [sig2, fp] = mread(dcm, fp, 1, 'uint32-le');                        % unknown uint32 (always 77 or 205?)
    
    if ((sig2 ~= 77) && (sig2 ~= 205))
        error('unrecognized format (%d following NoOfItems)', sig2);
    end

    if (debugOutput), str_data = ''; end
    val_data = cell(NoOfItems);
    
    for z=1:NoOfItems
        [field_info, fp] = mread(dcm, fp, 4, 'uint32-le');              % field length info uint32[4]
        
        % of these 4 uint32: #0,#1,#3 should be the same (duplicates of field width); #2 is always 77/205
        if ((field_info(1) ~= field_info(2)) || (field_info(1) ~= field_info(4)))
            error('field width inconsistency (%d %d %d)', field_info(1), field_info(2), field_info(4));
        end
        if ((field_info(3) ~= 77) && (field_info(3) ~= 205))
            error('unrecognized format (%d following field width)', field_info(3));
        end
        
        % data field is padded to multiple of 4 chars
        field_width = field_info(1);
        field_padding = mod(4 - mod(field_width, 4), 4);
        field_alloc = field_width + field_padding;
        
        if (field_width > 0)
            [tmp_data, ~] = mread(dcm, fp, field_width, 'c_str');
            if (debugOutput)
                str_data = [str_data tmp_data]; %#ok<AGROW>
                if (z < NoOfItems), str_data = [str_data '\']; end %#ok<AGROW>
            end
            
            switch vr
                case {'AE','AS','CS','DA','DT','LO','LT','OB','OW','PN','SH','SQ','ST','TM','UI','UN','UT'}
                    % these are string values
                    %fprintf('String VR %s, data = %s\n',vr,str_data);
                    if (z == 1), val_data = cell(vm,1); end
                    val_data{z} = tmp_data;
                case {'IS','SL','SS','UL','US'}
                    % these are int/long values
                    %fprintf('%s: Int/Long VM %d, VR %s, data = %s, val = %d\n',tag,vm,vr,tmp_data,str2num(tmp_data));
                    if (z == 1), val_data = zeros(vm,1); end
                    if (size(tmp_data,2) > 0), val_data(z) = str2double(tmp_data); end
                case {'DS','FL','FD'}
                    % these are floating point values
                    %fprintf('%s: Float/double VM %d, VR %s, data = %s, val = %.8f\n',tag,vm,vr,tmp_data,str2num(tmp_data));
                    if (z == 1), val_data = zeros(vm,1); end
                    if (size(tmp_data,2) > 0), val_data(z) = str2double(tmp_data); end
                otherwise % just assume string
                    %error('Unknown VR = %s found!\n',vr);
                    %fprintf('Unknown VR %s, data = %s\n',vr,str_data);
                    if (z == 1), val_data = cell(vm,1); end
                    val_data{z} = tmp_data;
            end
        end
        
        fp = fp + field_alloc; % skip padding at the end of this field
    end

    if (debugOutput)
        fprintf('%2d - ''%s''\tVM %d, VR %s, SyngoDT %d, NoOfItems %d, Data',y-1, tag, vm, vr, SyngoDT, NoOfItems);
        if (size(str_data))
          fprintf(' ''%s''', str_data);
        end
        fprintf('\n');
    end

    hdr.(tag) = val_data;
end

%--------------------------------------------------------------------------

function outdata = private_field_str_fix(indata)
% internal function to convert uint8 data to to char when reading private
% dicom fields. this is sometimes necessary if the dicom file has passed
% through a 3rd-party dicom server which does not have the private fields
% in its dictionary. the fields are then converted to an unknown type,
% which is usually uint8 instead of char.

if (ischar(indata))
    outdata = indata;
elseif (isa(indata,'uint8'))
    outdata = char(indata');
else
    error('unexpected data type %s - should be char or uint8!', class(indata));
end
