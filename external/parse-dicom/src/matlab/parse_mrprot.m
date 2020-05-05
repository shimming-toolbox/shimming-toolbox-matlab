function mrprot = parse_mrprot(arr)
% parse mrprot (text) structure
%   E. Auerbach, CMRR, Univ. of Minnesota, 2013

% isolate the meas.asc portion only (strip the meas.evp part if present)
if (size(arr,1) > 1), arr = arr'; end
spos = strfind(arr,'### ASCCONV BEGIN ###');
epos = strfind(arr,'### ASCCONV END ###');
if (~isempty(spos) && ~isempty(epos))
    arr = arr(spos:epos-1);
end

% mgetl returns a cell array where each cell is a line of text
[lines, numlines] = mgetl(arr);

mrprot = [];
for curline=1:numlines
    % work on one line at a time
    line = lines{curline};
    
    % strip any comments
    ctest = strfind(line, '#');
    if (ctest > 1)
        line = line(1:ctest-1);
    end
    
    if ( size(line,2) < 1 )
        % blank line -- skip
    elseif ( ~isletter(line) | ~isletter(line(1)) )  %#ok<OR2> % note: leave as | (not scalar)
        % blank line or comment -- skip
    else
        skip_this = false;
        
        % first, take the variable name -- everything before the '='
        [varname, stub] = strtok(line, '=');
        
        % break it down into levels delimited by '.'
        varnamearr = textscan(varname, '%s', 'delimiter', '.');
        varnamearr = varnamearr{:}';
        lvls = size(varnamearr,2);
        
        index = zeros(lvls,2);
        for x=1:lvls
            % look for brackets indicating array
            workstr = strtrim(varnamearr{x});

            btest1 = strfind(workstr, '[');
            btest2 = strfind(workstr, '][');
            
            if (btest2 > 0)         % found a 2D array
                btest1 = btest1(1);
                index(x,1) = str2double(getBracketString(workstr));
                index(x,2) = str2double(getBracketString(workstr(btest2+1:end)));
                varnamearr{x} = [workstr(1:btest1-1) blanks(size(workstr,2) - btest1 + 1)];
            elseif (btest1 > 0)     % found a 1D array
                index(x,1) = str2double(getBracketString(workstr));
                index(x,2) = 0;
                varnamearr{x} = [workstr(1:btest1-1) blanks(size(workstr,2) - btest1 + 1)];
            else
                index(x,1:2) = [0 0];
            end
        end
        
        fields = {{1}};
        for x=1:lvls
            varname = make_safe_fieldname(varnamearr{x});
            
            if (strncmp(varname, 'x__attribute__', 14)) % new for VD11 (note: after VD13 mod, make_safe_fieldname() prepends 'x')
                % skip it
                skip_this = true;
                break;
            end
            
            if (strncmp(varname, 'm_', 2)) % fix VAxx YAPS naming
                varname = varname(3:end);
                old_version = true;
                
                %fix VAxx array naming, since they are equivalent to
                %VBxx non-arrayed parameters
                if ( (strncmp(varname, 'al', 2)) || (strncmp(varname, 'afl', 2)) )
                    test_varname = varname(2:end);
                else
                    test_varname = varname;
                end
            else
                old_version = false;
                test_varname = varname;
            end
            
            % these changed in VD11 for some reason; rename for backward compatibility
            if (strncmp(varname, 'sWipMemBlock', 12))
                varname = 'sWiPMemBlock';
            end
            
            idx1 = index(x,1) + 1;
            idx2 = index(x,2) + 1;
            
            if (x < lvls)
                fields = {fields{:}, varname, {idx1, idx2}};
            else
                % last one
                stub = getYaPSBracketString(getEqualString(stub));
                
                if ( (strncmp(test_varname, 'afl', 3)) || ... % array of floats (numeric)
                        (strncmp(test_varname, 'ac', 2))  || ... % array of signed chars (?) (numeric)
                        (strncmp(test_varname, 'ad', 2))  || ... % array of doubles (numeric)
                        (strncmp(test_varname, 'al', 2))  || ... % array of longs (numeric)
                        (strncmp(test_varname, 'ax', 2))  || ... % array of complex values (numeric)
                        (strncmp(test_varname, 'MaxOnlineTxAmpl', 15))  || ... % array of floats
                        (strncmp(test_varname, 'MaxOfflineTxAmpl', 16))  || ... % array of floats
                        (strncmp(test_varname, 'IdPart', 6))  || ... % array of bytes (?)
                        (strncmp(test_varname, 'an', 2)) )      % array of signed value (long?) (numeric)
                    fields = {fields{:}, varname, {idx1, idx2}, str2num(stub)};  %#ok<*CCAT,ST2NM> % (use str2num here for arrays)
                    
                elseif ( (strncmp(test_varname, 'aui', 3)) || ... % array of unsigned int (hex)
                        (strncmp(test_varname, 'aul', 3)) )     % array of unsigned long (hex)
                    fields = {fields{:}, varname, {idx1, idx2}, getHexVal(stub)};
                    
                elseif (strncmp(test_varname, 'at', 2)) % array of text strings
                    fields = {fields{:}, varname, {idx1, idx2}, cellstr(getQuotString(stub))};
                    
                elseif ( (strncmp(test_varname, 'fl', 2)) || ... % float (numeric)
                        (strncmp(test_varname, 'l', 1))  || ... % long value (signed) (numeric)
                        (strncmp(test_varname, 'd', 1))  || ... % double value (numeric)
                        (strncmp(test_varname, 'n', 1))  || ... % signed value (long?) (numeric)
                        (strncmp(test_varname, 'b', 1))  || ... % bool (numeric)
                        (strncmp(test_varname, 'e', 1))  || ... % enum (long?) (numeric)
                        (strncmp(test_varname, 'WorstCase', 9))  || ... % float
                        (strncmp(test_varname, 'Nucleus', 7))  || ... % float
                        (strncmp(test_varname, 'BCC', 3))  || ... % int (numeric)
                        (strncmp(test_varname, 'WaitForUserStart', 16))  || ... % bool (numeric)
                        (strncmp(test_varname, 'i', 1)) )      % int (numeric)
                    fields = {fields{:}, varname, str2double(stub)};
                    
                elseif ( (strncmp(test_varname, 'ush', 3)) || ... % unsigned short (hex)
                        (strncmp(test_varname, 'ul', 2))  || ... % unsigned long (hex)
                        (strncmp(test_varname, 'ui', 2))  || ... % unsigned int (hex)
                        (strncmp(test_varname, 'un', 2))  || ... % unsigned value (hex)
                        (strncmp(test_varname, 'UseDouble', 9))  || ... % unsigned value (hex)
                        (strncmp(test_varname, 'Ignore', 6))  || ... % unsigned value (hex)
                        (strncmp(test_varname, 'uc', 2)) )      % unsigned char (hex)
                    fields = {fields{:}, varname, getHexVal(stub)};
                    
                elseif ( (strncmp(test_varname, 't', 1)) || ... % text string
                         (strncmp(test_varname, 's', 1)) )
                    fields = {fields{:}, varname, getQuotString(stub)};
                    
                elseif ( (old_version) && (strcmp(test_varname, 'SecUserToken')) )
                    % don't know what this is
                    fields = {fields{:}, varname, stub};
                    
                else
                    skip_this = true;
                    % 20180501::ryan.topfer@polymtl.ca changed:
                    % fprintf('parse_mrprot WARNING: unknown data type for %s (value = %s), discarding this line!!',varname, stub);
                    % in order to suppress the constant warning messages from within dicominfosiemens.m
                    % The format of the displayed warning changes slightly from
                    % the original, but by default it will still appear unless
                    % dicominfosiemens() is called first (i.e. since
                    % these warnings might still be useful -- better not to
                    % suppress them in general!)
                    warningStr = sprintf('parse_mrprot WARNING: unknown data type for %s (value = %s), discarding this line!!',varname, stub);
                    warningId  = 'PARSE_MRPROT:unknownDataType' ;
                    warning( warningId, warningStr ) ;
                end
            end
        end
        
        if (~skip_this), mrprot = setfield(mrprot,fields{:}); end
    end
end



%--------------------------------------------------------------------------

function stvar = getEqualString(text)
% strips leading whitespace and '=' character to find assignment string

stvar = text(strfind(text,'=')+1:end);
stvar = strtrim(stvar);

%--------------------------------------------------------------------------

function stvar = getYaPSBracketString(text)
% extracts text string from YAPS output in meas.asc
% sometimes brackets [] are used (old syngo), sometimes not (new syngo)

a = strfind(text,'[');
b = strfind(text,']');

if ( isempty(a) || isempty(b) )
    % no brackets found (>=VB12T)
    a = strfind(text,'=');
    if (a > 0) % found = sign
        stvar = strtrim(text(a+1:end));
    else
        stvar = text;
    end
elseif ((b - a) < 2)
    stvar = '';
else
    stvar = text(a+1:b-1);
end


%--------------------------------------------------------------------------

function stvar = getBracketString(text)
% extracts text string from within [] brackets

a = strfind(text,'[');
b = strfind(text,']');

if ( isempty(a) || isempty(b) )
    stvar = '';
elseif ((b - a) < 2)
    stvar = '';
else
    stvar = text(a+1:b-1);
end

%--------------------------------------------------------------------------

function stvar = getQuotString(text)
% extracts string between double quotes, e.g. "string"
%  also works with double-double quotes, e.g. ""string""

idx = strfind(text,'"');

if ( (length(idx) == 4) && (idx(1)+1 == idx(2)) && (idx(3)+1 == idx(4)) ) % double-double quotes
    stvar = text(idx(2)+1:idx(3)-1);
elseif (length(idx) >= 2) % double quotes, or ??? just extract between first and last quotes
    stvar = text(idx(1)+1:idx(end)-1);
else % malformed?
    stvar = text;
end

%--------------------------------------------------------------------------

function hval = getHexVal(text)
% gets C++ style hexadecimal value from a text string
%  (assumes nothing follows the hex string)

hexst = strfind(text,'0x');
if (~isempty(hexst))
    tmp = text(strfind(text,'0x')+2:end);
    hval = hex2dec(tmp);
else
    % in VD11 the hex notation is inconsistent...sometimes these are dec
    hval = str2double(strtrim(text));
end

%--------------------------------------------------------------------------

function [lines, numlines] = mgetl(arr)
% mgetl: parse an entire text file into cell array of strings where each
% cell is one line. recognizes dos and unix file formats.

arr = char(arr);
if (size(arr,1) > size(arr,2)), arr = arr'; end
arr_len = length(arr);
lf = strfind(arr, char(10));
numlines = length(lf);
if (lf(numlines) < arr_len)
    numlines = numlines + 1;
    lf(numlines) = arr_len;
end

lines = cell(numlines, 1);

stpos = 1;
for x=1:numlines
    if (x > 1), stpos = lf(x-1)+1; end
    endpos = lf(x) - 1;
    if (lf(x) > 1)
        if (arr(endpos) == char(13))
            endpos = endpos - 1; % strip cr and lf if both present
        end
    end
    if (endpos >= stpos), lines{x} = arr(stpos:endpos); end
end

%--------------------------------------------------------------------------

function stvar = make_safe_fieldname(tagStr)
% this function checks potential fieldnames and makes sure they are valid
% for MATLAB syntax, e.g. 2DInterpolation -> x2DInterpolation (must begin
% with a letter)

tagStr = strtrim(tagStr);

if (isletter(tagStr(1)))
    stvar = tagStr;
else
    stvar = strcat('x', tagStr);
end

if strfind(stvar, ';'), stvar = strrep(stvar, ';', '_'); end
if strfind(stvar, '@'), stvar = strrep(stvar, '@', '_'); end % VD13
