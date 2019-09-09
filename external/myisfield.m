function isFieldResult = myisfield (StructObjIn, fieldName )
%MYISFIELD
% myisfield( Obj, fieldName )
%
% Obj is the name of the structure, object, or an array of structures to search
% fieldName is the name of the field for which the function searches
% -> Returns TRUE if fieldName exists
% -> Returns FALSE otherwise
%
% from http://www.mathworks.com/matlabcentral/fileexchange/36862-viewprofiles/content/viewProfiles/myIsField.m
%
% Copyright (c) 2012, Brandon Armstrong
% All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:

%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% Updated::ryan.topfer@polymtl.ca

if ( nargin ~= 2 ) || ~( isstruct( StructObjIn ) || isobject( StructObjIn ) ) || ~ischar( fieldName )
    help(mfilename); 
    return; 
end

isFieldResult = false;

if isempty(StructObjIn)
    return
else
    fieldNames = fieldnames(StructObjIn(1));

    for i=1:length(fieldNames)
        if strcmp( fieldNames{i}, strtrim(fieldName) )
            isFieldResult = true;
            return;
        elseif isstruct( StructObjIn(1).(fieldNames{i}) )
            isFieldResult = myisfield( StructObjIn(1).(fieldNames{i}), fieldName );
            if isFieldResult
                return;
            end
        end
    end
end
