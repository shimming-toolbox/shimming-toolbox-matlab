function [ docFile, errMsg ] = printdoc( Dr, Options )
%PRINTDOC Write documentation to file
%    
%    [docFile, errMsg] = Dr.printdoc( ) 
%    [docFile, errMsg] = Dr.printdoc( "isOverwriting", true ) 
%
% Publishes the strings assigned to `Dr.docContent` to the filepaths assigned
% to `Dr.docFile`. 
%
% To overwrite existing files, use the above name-value argument pair.
%
% Output file paths and error messages are returned as string vectors to output
% arguments 1 and 2. Entries of `errMsg` are blank "" if printing succeeded without error.
    arguments
        Dr Documentor ;
        Options.isOverwriting {mustBeBoolean} = false( size(Dr) ) ;
    end

%% Check inputs 
assert( numel(unique([Dr.docFile]')) == numel(Dr), ... 
    '"docFile" property values (i.e. paths to the printed documentation files) must be unique.' ) ;

if ( numel(Options.isOverwriting) == 1 ) && ( numel(Dr) > 1 )
    Options.isOverwriting = repmat(Options.isOverwriting, size(Dr)) ;

elseif ( numel(Options.isOverwriting) ~= numel(Dr) )
    error( ['Input value for "isOverwriting" must be a logical scalar (e.g. `true`, to apply to all Documentor entries) ' ...
            'or a logical array with the same size as the given Documentor object array.'])
end

%% Initialize outputs 
docFile = [Dr.docFile]' ;
errMsg  = strings( size(docFile) ) ;

%% Print 
for iM = 1 : numel(Dr)

    if isfile( docFile(iM) ) &&  ~Options.isOverwriting(iM)
        errMsg(iM) = ['docFile already exists. To overwrite files, use the name-value argument pair '...
                      'printdoc("isOverwriting", true). Alternatively, reassign the output paths via the `docFile` property.' ] ;
    else
        try
            docFolder = fileparts( docFile(iM) ) ;
            
            if isfolder( docFolder ) 
                [isValidPath, Values, msgId] = fileattrib( docFolder ) ;
                isWritable = Values.UserWrite ;
                if ~isWritable
                    [isWritable,errMsg,~] = fileattrib( docFolder, '+w' ) ;
                    assert( isWritable, ['Cannot write to the path assigned to `docFile`. ' ...
                        'Change the path permissions or reassign `docFile`.'] ) ;
                end
            else
                [isMade, errMsg(iM), msgId] = mkdir( docFolder ) ;
                assert( isMade, msgId, errMsg(iM) ) ;
            end    

            [fid, errMsg(iM)] = fopen( Dr(iM).docFile, 'w+' ) ;
            assert( fid~=-1, errMsg(iM) ) ;

            fprintf( fid, '%s\n', Dr(iM).docContent ) ;
            fclose(fid);

        catch Me
            errMsg(iM) = Me.message ;
        end
    end

end


end
