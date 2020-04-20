function [] = autoassigndocfiles( Dr, outputDir )
%AUTOASSIGNDOCFILES Autogenerates and assigns filenames for the output documentation
%     
%     [] = Dr.autoassigndocfiles( )
%     [] = Dr.autoassigndocfiles( outputDir )
% 
% Sets the `docFile` property. 
%
% If `Dr` is a Documentor object array, `autoassigndocfiles` attempts to
% recreate within `outputDir` the original subdirectory structure among the
% source files of `Dr.mFile`.
%
% If `Dr` is a single object, subdirectories aren't considered: documentation
% will be printed directly into `outputDir`.
%
% If `outputDir` is not provided as an argument, by default (and when
% possible), documentation files will be printed to a folder called "docs", one
% level above the base directory common to the paths among `Dr.mFiles`
% (presuming this base dir is a source code-specific folder, though it needn't
% be). Failing that, "docs" will be in the temporary directory returned by
% MATLAB function `tempdir()`.

    mFiles = [Dr.mFile]' ;

    %% select output directory
    if nargin < 2 || isempty(outputDir)
        % parent folder one level above the base directory common to the paths among `Dr.mFiles` 
        targetParentDir = parent( parent( mFiles, 'common'), 'single' ) ;

        [isValidPath, Values, msgId] = fileattrib( targetParentDir ) ;
        
        if Values.directory && Values.UserWrite
            outputDir = strcat( targetParentDir, filesep, "docs" ) ;
        
        elseif isfolder( tempdir )
            outputDir = strcat( tempdir, "docs" ) ;
        end
    end
    
    mustBeStringScalarOrCharVector(outputDir) ;

    %% Recreate subdirectory structure in outputDir
    [folders, names, extensions] = arrayfun( @fileparts, mFiles ) ;

    % remove baseDir component
    subDirs = erase( folders, parent( mFiles, 'common' ) ) ; 

    % ignore @class directories: class documentation goes to parent folder
    isClassDir            = contains( subDirs, filesep + "@" ) ;
    subDirs( isClassDir ) = parent( subDirs( isClassDir ), 'single' ) ;

    % assign file name(s)
    docFiles = fullfile( outputDir, subDirs, names + ".md" ) ;
    
    for iD = 1 : numel( docFiles )
        Dr(iD).docFile = docFiles(iD) ; 
    end

    if numel( unique( docFiles ) ) ~= numel( Dr ) 
        warning( [ 'Default assignments for the output documentation lead to nonunique paths.\n' ...
                   'You will need to manually reassign unique values to the `docFile` entries ' ...
                   'of the returned array prior to printing.\n' ] ) ;
    end

end
