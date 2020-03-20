function [mType, mPath, mExist] = mfiletype( mFile )
%MFILETYPE Returns the type of .m file ["script"|"function"|"classdef"|"method"]
%    
%    [mType, mPath, mExist] = mfiletype( mFile )
%
% __DESCRIPTION__
% `mfiletype` checks the list of Matlab source files specified by the [string-,
% char-, or cellstr-] array of file paths `mFile` and returns the string vector `mType`
% with the i-th element containing
%
% `mType(i)`| `mFile(i)` is...
% ----------|-----------------------------------------------------------------|
% "script"  | a script file
% "function"| a function file (even one with void arguments)
% "classdef"| a class definition file
% "method"  | a class method file (any .m file in a folder beginning with "@")
%   "NA"    | an unimplemented or non-Matlab file with a .m file extension
%    ""     | an invalid file path (e.g. folder, non-Matlab, or non-existent)
%
% **Note** `mType(i)` will be "NA" whenever the file could not be assessed
% (e.g. a function that produces an error when called by `nargin`). An
% assignment other than "NA" however does not guarantee that a file *is*
% implemented, and it may still contain errors.
%
% Additional returns: 
%
% -`mPath` a string vector containing the full file paths to the source files
%
% -`mExist` a vector of doubles returned by calling the Matlab [exist] function
% for each of the files. Note that the call is made with the single path
% argument, so `mExist` will be "2" even for classdef files.  Also, note that
% `exist()` will return "2" even when the file merely ends in '.m' but is not
% actually a MATLAB file.
%
% __SEE__
% - [exist][https://www.mathworks.com/help/matlab/ref/exist.html]
%
% `mfiletype` implements the method outlined
% [here](https://blogs.mathworks.com/loren/2013/08/26/what-kind-of-matlab-file-is-this/)
%
% See also
% EXIST
    arguments
        mFile {mustBeA{ mFile,["string" "char" "cellstr"] }} ;
    end

mPath   = abspath( strip( string( mFile ) ) ) ; 
mExist  = arrayfun( @exist, mPath ) ; 
mType   = repmat( "NA", [length(mPath) 1 ]) ;

% Saving initial working dir as it will change for each file to ensure it
% is at the top of the path to avoid potential naming conflicts. (Could update
% the MATLAB path for each file and then undo the update but this seems easier
% and presumably equivalent to, if not faster than, the alternative.)
userDir = pwd ; 

for iM = 1 : numel( mPath ) 

    [ folder, name, ext ] = fileparts( mPath( iM ) ) ;

    cd( folder ) ; % ensures file is at the top of the MATLAB path 

    % Explicitly test classdef possibility first to bypass the other statements.
    % (i.e. exist() currently returns "2" even for classdef files when
    % called without the 2nd arg: this seems a bit odd so it might change?)
    if exist( name, 'class' )==8
        mType( iM ) = "classdef" ;

    elseif mExist( iM ) == 2
        try
            nArg = nargin( name ) ; % issues error if not a function/method name
            [~,parentName] = fileparts( folder ) ; % path-stripped folder name
            
            if startsWith( parentName, '@' ) 
% should be sufficient but double-checking could be done with meta.class if necessary
                mType( iM ) = "method" ;
            else
                mType( iM ) = "function" ;
            end

        catch ME
            if strcmp( ME.identifier, 'MATLAB:nargin:isScript' ) 
                mType( iM ) = "script" ;
            end
        end
    end
end

cd( userDir ) ; % return to initial working directory

end
