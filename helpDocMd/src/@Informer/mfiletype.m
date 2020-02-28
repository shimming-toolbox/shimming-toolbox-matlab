function [mType, mPath, mExist] = mfiletype( mFile )
%MFILETYPE Return type of .m file ("script"|"function"|"classdef"|"method")
% 
% MFILETYPE returns the type of .m file(s) as a string (or string vector). 
%
% ### Syntax ###
%
% [mType] = MFILETYPE( mFile )
% [mType, mPath] = MFILETYPE( mFile )
% [mType, mPath, mExist] = MFILETYPE( mFile )
% 
% ### Usage ###
%
% [mType] = MFILETYPE( mFile )
%
% mFile must be a valid file path (or set of file paths) as a string array,
% char array, or cell of character vectors (aka "cell-string").
% Possibilities for the return mType are: 
%
% - "script"   : a script
% - "function" : a function, even one with void arguments
% - "classdef" : a class definition file
% - "method"   : a class member method, i.e. any function within a folder beginning with "@"
% - "NA"       : not a MATLAB file, or not implemented (e.g. function produces an error when called by nargin())
%
% [mType, mPath] = MFILETYPE( mFile )
% 
% Returns string vector 'mPath' of the full file paths.
%
% [mType, mPath, mExist] = MFILETYPE( mFile )
%
% Returns vector of doubles 'mExist': the return values from MATLAB exist() for
% each of the files. Note that the call is made with the single path argument,
% so mExist will return "2" even for classdef files. (Also, note that exist()
% will return "2" for all of the above cases (even when the file merely ends in
% '.m' but is not a MATLAB file.) 
%
% ### References ###
%
% See also:
%
% - <https://www.mathworks.com/help/matlab/ref/exist.html exist>
%
% - This function follows the method outlined <https://blogs.mathworks.com/loren/2013/08/26/what-kind-of-matlab-file-is-this/ here>
    arguments
        mFile {mustBeStringOrCharOrCellstr} ;
    end

mPath   = abspath( strip( string( mFile ) ) ) ; % throws an error for invalid paths
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

