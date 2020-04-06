function [mType, mPath, mExist] = mfiletype( mFile )
%MFILETYPE Returns the type of .m file ["script"|"function"|"classdef"|"method"]
%    
%    [mType, mPath, mExist] = mfiletype( mFile )
%
% __DESCRIPTION__
% `mfiletype` checks the list of Matlab source files specified by the [string-,
% char-, or cellstr-] array of file paths `mFile` and returns:
% 
% 1. `mType` a string vector with entries
%
% `mType(i)`| when input `mFile(i)` is...  
% ----------|-----------------------------------------------------------------|
% "script"  | a script file 
% "function"| a function file (even one with void arguments)
% "classdef"| a class definition file
% "method"  | a .m function (~=constructor) in a folder beginning with "@"
%   "NA"    | an unimplemented or non-Matlab file with a .m file extension
%    ""     | an invalid file path (e.g. folder, non-Matlab, or non-existent)
% 
% 2.`mPath` a string vector of the full file system paths.
%
% 3. `mExist` a vector of return values (doubles) from Matlab function `exist` [2], 
% namely: `mExist(i) = exist( mFile(i) );` the value *should* be == 2.
%
% __NOTE__
% Although `mType(i)` will be "NA" whenever file `i` could not be assessed
% (e.g. a function that produces an error when called by `nargin`), an
% assignment other than "NA" does *not* guarantee the file is in
% working condition, and it may still contain errors.
%
% This owes to `mfiletype`'s simplistic implementation [1]: The approach
% relies entirely on MATLAB's `exist` function [2] which, as the official
% documentation states, does not check the contents of '.m' files.  
%
% [1]: (https://blogs.mathworks.com/loren/2013/08/26/what-kind-of-matlab-file-is-this/)
% [2]: (https://www.mathworks.com/help/matlab/ref/exist.html/)
%
% See also
% EXIST
    arguments
        mFile { mustBeA( mFile, ["string" "char" "cellstr"] ) } ;
    end

% pos. TODO? Enable optional keyword input (e.g. function name) and defer to it
% if abspath etc. fail to find anything assuming file path paths? ("+")

Paths  = Pathologist( mFile ) ;
mPath  = Paths.abs( mFile ) ;

mExist = Paths.exist() ;
mType  = repmat( "", size( mPath ) ) ;
mType( endsWith( mPath, ".m" ) ) = "NA" ;
mType( [isfolder(mPath) | ~isfile(mPath)] ) = "" ;

% mPath indices corresponding to .m files
iMFiles = find( [ isfile( mPath ) & endsWith( mPath, ".m" ) ] ) ;

% Note: exist( mPath(i) , "class" ) will returns 0 even when mPath(i) is a
% class file. Instead, one needs to use the class *name*. To hopefully bypass
% path precedence and namespace issues: cd to each file to ensure it has
% precedence on the Matlab path, avoiding naming conflicts. (Alt., the MATLAB
% path could be updated for each file, then reverted, but that may be slower.)
userDir = pwd ;

for iFile = 1 : nnz( iMFiles ) 
  
    iM = iMFiles( iFile ) ;

    [ folder, name, ext ] = fileparts( mPath( iM ) ) ;
    cd( folder ) ;  

    % Explicitly test classdef possibility first to bypass the other statements.
    % (i.e. exist() currently returns "2" even for classdef files when
    % called without the 2nd arg: this seems a bit odd so it might change?)
    if exist( name, 'class' )==8
        mType( iM ) = "classdef" ;

    elseif mExist( iM ) == 2
        try
            % issues error if not a function/method name
            nArg = nargin( name ) ; 
            
            % path-stripped folder name
            [~, parentName] = fileparts( folder ) ; 
                
            if startsWith( parentName, '@' ) 
            % should be sufficient, but meta.class.fromName could be used to double-check. 
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
