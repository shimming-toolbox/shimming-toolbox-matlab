classdef (Abstract) Advisor
%ADVISOR Templates and tools for documenting code
% 
% Advisor provides templates for writing help entries to new MATLAB files
%
% WIP

properties
    isCopying=false;
end
% =========================================================================
% =========================================================================    
methods

function [ Ad ] = Advisor(  )
    return ;
end

end
% =========================================================================
% =========================================================================    
% methods( Static )
methods( Static )

function [  ] = template(  )
%TEMPLATE Templates for documentting new .m files
%
% should return (more importantly) + *define* a template for .m help documentation:
% for basic scripts, functions, class definitions, etc. 
%
% should also be able print to a new file.m (w/out overwrite option)
% 
% DEFAULT.link[1]: http://en.wikipedia.org/wiki/Markdown        "Markdown"
DEFAULT_SECTIONS = { 'NAME', 'SYNTAX', 'DESCRIPTION', 'INPUTS', 'OUTPUTS', 'EXAMPLES', 'NOTES', 'SEE ALSO' } ;
%?

% FUNCTION_NAME brief summary (:
%     
%    [paths, List] = findfiles( sFolder, sPattern, isRecursive, isExcludingHidden, returnType )
% 
% Looks for files and/or subfolders in `sFolder` with names matching `sPattern`
% by calling Matlab function [dir] and returns the file paths as elements of a
% string column vector `paths`. The single-element structs output by `dir()`
% are arrayed and returned as `List.
% 
% __INPUTS__ 
%     
%   sFolder=["."]
%     The base directory of the search.    
%
%   sPattern=["*.*"] 
%     The searchPattern of interest. If provided as a string array, patterns
%     are searched successively. (The default corresponds to including  all
%     files with explicit with explicit file extensions.)
%
%   isRecursive=[true|1]
%     Toggle to include (1) or exclude (0) subdirectories in the search.
%
%   isExcludingHidden=[true|1]
%     Toggle to include (1) or exclude (0) hidden files (i.e. for
%     Unix: filenames beginning with ".")
%
%   returnType=["files"]
%     Selects what type of path elements are retained in the two outputs : 
%     Options are: "files", "folders", or "both". 
% 
% ETC 
%
%   For more info, refer to the documentation for
%   [dir](https://www.mathworks.com/help/matlab/ref/dir.html)
%
% See also DIR

end


function [  ] = recommendations(  )
%RECOMMENDATIONS On technical writing 
%
% ... 

end

function [  ] = references(  )
%REFERENCES Suggested references 
%
% ...

end

end
% =========================================================================
% =========================================================================    
methods( Static )
    %.....
    [txt] = readme( isCopying, S )
    %.....
    [txt] = txttable( varargin )
    %..... 
    [txt] = function_template( pathIn )
end
% =========================================================================
% =========================================================================    

end
