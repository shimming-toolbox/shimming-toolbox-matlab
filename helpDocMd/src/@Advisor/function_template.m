function [txt] = function_template( )
% FUNCTION_TEMPLATE Provides a template for function help entries
%     
%     [txt] = function_template( )
%
% `function_template` returns the template as a string vector `txt`.
%
% 
%

txt = [] ;
% should also be able print to a new file.m (w/out overwrite option)
% 
% DEFAULT.link[1]: http://en.wikipedia.org/wiki/Markdown        "Markdown"
% DEFAULT_SECTIONS = { 'NAME', 'SYNTAX', 'DESCRIPTION', 'INPUTS', 'OUTPUTS', 'EXAMPLES', 'NOTES', 'SEE ALSO' } ;
%?

% FUNCTION_NAME **Brief summary goes here. <80 char incl. function name** (See Note [1])
% 
% __SYNTAX__ <-*(MATLAB convention for the 1st header title. See [Note 2])*
% >> 4x spaces from margin + 1 empty line formats the following `as code` in [Markdown].
% ....(See [Note 3]) 
%    
%    [out1] = FUNCTION_NAME( mandatoryInput1, mandIn2_A )
%    [out1] = FUNCTION_NAME( mandatoryInput1, mandIn2_B )
%    [out1] = FUNCTION_NAME( mandatoryInput1, mandIn2_C, optionalInput1 )
%    [out2] = FUNCTION_NAME( mandatoryInput1, mandIn2_C, optionalInput1, optionalInput2 )
% 
% __DESCRIPTION__ <-*(2x underscores format to bold in Markdown, but italics in MATLAB Markup...)*
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
% __NOTES__ 
%
% [1] MATLAB convention is to use ALL CAPS for all instances of the function name.
% If the user queries `help` running the MATLAB GUI, FUNCTION_NAME is automatically reformated to the
% proper casePersonally, I would recommend not using them except maybe for the
% very first
% instance because it stands to mislead anyone who happens to be reading it
% from while running MATLAB in a terminal session. (.
%
% The leading "H1" line is special as it is searchable via the [lookfor] function, so
% make sure it contains relevant terms! Also, keep it brief or it risks being
% truncated when displayed by `lookfor`. 
%
% Also note, MATLAB convention is to use
% ALL CAPS for the function name. 
%
% [H1](https://www.mathworks.com/help/matlab/matlab_prog/add-help-for-your-program.html)
% [lookfor](https://www.mathworks.com/help/matlab/ref/lookfor.html)
%
% [Markdown](https://daringfireball.net/projects/markdown/syntax)
%   For more info, refer to the documentation for
%   [dir](https://www.mathworks.com/help/matlab/ref/dir.html)
% [UNIX](https://en.wikipedia.org/wiki/Man_page) favors "SYNOPSIS", others use "USAGE"
% See also 
% LOOKFOR

end

