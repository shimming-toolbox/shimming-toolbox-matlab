classdef (Abstract) Helper
%HELPER Templates and tools for documenting code
%alt.name: GUIDEDOC? DOCASSISTANT? HELPDOC?
%
% 
% HELPER provides templates for writing help entries to new MATLAB files
%
% WIP

properties
    isCopying=false;
end
% =========================================================================
% =========================================================================    
methods

function [ Hlp ] = Helper(  )
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
    [txt] = table( varargin )
end
% =========================================================================
% =========================================================================    

end
