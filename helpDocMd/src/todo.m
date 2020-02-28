% TODO
% 
% 1. Parsing functions' and methods' arguments blocks to enable auto-documenting of
% input constraints (size, type, validation functions) and default assignments.
% (for defaults, it might be more sensible generally to copy the assignment
% portion of the code itself (e.g. "= zeros(100,100)") rather than the value itself.)
%   
%   Notes: 
%   - (See related item (9?) down the list re: adding the option to include
%   implementation details in the documentation)
%
%   - nargin and nargout are readily obtained for functions + methods, but
%   getting the remaining details will require something more involved (e.g.
%   parse return from a call to MATLAB's TEXTSCAN()) 
%   (probably not that difficult given the arguments block needs to fall between
%   the leading documentation comments and the beginning of the function code.)
%
%   - MATLAB meta.class objects provide a lot of useful material to document
%   classes and their members. E.g.: 
%   ```
%   Mc = meta.class.fromName( 'Documentor' ) ; 
%   ```
%   Mc.MethodList will then contain info on Documentor's methods -- including
%   the input names (counting them then gives the number of arguments).
%   Curiously however, it seems that if you add an arguments block to a class
%   method---even if you don't specify default values, which normally renders
%   arguments optional---suddenly the Input names disappear and are replaced by
%   an uninformative 'varargin'?? 
%
%   It generally seems that the whole MATLAB meta.class library isn't
%   fully/coherently implemented. (e.g. many properties that sound as though
%   they would be very useful (e.g. Descrition, DetailedDescription) are not
%   currently used (according to their own documentation
%   <https://www.mathworks.com/help/matlab/ref/meta.class.html meta.class>
%   Further, if you convert call 
%   ```
%   struct( Mc.MethodsList )
%   ```
%   hidden properties are displayed (e.g. InputTypes) that clearly haven't been
%   implemented (e.g InputTypes always seems to be {'any'} even when an
%   arguments block clearly states what is permissible)
%
% 2. Make (or identify existing) templates for .m files to facilitate
% auto-documentation (i.e. standard formating for help section of scripts,
% functions, classes + class members).  E.g. Should there be typically sections
% like ### Summary ###, ### Inputs ###, ### References ###, etc.? 
%
% 3. Configure the output format for documentation (typesetting, indentation, etc.)
%   e.g. what should generally be emboldened?  italicized?
%
% 4. Function to write or modify a mkdocs.yml configuration file according to
% new documentation file/directory structure
%
% 5. Source control (git) + host (e.g. readthedocs) integration:
%   when/when not to rebuild doc files and push?
%
% 6. Adding+formatting <https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet#links links>
%
%   - Handing differences between Matlab markup and Markdown syntaxes?
%   
%   - Adding hyperlinks in the documentation to the corresponding code on github?
%   Is this something readthedocs can do automatically? If the file structure
%   is kept the same (i.e. doc directory tree mirrors that of the code) this
%   should be straightforward?
% 
%   Decide among the approaches listed here which one would be easiest to use : 
%   <https://www.mkdocs.org/user-guide/writing-your-docs/#internal-links>
%
% 7. Git organization: e.g. make helpMeDoc a package in a separate repo in NeuroPoly and
%    include as a 'submodule' in realtime_shimming? It's not specific to the
%    shim project and---while the world waits for MathWorks to provide this sort
%    of functionality---it could probably be useful to others
%
% 8. Unit testing?
% 
% 9. Incorporating options generally into the documentation process (e.g. as a
% parameters struct input to Documentor(), a configurable text file, or
% something else).  As it stands, Documentor will include pretty much all the
% available info; however, having the option of simplified 'user' documentation
% would likely be useful.
%
%  - Adding option (perhaps as a logical called "isGeneratingOverview") to
%  generate an overview/table of contents file, consisting of the code names
%  and the summary header line (i.e. from Informer.Attributes.Description).
%  This might be similar to MATLAB's Contents.m files, or the displayed return
%  from LOOKFOR (e.g. type: lookfor image) e.g.
%  <https://octave.sourceforge.io/octave/overview.html>
%
%  - Adding the option of including implementation details in the
%  documentation?  e.g. what internal calls to custom (or matlab) functions are
%  made in a given piece of code?  (e.g. see
%  <https://www.mathworks.com/help/matlab/ref/checkcode.html codecheck>)

    
 

