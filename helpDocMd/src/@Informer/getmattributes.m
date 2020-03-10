function [Att] = getmattributes( mFile )
%GETMATTRIBUTES Return functional description of a .m file
% 
% ### Syntax
%    
%    Att = GETMATTRIBUTES( mFile )
% 
% ### Inputs
%
% - mFile: a path string to a single .m file (a script, function, classdef, or class method)
%
% ### Outputs
%
% - Att: a struct of .m file attributes containing the following basic fields:
% 
%   - .mType: the type of .m file
%
%   - .Name: Name of the script, function, class or class method 
%       
%   - .Description: Header line of help-text (string vector returned from
%   Informer.extracthelpheader) 
%
%   - .DetailedDescription: Body of help-text (string vector returned from
%   Informer.extracthelpbody) 
%
% ### Function files
%
% If the .m file is a function, Att additionally contains fields: 
%
% - .nInputs: Number of input arguments 
%
% - .nOutputs: Number of output arguments
%
% *TODO*: implement means of adding (at least for inputs): argument names,
% types, validation functions, and defaults, when declared in an
% arguments block, and add these as fields to Attributes)
%
% ### Classdef files
%
% If the .m file is a class definition, Attributes derives the following 
% additional fields from an instance of the associated meta.class object:
% (For more info, refer to the MATLAB documentation: <https://www.mathworks.com/help/matlab/ref/meta.class.html>)
% 
% _Fields containing logical scalars_:
% 
% - .Hidden 
% - .Sealed
% - .Abstract
% - .Enumeration
% - .ConstructOnLoad
% - .HandleCompatible
% - .RestrictsSubclassing
%
% _Fields containing string arrays_:
%
% - .SuperclassList (the names of superclasses)
% - .InferiorClasses (the names of deriving inferior classes)
% - .ContainingPackage (the name of the containing package, if applicable, as a string-scalar)
%
% _Fields containing struct arrays_:
%
% - .MethodList (class methods derived from meta.method objects)
%
%   Elements of MethodList possess the following fields, with all but the final 4
%   (which contain string arrays) containing scalar logicals:
%   - .Static 
%   - .Abstract 
%   - .Static 
%   - .ExplicitConversion
%   - .Sealed
%   - .Hidden
%   - .Abstract
%   - .Access 
%   - .InputNames
%   - .OutputNames
%   - .DefiningClass 
%
% - .PropertyList (Class properties derived from meta.property objects)
%
% TODO: elaborate...
%
%
% TODO: add events and enumerations substructs from meta.class:
% - .EventList
% - .EnumerationMemberList
%
% #### Notes 
%
% Despite _sounding_ very useful based on the meta.class property names
% (including those normally hidden but made visible upon converting a meta.class
% object to a struct!), unfortunately, it seems that MathWorks hasn't yet
% bothered to implement anything for many of the meta.class properties that
% would actually make them meaningful (at least, not as of MATLAB 2019b). 
% 
% (Hence, for now, although the .Description and .DetailedDescription fields of
% the returned attributes struct borrow their field names from (nominal)
% meta.class properties, these are fields are actually assigned independently
% of meta.class by calls to Informer methods: gethelptext(), extracthelpheader()
% and extracthelpbody().)
%
% Furthermore, the current implementation of meta.class exhibits some curious behaviour:
% e.g. it seems that if you add an arguments block to a class method---even if
% you do not specify default values, which would normally render the arguments
% optional---suddenly the InputNames of the corresponding meta.method object
% disappear and are replaced by an uninformative 'varargin'??
%
% ### References
%
% See also 
%
% - Informer.mfiletype
% - Informer.gethelptest
% - Informer.gethelpbody
% - Informer.gethelpheader
% - <https://www.mathworks.com/help/matlab/ref/meta.class.html meta.class>
% - <https://www.mathworks.com/help/matlab/ref/meta.property.html>
% - <https://www.mathworks.com/help/matlab/ref/meta.validation-class.html meta.Validation>
    arguments
        mFile {mustBeStringOrChar, mustBeFile} = which("getmattributes.m") ;  
    end

[mFolder, Att.Name]        = fileparts( mFile ) ;
[Att.mType, mFile, mExist] = Informer.mfiletype( mFile ) ;

if strcmp( Att.mType, "NA" ) 
    warning('Not a valid .m file.') ;
    return ;
end

userDir = pwd ;
cd( mFolder ) ; 

try % Return to userDir if an error occurs

    switch Att.mType

        case { "script" }
            Att = addbasicdescription( Att ) ;

        case { "function" }
            Att = addfunctionattributes( Att ) ;
        
        case { "method" } 
            Att = addmethodattributes( Att ) ;
        
        case { "classdef" }
            Att = addclassattributes( Att ) ;

         otherwise
             error('Unexpected result. See code') ;
    
    end

    cd( userDir ) ;

catch Me

    cd( userDir ) ;
    Me.rethrow

end

end %getmattributes()

% -----
%% Local functions
% -----
function [Att] = addbasicdescription( Att )

    mHelp                   = Informer.gethelptext( Att.Name ) ;
    Att.Description         = Informer.extracthelpheader( mHelp, Att.Name ) ;
    Att.DetailedDescription = Informer.extracthelpbody( mHelp ) ;

end

function [Att] = addfunctionattributes( Att )
% TODO: elaborate function Attributes (e.g. with info derived from
% parsing arguments block -- TODO applies to class methods as well)

    Att          = addbasicdescription( Att ) ; 
    Att.nInputs  = nargin( Att.Name ) ;
    Att.nOutputs = nargout( Att.Name ) ;

end

function [Att] = addclassattributes( Att )
    
    Mc        = meta.class.fromName( Att.Name ) ;
    Att       = Informer.metainfo( Mc ) ;
    Att.mType = "classdef" ; % add .mType again since Att is overwritten in the above call

end

function [Att] = addmethodattributes( Att )
%ADDMETHODATTRIBUTES Return struct of method information
%
% Retrieves class method info from an associated meta.method object. (This
% requires first instantiating the meta.class object of the defining class, and
% then picking out the specific meta.method object of interest)
%
% TODO: elaborate attributes (e.g. with info derived from parsing arguments
% block)

    %% Create meta.class object pertaining to the defining class:

    % NOTE: the switch to the defining class folder should have occured above in main function 
    [~, classFolderName] = fileparts( string( pwd ) ) ;
    assert( startsWith( classFolderName, "@" ), 'Unexpected result. See code.' )
    
    className = erase( classFolderName, "@" ) ;
    Mc        = meta.class.fromName( className ) ;

    %% Get the meta.method object pertaining to the .m method file:
    iMethod   = [ string( {Mc.MethodList.Name} ) == Att.Name ] ;
    MM        = Mc.MethodList( iMethod ) ;
    
    %% Call Informer.metainfo( ) to get attributes struct:
    Att       = Informer.metainfo( MM ) ;
    Att.mType = "method" ; % add .mType again since Att is overwritten in the above call

end
