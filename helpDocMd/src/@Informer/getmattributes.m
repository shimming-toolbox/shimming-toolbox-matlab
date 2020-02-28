function [Att] = getmattributes( mFile )
%GETMATTRIBUTES Return functional description of a .m file
% 
% ### Syntax ###
%
% Att = GETMATTRIBUTES( mFile )
% 
% ### Inputs ###
%
% - mFile: a path string to a single .m file 
%
% ### Outputs ###
%
% - Att: a struct of .m file attributes containing the following basic fields:
% 
%   - mType: the type of .m file
%
%   - .Name: Name of the script, function, class or class method 
%       
%   - .Description: Header line of help-text (string vector returned from
%   Informer.extracthelpheader) 
%
%   - .DetailedDescription: Body of help-text (string vector returned from
%   Informer.extracthelpbody) 
%
% ### Function and method files ###
%
% If the .m file is a function or class method, Att additionally
% contains fields: 
%
% - .nInputs: Number of input arguments 
%
% - .nOutputs: Number of output arguments
%
%   (TODO: implement means of adding (at least for inputs): argument names,
%   types, validation functions, and defaults, when declared in an
%   arguments block, and add these as fields to Attributes)
%
% ### Classdef files ###
%
% If the .m file is a class definition, Attributes copies the following
% additional fields from an instance of the associated meta.class object
% (itself returned as the second output argument 'Mc'):
%
%   [Att,Mc] = getmattributes( mFile )
% 
% Logicals describing class attributes:
% .Hidden 
% .Sealed
% .Abstract
% .Enumeration
% .ConstructOnLoad
% .HandleCompatible
% .RestrictsSubclassing
%
% .SuperclassList: String array of superclass names
% .InferiorClasses: String array of inferior class names
% .ContainingPackage: Name of the containing package as a string-scalar.    
%
% For more information on meta.class objects, refer to the MATLAB documentation 
% <https://www.mathworks.com/help/matlab/ref/meta.class.html>
%
% #### Note #### 
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
% ### References ###
%
% See also 
%
% - Informer.mfiletype
% - Informer.gethelptest
% - Informer.gethelpbody
% - Informer.gethelpheader
% - <https://www.mathworks.com/help/matlab/ref/meta.class.html meta.class>
%
%    
% % Methods struct derived from meta.method object (applies to mType="classdef" only, otherwise empty)
    % Methods = [] ;
    %
    % % mType="classdef" only: Properties struct derived from meta.property object. (applies to mType="classdef" only, otherwise empty)
    % Properties = [] ;
    % % unimplemented/TODO 
    % Enumerations struct ;
    % % unimplemented/TODO 
    % Events struct ;
    % %    
    % See also:
    % <https://www.mathworks.com/help/matlab/ref/meta.property.html>
    % <https://www.mathworks.com/help/matlab/ref/meta.validation-class.html meta.Validation>
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
                Att = adddescription( Att ) ;
            case { "function" }
            % TODO: elaborate function Attributes (e.g. info derived from
            % parsing arguments block) 
                Att = addfunctionattributes( Att ) ;
            
            % TODO: refactor the classdef/method cases, along with
            % Informer.metainfo() (which they both call and which is itself
            % pretty ugly...)
            case { "classdef" }
                Mc  = meta.class.fromName( Att.Name ) ;
                Att = Informer.metainfo( Mc ) ;
                Att.mType = "classdef" ;

            case { "method" } 
                % Get a meta.method object pertaining to the method:
                % this requires instantiating a _meta.class_ object, and
                % then picking out the specific meta.method object of interest: 
                [~, classFolderName] = fileparts( mFolder ) ;
                assert( startsWith( classFolderName, "@" ), 'Unexpected result. See code.' )
                
                className = erase( classFolderName, "@" ) ;
                Mc        = meta.class.fromName( className ) ;
                % Mc may have many meta.method objects. Extract the one
                % pertaining to the given .m method file & input it to
                % Informer.metainfo.
                iMethod   = [ string( {Mc.MethodList.Name} ) == BasicAtt.Name ] ;
                MM        = Mc.MethodList( iMethod ) ;

                Att       = Informer.metainfo( MM ) ;
                Att.mType = "method" ; 
 
             otherwise
        
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
function [Att] = adddescription( Att )
    
    mHelp                   = Informer.gethelptext( Att.Name ) ;
    Att.Description         = Informer.extracthelpheader( mHelp, Att.Name ) ;
    Att.DetailedDescription = Informer.extracthelpbody( mHelp ) ;

end

function [Att] = addfunctionattributes( Att )

    Att          = adddescription( Att ) ; 
    Att.nInputs  = nargin( Att.Name ) ;
    Att.nOutputs = nargout( Att.Name ) ;

end

%NOT used: (replaced with Informer.metainfo() )
function [Att] = getmethodattributes( BasicAtt )
%GETMETHODATTRIBUTES Return struct of method attributes from meta.class object
%
% ### Syntax ###
%
%  [Att] = GETMETHODATTRIBUTES( BasicAtt )
% 
% Copies elements of the following meta.class object properties to the fields
% of attributes struct 'Att':
%
% Access
% Static
% Abstract
% Sealed
% ExplicitConversion
% Hidden
% InputNames
% OutputNames
% DefiningClass (class names as strings) 
% 
% ### Notes ### 
%
% 1. Currently, all the hidden properties (and some of the non-hidden
% ones) of meta.method objects appear to be quite useless--despite their
% promising sounding names! Therefore, these are not copied to Att. 
%
% 2. For now, to (hopefully) avoid any strange bugs/behaviour from the
% meta.class package, the list of properties to copy from the meta.method
% instance is explictly/hard-coded rather than copying all the properties
% wholesale.  
 
    Att = BasicAtt ;

    Mc  = meta.class.fromName( Att.Name ) ;
    McM = Mc.MethodList( [string({Mc.MethodList.Name}) == string( Att.Name )] ) ;

    % for McM.DefiningClass, just keep the name:
    Att.DefiningClass = string( McM.DefiningClass.Name ) ;

    fields = [ "Access" ;
               "Static" ;
               "Abstract" ;
               "Sealed" ;
               "ExplicitConversion" ;
               "Hidden" ;
               "InputNames" ;
               "OutputNames" ; ] ;
    
    for iField = 1 : numel( fields )
        Att.( fields(iField) ) = McM( iField ) ;
    end
    
    Att.nInputs  = numel( Att.InputNames ) ;    
    Att.nOutputs = numel( Att.OutputNames ) ;    

end

function [Att] = getclassattributes( BasicAtt )
%GETCLASSATTRIBUTES Return struct of class attributes from meta.class object
% 
% ### Syntax ###
%
%  [Att] = GETCLASSATTRIBUTES( BasicAtt )
% 
% Copies elements of the following meta.class object properties to the fields
% of attributes struct 'Att':
%
% Hidden
% Sealed
% Abstract
% Enumeration
% ConstructOnLoad
% HandleCompatible
% InferiorClasses
% ContainingPackage
% RestrictsSubclassing
% PropertyList
% MethodList
% EventList
% EnumerationMemberList
% SuperclassList
    Att = BasicAtt ;
    Mc  = meta.class.fromName( Att.Name ) ;

    fields = string( fieldnames( Mc ) ) ;

    %% Basic attributes:
    for iField = 1 : length( fields )
        if islogical( Mc.( fields(iField) ) ) 
            % Applies to: "Hidden", "Sealed", "Abstract", "Enumeration",
            % "ConstructOnLoad", "HandleCompatible", "RestrictsSubclassing"
            Att.(fields{iField}) = Mc.( fields(iField) ) ;
        end
    end

    %% Inheritance and package info:

    if isempty( Mc.SuperclassList ) 
        Att.SuperclassList = "N/A" ;
    else
        Att.SuperclassList = string( {Mc.SuperclassList.Name} ) ;
    end

    if isempty( Mc.InferiorClasses ) 
        Att.InferiorClasses = "N/A" ;
    else 
        % NOTE: 'Mc.InferiorClasses' is a cell array whereas Mc.SuperclassList
        % is an object array, hence the for-loop: 
        names = string( Mc.InferiorClasses{1}.Name ) ;

        for iClass = 2 : numel(Mc.InferiorClasses)
            names = [ names string( Mc.InferiorClasses{ iClass }.Name ) ] ;
        end

        Att.InferiorClasses = names ;
    end

    if isempty( Mc.ContainingPackage ) 
        Att.ContainingPackage = "N/A" ;
    else
        Att.ContainingPackage = string( Mc.ContainingPackage.Name ) ;
    end
    
end %getclassattributes()
