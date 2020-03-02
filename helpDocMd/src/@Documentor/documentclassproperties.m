function docStr = documentclassproperties( Dr )
%DOCUMENTCLASSPROPERTIES Return string vector of class property documentation
    
% { 'GetAccess', 'SetAccess', 'Dependent', 'Constant', 'Abstract', 'Transient', 'Hidden',
%   'GetObservable', 'SetObservable', 'AbortSet', 'NonCopyable', 'GetMethod', 'SetMethod', 'HasDefault' } 

Info = Dr.Info.Attributes ;

assert( strcmp(Info.mType, "classdef"), 'mFile is not a class' ) ;

docStr = [ "---" ; "### Properties ###" ; "" ] ;

if isempty( Info.PropertyList )
    docStr = [ docStr ; "" ; "[No Properties]"; "" ] ;
    return ;
else
    propFields = fieldnames( Info.PropertyList ) ;
     
    for iProp = 1 : numel( Info.PropertyList ) 
        docStr = [ docStr ; "" ] ;

        Prop = Info.PropertyList( iProp ) ; 
        
        if ~isempty( Prop.Description )
            nameAndDescription = strcat( "*", Prop.Name, "*", " : ", "_", Prop.Description, "_" )
            docStr(end+1) = nameAndDescription;
        else % name only
            docStr(end+1) = strcat( "*", Prop.Name, "*" );
        end

        if ~isempty( Prop.DetailedDescription )
            % NOTE: could change Informer so it doesn't fill DetailedDescription with the copy of Description
            if ~strcmp( Prop.DetailedDescription, Prop.Description )                         
                docStr = [ docStr ; "Description: " ; Prop.DetailedDescription ] ;
            end
        end
        
        for iField = 4 : length( propFields ) % start at 4 to skip name, description, detailed description:
            
            field = string( propFields( iField ) ) ;

            if isempty( Prop.(field) )
                docStr(end+1) = strcat( "- ", field, " : [N/A] " ) ;
            
            elseif strcmp( field, "Validation" )

                docStr(end+1) = "- Validation: " ;

                if Prop.Validation.Class ~= "" 
                    docStr(end+1) = strcat( "Class: ", Prop.Validation.Class ) ;
                end

                docStr(end+1) = strcat( "Validator functions: ", strjoin( Prop.Validation.ValidatorFunctions, "," ) ) ;
                
                %TODO add size constraints
            else 
                docStr(end+1) = strcat( "- ", field, " : ", string( Prop.( field ) ) ) ;
            end
        end 
    end
end

end
    % if strcmp( mfiletype( Dr.pathIn ), 'classdef' )
    %
    % mdDoc = Dr.mHelp ;
    %     %% Class members: Properties    
    %     mdDoc = [ mdDoc; "" ; "### Members ###" ; ""]; 
    %
    %     Props = Mc.PropertyList
    %
    %     mdDoc = [ mdDoc; "#### Properties ####" ] ;
    %
    %     if isempty( Props )
    %
    %         mdDoc = [ mdDoc; "" ; "_(None)_" ; "" ] ;
    %     else
    %
    %         for iProp = 1 : length( Props )
    %
    %             Prop   = Props( iProp )
    %             mdDoc  = [ mdDoc ; strcat( "##### ", Prop.Name, "##### " ) ; "" ] ;
    %             fields = fieldnames( Prop ) ;
    %
    %             for iField = 1 : length( fields )   
    %                 field = fields{ iField } ;
    %
    %                 switch field
    %                     case { 'GetAccess', 'SetAccess', 'Dependent', 'Constant', 'Abstract', 'Transient', 'Hidden',
    %                            'GetObservable', 'SetObservable', 'AbortSet', 'NonCopyable', 'GetMethod', 'SetMethod', 'HasDefault' } 
    %                         if isempty( Prop.(field) )
    %                             entry = "_(None)_" ;
    %                         else
    %                             entry = join( string( Prop.( field ) ) ) ;
    %                         end
    %
    %                         mdDoc = [ mdDoc ; strcat( "- ", fields{iField}, ": ", entry ) ] ;
    %
    %                     otherwise
    %                         % do nothing
    %                 end
    %             end
    %     
    %             % TODO : Printing defaults: what is suitable to print (e.g. if
    %             % default is rand(100,100,100) clearly it would be preferable
    %             % to print the function call itself rather than the value).
    %             % might need to parse the code text itself. 
    %            % if Prop.HasDefault
    %            %        mdDoc = [ mdDoc ; "- Default: " ] ;
    %            %        isPrintingDefault=false ;
    %            %     try
    %            %        defaultString = splitlines( string( Prop.Default ) ) ;
    %            %        mdDoc = [mdDoc ; defaultString ] ; 
    %            %      catch Me
    %            %          defaultString = "_(Not printable)_" ;
    %            %    end
    %            %
    %            %  else
    %            %      mdDoc = [ mdDoc ; "- Default: _(None)_" ] ;
    %            % end
    %
    %
    %             end
    %
    %         end 
    %     end
