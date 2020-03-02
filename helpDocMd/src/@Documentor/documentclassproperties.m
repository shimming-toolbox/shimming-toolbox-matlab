function docStr = documentclassproperties( Dr )
%DOCUMENTCLASSPROPERTIES Return string vector of class property documentation

Info = Dr.Info.Attributes ;

assert( strcmp(Info.mType, "classdef"), 'mFile is not a class' ) ;

docStr = [ "- - -" ; "### Properties ###" ; "" ] ;

if isempty( Info.PropertyList )
    docStr = [ docStr ; "" ; "[No Properties]"; "" ] ;
    return ;
else
    for iProp = 1 : numel( Info.PropertyList ) 
        docStr = [ docStr ; "" ] ;

        Prop = Info.PropertyList( iProp ) ;
        
        if ~isempty( Prop.Description )
            docStr = [ docStr ; strcat( "##### ", Prop.Name, " #####" ) ; "" ; strcat( "_", Prop.Description, "_" ) ] ;
        else % name only
            docStr = [ docStr ; strcat( "##### ", Prop.Name, " #####" ) ; "" ] ;
        end
        
        if ~strcmp( Prop.DetailedDescription, Prop.Description )                         
        % NOTE: could change Informer so it doesn't fill DetailedDescription with the copy of Description
            docStr = [ docStr ; "" ; Prop.DetailedDescription ] ;
        end
            
        Prop = rmfield( Prop, {'Name'; 'Description' ; 'DetailedDescription'} ) ;
        
        if Dr.isDetailed

            propFields  = fieldnames( Prop ) ;

            %% Place basic (logical) attributes into an HTML table
            % i.e.:
            % basicFields = { 'Dependent', 'Constant', 'Abstract', 'Transient', 'Hidden' ; ...
            %             'GetObservable', 'SetObservable', 'AbortSet', 'NonCopyable', 'HasDefault' } ;
            basicFields = propFields( structfun( @islogical, Prop ) ) ;
            
            for iF = 1 : numel( basicFields )
               BasicAttributes.( basicFields{iF} ) = Prop.( basicFields{iF} ) ;
            end
            
            docStr = [ docStr ; "" ; Dr.tableattributes( BasicAttributes ) ; "" ] ;
            clear BasicAttributes ;
            Prop = rmfield( Prop, basicFields ) ;
            
            %% 
            propFields = fieldnames( Prop ) ;
         
            for iField = 1 : length( propFields ) 

                field = string( propFields( iField ) ) ;

                if isempty( Prop.(field) )
                    docStr(end+1) = strcat( "- ", field, " : [N/A] " ) ;

                elseif strcmp( field, "Validation" )

                    docStr(end+1) = "- Validation: " ;

                    if Prop.Validation.Class ~= "" 
                        docStr(end+1) = strcat( "Class: ", Prop.Validation.Class, "" ) ;
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
