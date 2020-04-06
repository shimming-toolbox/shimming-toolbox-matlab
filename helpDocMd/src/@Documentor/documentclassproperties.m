function docStr = documentclassproperties( Info, isDetailed )
%DOCUMENTCLASSPROPERTIES Return string vector of class property documentation
    arguments
        Info struct ;
        isDetailed {mustBeBoolean} = true ;
    end

% assert( strcmp(Info.mType, "classdef"), 'mFile is not a class' ) ;

docStr = [ "- - -" ; "## Properties" ; "" ] ;

if isempty( Info.PropertyList )
    docStr = [ docStr ; "" ; "[No Properties]"; "" ] ;
    return ;
else
    for iProp = 1 : numel( Info.PropertyList ) 

        Prop = Info.PropertyList( iProp ) ;
        
        if isDetailed || ( strcmp( Prop.GetAccess, 'public' ) && ~Prop.Hidden )
            
            docStr = [ docStr ; "" ; Documentor.documentbasic( Prop, 3 ) ] ;

            % remove fields addressed already in documentbasic
            Prop = rmfield( Prop, {'Name'; 'Description' ; 'DetailedDescription'} ) ;

            [attTable, Prop] = tablelogicalattributes( Prop ) ;
            docStr     = [ docStr ; attTable ; "" ] ;
            
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

function [ attTable, Prop] = tablelogicalattributes( Prop )
%TABLELOGICALATTRIBUTES Document logical attributes as an HTML table
%
% document logical attributes as an HTML table (moved here to trim down ugly for loop above)

    propFields = fieldnames( Prop ) ;
    boolFields = propFields( structfun( @islogical, Prop ) ) ;

    for iF = 1 : numel( boolFields )
       BoolAttributes.( boolFields{iF} ) = Prop.( boolFields{iF} ) ;
    end

    attTable = Documentor.tableattributes( BoolAttributes ) ;
    Prop     = rmfield( Prop, boolFields ) ;

end 

end
