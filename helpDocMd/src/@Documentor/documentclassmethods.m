function docStr = documentclassmethods( Dr )
%DOCUMENTCLASSMETHODS Return string vector of class method documentation

Info = Dr.Info.Attributes ;

assert( strcmp(Info.mType, "classdef"), 'mFile is not a class' ) ;

docStr = [ "---" ; "### Methods ###" ; "" ] ;

if isempty( Info.MethodList )
    docStr = [ docStr ; "" ; "[No Methods]"; "" ] ;
    return ;
else
    mthdFields = fieldnames( Info.MethodList ) ;
     
    for iMthd = 1 : numel( Info.MethodList ) 
        docStr = [ docStr ; "" ; "---" ; "" ;] ;

        Mthd = Info.MethodList( iMthd ) ; 
         
        if contains( which( Mthd.DefiningClass ), 'built-in' ) 
        % MATLAB built-in method: include name and defining class only:
            docStr(end+1) = strcat( "### ", Mthd.Name ) ;
            docStr(end+1) = [ "[ _built-in method derived from *" + Mthd.DefiningClass + "* class_ ]" ];
            docStr(end+1) = strcat( "For more info, see MATLAB documentation]" ) ;
        else
            if ~isempty( Mthd.Description )
                docStr = [ docStr ; "" ; strcat( "### ", Mthd.Name ) ; "" ] ; 
                docStr = [ docStr ; "" ; strcat( " _", Mthd.Description, "_ " ) ; "" ] ;
            else 
                docStr = [ docStr ; "" ; strcat( "### ", Mthd.Name ) ; "" ] ;
            end

            if ~isempty( Mthd.DetailedDescription )
                % NOTE: could change Informer so it doesn't fill DetailedDescription with the copy of Description
                if ~strcmp( Mthd.DetailedDescription, Mthd.Description )                         
                    docStr = [ docStr ; "Description: " ; Mthd.DetailedDescription ] ;
                end
            end
        
            % start at 4 to skip name, description, detailed description:
            for iField = 4 : length( mthdFields ) 
            
                field = string( mthdFields( iField ) ) ;

                if isempty( Mthd.(field) )
                    docStr(end+1) = strcat( "- ", field, " : [N/A] " ) ;
                else 
                    docStr(end+1) = strcat( "- ", field, " : ", strjoin( string( Mthd.( field ) ), ", " ) ) ;
                end
            end
        end
    end
end

end
