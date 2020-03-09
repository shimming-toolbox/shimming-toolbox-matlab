function [Att] = metainfo( MetaObj )
%EXTRACTMETAINFO Return info (struct) derived from a meta.object
% 
% ### Syntax
%    
%    Att = metainfo( MetaObj )
%
% Useful properties of meta[.class/.property/.method] object 'MetaObj' are
% copied to the fields of struct Att
%
%TODO : restructuring of meta.Validation.Size

SUPPORTED_TYPES = ["meta.class" ; "meta.method" ; "meta.property"] ;

narginchk(1,1) ;

metaObjType = unique( string( class( MetaObj ) ) ) ;

assert( ( numel( metaObjType ) == 1 ) && ismember( metaObjType, SUPPORTED_TYPES ), ...
    ['Array elements of input MetaObj should all be of the same type\n' ...
    'i.e.: ' char( join( SUPPORTED_TYPES, " or " ) ) ], '%s' ) ; 

names = string( properties( metaObjType ) ) ;

%% Initialize empty struct with field names = property names
for iName = 1 : numel( names ) ;
    Att.( names( iName ) ) = [] ;
end

%% Copy property values

for iObj = 1 : numel( MetaObj )
    for iName = 1 : length( fieldnames( MetaObj( iObj ) ) )

        name = names( iName ) ;

        try
            Prop = MetaObj( iObj ).( name ) ;
        catch Me
        % 'try' throws an error if Prop = meta.property.DefaultValue is called
        % when no DefaultValue has been set. (there may be other cases?) 
            Prop = "" ;
        end

        propClass = class( Prop ) ;

        switch propClass 
            case 'logical' 
            % applies to several attributes for each of the SUPPORTED_TYPES
                Att(iObj).( name ) = Prop ;

            case { 'char' ; 'string' } 
                Att(iObj).( name ) = string( Prop ) ;

            case 'function_handle'
            % *potentially applies to Access, GetAccess, SetAccess properties
            
                if isempty( Prop )
                    Att(iObj).( name ) = "" ;
                else
                    s = string( func2str( Prop ) ) ;
                    if contains(s, ".m>")
                        s = extractAfter( s, ".m>" ) ; % remove local path
                    end
                    Att(iObj).( name ) = s ;
                end

            case { 'meta.method' ; 'meta.property' } 
                Att(iObj).( name ) = Informer.metainfo( Prop ) ;
            
            case 'meta.class'
            % copy the .Name entry only 
                Att(iObj).( name ) = string( { Prop.Name } ) ;
            
            case 'meta.Validation'
                % Att(iObj).( name ) = getvalidationinfo()
                if isempty( Prop )
                    Att(iObj).( name ) = [] ;
                else
                    V = struct( 'Class', "", 'Size', [], 'ValidatorFunctions', "" ) ;
                    if ~isempty( Prop.Class )
                        V.Class = string( Prop.Class.Name ) ;
                    end

                    if ~isempty( Prop.Size )
                        V.Size = Prop.Size ; %TODO : restructuring of meta.Validation.Size
                    end

                    if ~isempty(Prop.ValidatorFunctions)
                        for iF = 1 : numel( Prop.ValidatorFunctions )
                            V.ValidatorFunctions(iF) = string( func2str( Prop.ValidatorFunctions{iF} ) ) ;
                        end
                    end
                    Att(iObj).( name ) = V ;
                end
            
            case 'cell'
                if iscellstr( Prop ) 
                    Att(iObj).( name ) = string( Prop ) ;
           
                elseif isempty( Prop )
                    Att(iObj).( name ) = "" ;
           
                else
                    s = strings( [ 1 numel( Prop ) ] ) ;
                    for iP = 1 :numel( Prop )
                % applies to: meta.class.InferiorClasses
                % and (potentially) applies to:
                % - meta.method.Access
                % - meta.property.GetAccess, .SetAccess
                % --> copy class names

                        if isa( Prop{iP}, 'meta.class' )
                            s(iP) = string( Prop{iP}.Name ) ;
                        else
                            try 
                                sp(iP) = string( Prop{iP}.Name ) ;
                            catch Me
                                Me.rethrow ;
                            end
                        end
                    
                        Att(iObj).( name ) = s ;
                    end
                end
            end
    end
   
    %% Fill the otherwise empty 'Description' fields 
    switch char(metaObjType)
        case 'meta.class'
            mHelp = Informer.gethelptext( Att(iObj).Name ) ;
        otherwise
            mHelp = Informer.gethelptext( [ strcat( Att(iObj).DefiningClass, ".", Att(iObj).Name ) ] ) ;
    end
    Att(iObj).Description         = Informer.extracthelpheader( mHelp, Att(iObj).Name ) ;
    Att(iObj).DetailedDescription = Informer.extracthelpbody( mHelp ) ;

end
