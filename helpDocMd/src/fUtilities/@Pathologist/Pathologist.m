classdef Pathologist 
% Pathologist File path utility class
% 
% ## Purpose
%
% Pathologist contains methods for handling paths to files and folders. 
%
% Since a path designation in MATLAB can take one of three forms (char-,
% string-, or cellstr-array), a primary purpose the class is merely to provide
% a convenient means of handling the different implementations. Namely, a set
% of input paths will always, first, be typecast as a string and stored as the
% property `data`. After processing via any `Pathologist` methods, paths are
% returned in the form of their input (see Pathologist.returnasinput).
%
% **NOTE** The name `pathtool` already being taken by a builtin MATLAB
% function, in keeping with the `Doc. Md.` theme and the unfortunately current
% zeitgeist due to COVID-19, "Pathologist" seems an appropriate moniker.

properties( Access=private )

    % Input data-type (aka 'class') 
    typeIn(1,1) string {mustBeMember(typeIn, ["string" "char" "cell"])} = "char" ;

    sizeIn = 1

    % Path data: real or fictitious files or folders
    data {mustBeString} = "./" ;

end

properties( Dependent )

    % Number of input paths
    nPaths ;

    % Path(s) to immediate parent folder(s) (single level)
    % Main doc. entry (for now?): parent.m
    parentDir ;

    % Deepest common parent folder (the "fork")
    % `if Path.nPaths == 1, Path.baseDir = Path.parentDir`
    % Main doc. entry (for now?): parent.m
    baseDir(1,:) ;

end

properties( Access=private, Dependent)

    % isExisting(1,1)
    % areAllExisting(1,1)
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Path = Pathologist( pathIn )
    
    if nargin == 0
        return ;
    end

    Path.typeIn = string( class( pathIn ) ) ;
    Path.sizeIn = size( pathIn ) ;
    Path.data   = Pathologist.typecast( pathIn, "string" ) ;

end
% =========================================================================    
function [baseDir] = get.baseDir( Path )

    baseDir = retracepath( Path.data(:) ) ;
    baseDir = Path.returnasinput( baseDir ) ;

    %% Retrace path 
    % NOTE this doesn't need to be a function if it stays here, but might be
    % more sensible as a static function...) 
    function [baseDir] = retracepath( pathIn )
    % Reassemble the shared path one sub dir at at time, starting at the root;
    % return where it diverges. 
        baseDir     = unique( extractBefore( Path.data, 2 ) ) ; % presumably == filesep
        nParentsMin = min( count( Path.data, filesep ) ) ;

        if numel( baseDir ) ~= 1 % nothing in common :(
            baseDir = "?" ; 
            return ;
        end

        iParent       = 1 ;
        pathRemaining = Path.data ;

        while all( strlength( pathRemaining)>0 ) && ( iParent < nParentsMin )

            iParent = iParent + 1 ;

            [subDir, pathRemaining] = strtok( pathRemaining, filesep ) ;
            uniqueSubDir            = unique( subDir ) ;

            if numel( uniqueSubDir ) == 1 
                baseDir = fullfile( baseDir, uniqueSubDir ) ;
            
            else % path forked at iParent-1; exit while
                pathRemaining = "" ;
            end
        end

    end % retracepath()

end
% =========================================================================    
function [nPaths] = get.nPaths( Path )
    
    nPaths = numel( Path.data ) ;
end
% =========================================================================    
function [parentDir] = get.parentDir( Path )
    
    %% Get parent directories (single level)
    if Path.nPaths == 1
        parentDir = fileparts( Path.data ) ;
        
        % 2x if input was of the form: '.../folder/'
        if endsWith( Path.data, filesep ) 
            parentDir = fileparts( parentDir ) ;
        end
    else
        parentDir = arrayfun( @fileparts, Path.data ) ;
    end

    parentDir = returnasinput( Path, parentDir ) ;

end
% =========================================================================    

end % methods
% =========================================================================    
% =========================================================================    
methods( Access=private )
% =========================================================================    
function [varargout] = returnasinput( Path, varargin ) 
% RETURNASINPUT Resize and recast a path value to correspond with user input
%    
%     [pathOut] = returnasinput( Path )
%     [pathOut] = returnasinput( Path, pathIn )
%     [pathOut1, pathOut2, ...] = returnasinput( Path, pathIn1, pathIn2,... )
%
% Reshapes `pathIn` according to `Path.sizeIn` when the 2 have the same number
% of elements and recasts it as `Path.typeIn` to return `pathOut`.
%
% When called with the single argument, `pathIn` is assigned the value of
% `Path.data`. 
%
% When called with > 2 inputs, resizing + recasting is applied successively to
% define the respective returns.

    if nargin == 1
        pathIn = Path.data ;
    elseif nargin == 2
        pathIn = varargin{1} ;
        if numel( pathIn ) > 1
            if strcmp( Path.typeIn, 'char' )
                pathIn = reshape( pathIn, [Path.sizeIn(1) 1 Path.sizeIn(3:end)] ) ;
            else
                pathIn = reshape( pathIn, Path.sizeIn ) ; 
            end
        end

        varargout{1} = Path.typecast( pathIn, Path.typeIn ) ;

    elseif nargin > 2
        for iP = 1 : numel( varargin )
            varargout{iP} = Path.returnasinput( varargin{iP} ) ;
        end
    end

end
% =========================================================================    

end % private methods
% =========================================================================    
% =========================================================================    
methods( Static )
    %..... 
    [pathOut, pathType] = abs( pathIn )
    %..... 
    [txtOut] = typecast( txtIn, castAs )
end
% =========================================================================    
% =========================================================================    

% =========================================================================    
%
end % classdef Pathologist
