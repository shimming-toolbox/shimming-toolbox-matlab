classdef ShimComAcdc % < ShimCom

    
    
    
    
    
    
    
    
    
methods
    
    function Shims = ShimComAcdc( Params )
        
        if nargin < 1 || isempty( Params )
            Params.dummy = [] ;
        end
        
        Shims.Specs = ShimSpecs( Params ) ;
        
        Shims.Cmd = ShimComAcdc.definecommands() ;
        
        Shims.Comport = ShimComAcdc.initializecomport( Shims.Specs) ;
        
        Shims.Data.Output = '' ;
        Shims.Data.Input = '' ;
        
        Shims.Params.nBytesToRead = [] ;
        Shims.Params.nAttempsBeforeError = 3 ;
        
    end
    % =====================================================================
    function Comport = initializeport( Shims )
    % INITIALIZEPORT
    % 
    %
    % Create the serial Object to communicate with the Arduino with the
    % selected port in Params
    
    
        DEFAULT_PORT = 'COM1' ;
        
        if ~myisfield( Params, 'PortName' ) || isempty(Params.PortName)
            Params.PortName = DEFAULT_PORT ;
        end
        
        Comport = serial(Params.PortName) ;
        % fopen(Comport) ; ?
    end
    % =====================================================================
    function [] = setandloadshim( Shims , iChannel , Current )
    % SETANDLOADALLSHIMS    
    %
    %
    % Set the selected channel to the selected current and print the
    % arduino's answer in the dialog box
        if iChannel > 0 && iChannel < Shims.Specs.Amp.nActiveChannels
            fopen( Shims.Comport ) ;
            
            fprintf( Shims.Comport, Shims.Cmd.SetAndLoadChannel{iChannel} + char(Current)) ;
            answer = fscanf(Shims.Comport, '%f', 1) ;
            ShimUse.display(answer) ;
        else
            error('Channel is not available, feel free to mount more on the coil ;) ') ;
        end
    end
    % =====================================================================
    function [] = setansloadallshimstoonevalue( Shims, Current )
        
        fopen( Shims.Comport ) ;
        
        fprintf( Shims.Comport, Shims.Cmd.setAndLoadAllShimsToOneValue + char(Current) ) ; 
        answer = fscanf( Shims.Comport, '%f', 1 ) ;
        
        ShimUse.display(answer) ;
        
        fclose( Shims.Comport ) ;
    end
    % =====================================================================
    function [] = setandloadallshim( Shims, Current )
        
        nChannels = Shims.Specs.Amp.nActiveChannels ;
        
        if length(Current) ~= nChannels
            error('Number of currents in input must be ' + num2str(nChannels)) ;
        end
        
        fopen( Shims.Comport ) ;
        
        for iChannel = 1:nChannels
            fprintf( Shims.Comport, Shims.Cmd.setAndLoadChannel{iChannel} + char( Current(iChannel) ) ) ;
            answer = fscanf( Shims.Comport, '%f', 1 ) ;
            
            ShimUse.display(answer) ;
        end
        
        fclose( Shims.Comport ) ;
    end
    
    % =====================================================================
    function [] = setallshimtozero ( Shims )
        
        fopen( Shims.Comport ) ;
        
        fprintf( Shims.Comport, Shims.Cmd.returnToZero ) ; 
        answer = fscanf( Shims.Comport, '%f', 1 ) ;
        
        ShimUse.display(answer) ;
        
        fclose( Shims.Comport ) ;
    end
    % =====================================================================
    function [] = loadallshims( Shims )
        fopen( Shims.Comport ) ;
        
        fprintf( Shims.Comport, Shims.Cmd.loadAllShims ) ; 
        answer = fscanf( Shims.Comport, '%f', 1 ) ;
        
        ShimUse.display(answer) ;
        
        fclose( Shims.Comport ) ;
    end
    % =====================================================================
    function [] = calibrateallshims( Shims )
        fopen( Shims.Comport ) ;
        
        fprintf( Shims.Comport, Shims.Cmd.calibrateAllShims ) ; 
        answer = fscanf( Shims.Comport, '%f', 1 ) ;
        
        ShimUse.display(answer) ;
        
        fclose( Shims.Comport ) ;
    end
    % =====================================================================
    
    function Cmd = definecommands( Shims )
        
        Cmd.setAndLoadAllShimsToOneValue = 'A' ;
        
        Cmd.calibrateAllShims = 'C' ;
        
        Cmd.returnToZero = 'Z' ;
        
        Cmd.LoadAllShims = 'V' ;
        
        Cmd.setAndLoadChannel{1} = 'a';
        Cmd.setAndLoadChannel{2} = 'b';
        Cmd.setAndLoadChannel{3} = 'c';
        Cmd.setAndLoadChannel{4} = 'd';
        Cmd.setAndLoadChannel{5} = 'e';
        Cmd.setAndLoadChannel{6} = 'f';
        Cmd.setAndLoadChannel{7} = 'g';
        Cmd.setAndLoadChannel{8} = 'h';
        
    end
end
end