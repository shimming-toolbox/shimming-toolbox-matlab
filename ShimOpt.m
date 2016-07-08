classdef ShimOpt < MaRdI
%SHIMOPT
%
% Shim Optimization, a MaRdI subclass (i.e. ShimOpt < MaRdI)
%
% .......
%
% Shim = ShimOpt( pathToCalibrationInfo )
%
%   Shim contains fields
%
%       .img
%           Shim reference maps
%
%       .Hdr
%           Info Re: calibration data
%           (e.g. Hdr.MaskingImage defines the spatial support of the ref maps)
%
%       .Field
%           Object of type MaRdI pertaining to field distribution to be shimmed
%
%       .Model
%           .currents  
%               Optimal shim current vector (i)
%           .field     
%               Optimal shim field from projection of i onto reference maps (Ai)
%
%       .Parameters
%
%       .Probe
%           Object of type ProbeTracking
%           (i.e. the respiratory probe for real-time shimming.)
%
% =========================================================================
% Part of series of classes pertaining to shimming:
%
%     ProbeTracking
%     ShimCom
%     ShimOpt
%     ShimSpecs
%     ShimUse
%     
% =========================================================================
% Updated::ryan.topfer@polymtl.ca::Mon 20 Jun 2016 18:15:05 EDT
% =========================================================================

% =========================================================================
% *** TODO 
% .....
%   Shim optimization and (more importantly) field map interpolation
%   would be faster were the field maps initially cropped to the approximate
%   dimensions of the shim VOI. Procedurally (re:scan) this may not make sense
%   if our protocol remains [1) acq. field maps 2) acq. anatomical for VOI]...
%
% .....
% OPTIMIZESHIMCURRENTS()
%   Run (unconstrained) CG solver;
%   If solution achievable given system constraints, return solution;
%   Else, run fmincon given constraints & return that solution instead;
%
% .....
% EXTENDHARMONICFIELD()
%   Write function to check if field map is exists + is reasonable over
%   the shim voi. if some portion is missing, a harmonic extension could
%   be performed to fill in some of the missing field.
%   
% .....
% extractharmonicfield() & optimizeshimcurrents()
%   
%   documentation + optional params
% .....
% RESLICEIMG()
%   griddata takes too f-ing long.
%   write interp function in cpp
% .....
% RE: Pressure reading (e.g. CALIBRATE...) 
%   Do not write to file until finish. (unless this can be done very rapidly!) 
%
% =========================================================================

properties
    Field ;
    Model ;
    Parameters ;
    Probe ;
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimOpt( pathToCalibrationInfo )
%SHIMOPT - Shim Optimization
%
    fprintf(['\n Preparing for shim \n\n #####  \n\n'...
            'Loading SpineShim calibration info from ' pathToCalibrationInfo '\n\n']) ;

    load( pathToCalibrationInfo ) ;

    Shim.Hdr  = SpineShim.Hdr ;

    nVoxels = numel( SpineShim.img(:,:,:, 1) ) ;

    Shim.Parameters.nActiveChannels = size( SpineShim.img, 4 ) ;
    % Parameters for conjugate-gradient optimization
    Shim.Parameters.CG.tolerance     = 1E-6 ;
    Shim.Parameters.CG.maxIterations = 10000 ;    

    %%-----
    % dBdI linear 'Current-to-Field' operator
    Shim.img = SpineShim.img ;
    Shim.Hdr.MaskingImage = SpineShim.mask ;

    Shim.Field = [ ] ;       
    Shim.Model = [ ] ; 
    Shim.Probe = [ ] ; 

end
% =========================================================================
function Shim = extractharmonicfield( Shim )
%EXTRACTHARMONICFIELD
% Extract (smooth) harmonic field via RESHARP (Sun, H. Magn Res Med, 2014)
%
% ------
[localField, reducedMask ] = resharp( Shim.Field.img, ...
                                 Shim.Field.Hdr.MaskingImage, ... 
                                 Shim.Field.getvoxelsize(), ...
                                 2*Shim.Field.getvoxelsize(), ...
                                 0) ;

reducedMask                 = shaver(reducedMask, 1) ;
Shim.Field.img              = reducedMask .* ( Shim.Field.img - localField ) ;
Shim.Field.Hdr.MaskingImage = reducedMask ;

end
% =========================================================================
function Shim = computecurrentsoffset( Shim, ...
                    currentsInspired, currentsExpired, ...
                    pressureInspired, pressureExpired )
%COMPUTECURRENTSOFFSET
% (DC offset // bias)
%
%   Syntax
%
%   Shim = COMPUTECURRENTSOFFSET( Shim, iIn, iEx, pIn, pEx  )
%
%       creates field Shim.Model.currentsOffset
% ------
Specs = ShimSpecs();

A = Shim.getshimoperator ;
M = Shim.gettruncationoperator ;

shimFieldOffset = ...
    A*(currentsInspired*pressureExpired - currentsExpired*pressureInspired)/(pressureInspired - pressureExpired);

% ------- 
% Least-squares solution via conjugate gradients
Shim.Model.currentsOffset = cgls( A'*M'*M*A, ... % least squares operator
                            A'*M'*M*shimFieldOffset, ... % effective solution vector
                            zeros( [Shim.Parameters.nActiveChannels 1] ), ... % initial model (currents) guess
                            Shim.Parameters.CG ) ;

end
% =========================================================================
function Shim = calibraterealtimeupdates( Shim, Params )
%CALIBRATEREALTIMEUPDATESPROBE
%
%   Syntax
%
%   Params.
%       .pressureLogFilenames
%           Cell array of 2 or more filenames to be written for each of the
%           pressure logs recorded during the scans  
%        
%        .threshold 
%           Used for phase unwrapping: voxels are excluded wherever the
%           corresponding magnitude image falls below .threshold times the
%           maximum recorded magnitude value
%           [default : 0.01]
%
%       .maxRunTime 
%           Max duration of pressure recording
%           [default : 30 s]
% ------
fprintf('\n-----\n \t Calibrating respiratory probe.\n')

DEFAULT_MAXRUNTIME = 30 ; % [units : s]
% ------
if  ~myisfield( Params, 'pressureLogFilenames' ) ...
        || max( size( pressureLogFilenames ) ) < 2
    error('Params.pressureLogFilenames must contain at least 2 filenames');
end

if  ~myisfield( Params, 'maxRunTime' ) || isempty( Params.maxRunTime )  
    Params.maxRunTime = DEFAULT_MAXRUNTIME ;
end

if  ~myisfield( Params, 'threshold' ) || isempty( Params.threshold )  
    Params.threshold = DEFAULT_THRESHOLD ;
end

% ------
if  ~myisfield( Shim, 'Probe' ) || isempty( Shim.Probe )  
    Shim.Probe = ProbeTracking ;
end


nCalibrationScans = max( size(Params.pressureLogFilenames) ) ;
pressureLogs      = cell{ nCalibrationScans, 1 } ;
medianPressures   = zeros( nCalibrationScans, 1 ) ;

% -------
% Record pressure    
for iCalibrationScan = 1 : nCalibrationScans
    
    % isUserSatisfied = false ;
    %
    % while ~isUserSatisfied
    %     
    %     disp(' ')
    %     msg = ['Press [Enter] to begin recording calibration pressure log ' ...
    %           num2str(iCalibrationScan) ' of ' num2str(nCalibrationScans) '.']
    %     assert(isempty(msg),'Cancelled calibration.')
    %     
    %     Shim.Probe.Data.pressure = 0 ;
    %     Shim.Probe = Shim.Probe.recordpressurelog( Params ) ;
    %
    %     % ------
    %     % save p-log
    %     pressureLogFid = fopen( ...
    %         Params.pressureLogFilename{iCalibrationScan}, 'w+' ) ;
    %     fwrite( pressureLogFid, Shim.Probe.Data.pressure, 'double' ) ;
    %     fclose( pressureLogFid );
    %
    %     figure ;
    %     plot( pressureLog, '+' ) ;
    %     title( Params.pressureLogFilename{iCalibrationScan} ) ;
    %     xlabel('Sample index');
    %     ylabel('Pressure (0.01 mbar)');
    %     
    %     response = input(['Is the current pressure log satisfactory?' ...
    %         'Enter 0 to rerecord; 1 to continue']) ;
    %
    %     isUserSatisfied = logical(response) ;
    %
    % end
    %
    % pressureLogs{ iCalibrationScan } = pressureLog ;
end

% -------
% Determine median apnea-pressure
for iCalibrationScan = 1 : nCalibrationScans

    fprintf(['\n Determine median pressure during apnea : \n \t ' ...
          num2str(iCalibrationScan) ' of ' num2str(nCalibrationScans) '.\n']) ;
    
    medianPressures{ iCalibrationScan } = ...
        ProbeTracking.userselectmedianpressure(  pressureLogs{ iCalibrationScan } ) ; 

end


% ------
fprintf('\n Processing associated GRE data...\n')

for iCalibrationScan = 1 : nCalibrationScans

    % msg = ['Enter directory 
    % of begin recording calibration pressure log ' ...
    %       num2str(iCalibrationScan) ' of num2str(nCalibrationScans)])
    % assert(isempty(msg),'Cancelled calibration.')
    Params.threshold = 0.01 ; % as percent of max measured intensity.

    Params.pathToShimReferenceMaps = '~/Projects/Shimming/RRI/data/SpineShimReferenceMaps.mat' ;

    Params.dataLoadDir  = '/Users/ryan/Projects/Shimming/Static/SpineShimMrm2016/data/24112015/shim_007/' ;

    Params.Path.Mag.echo1  = [Params.dataLoadDir '03-gre_fieldmapping/echo_4.92' ] ;
    Params.Path.Mag.echo2  = [Params.dataLoadDir '03-gre_fieldmapping/echo_7.38' ] ;

    Params.Path.Phase.echo1 = [Params.dataLoadDir '04-gre_fieldmapping/echo_4.92' ] ;
    Params.Path.Phase.echo2 = [Params.dataLoadDir '04-gre_fieldmapping/echo_7.38' ] ;


    % =========================================================================
    % PROCESS GRE FIELD MAP
    % =========================================================================
    [Field,Extras] = ShimOpt.mapfield( Params )

    Extras.PhaseEcho2.img = Extras.PhaseEcho2.Hdr.MaskingImage .* ( Extras.PhaseEcho2.img + 2*pi ) ;

    Field = Extras.PhaseEcho1.mapfrequencydifference( Extras.PhaseEcho2 ) ;

    Field.Hdr.MaskingImage = Extras.PhaseEcho1.Hdr.MaskingImage .* Extras.PhaseEcho2.Hdr.MaskingImage ; 
end

% ------
fprintf(['\n Determining optimal pressure-to-field coupling coefficients' ...
    '(May take several minutes)'])

% Specs = ShimSpecs();
%
%
% A = Shim.getshimoperator ;
%
% Shim.Model.couplingCoefficients = ...
%     A*(currentsExpired - currentsInspired)/(pressureInspired - pressureExpired) ;
%

end
% =========================================================================
function Shim = setcouplingcoefficients( Shim, ...
                    currentsInspired, currentsExpired, ...
                    pressureInspired, pressureExpired )
%SETCOUPLINGCOEFFICIENTS
%
%   Syntax
%
%   Shim = SETCOUPLINGCOEFFICIENTS( Shim, iIn, iEx, pIn, pEx )
%
%   iIn & iEx are vectors of optimal currents for inspired and expired fields.
%   pIn & pEx the associated pressures (scalars)
%
%       creates field Shim.Model.couplingCoefficients
% ------

A = Shim.getshimoperator ;

Shim.Model.couplingCoefficients = ...
    A*(currentsExpired - currentsInspired)/(pressureInspired - pressureExpired) ;

end
% =========================================================================
function [Shim] = optimizeshimcurrents( Shim, Params )
%OPTIMIZESHIMCURRENTS 
%
%   OPTIMIZESHIMCURRENTS( Shim, Params )
%   
% Params can have the following fields 
%   
%   .maxCurrentPerChannel
%       [default: 4 A,  determined by class ShimSpecs.Amp.maxCurrentPerChannel]
%
% ------
Specs = ShimSpecs();

if nargin < 1
    error('Function requires at least 1 argument of type ShimOpt')
elseif nargin == 1
    Params.dummy = [];
end

if ~myisfield(Params, 'maxCurrentPerChannel') || isempty( Params.maxCurrentPerChannel ) 
    Params.maxCurrentPerChannel = Specs.Amp.maxCurrentPerChannel ; 
end


A = Shim.getshimoperator ;
M = Shim.gettruncationoperator ;

b = M*(-Shim.Field.img(:)) ;

% ------- 
% Least-squares solution via conjugate gradients
Shim.Model.currents = cgls( A'*M'*M*A, ... % least squares operator
                            A'*M'*b, ... % effective solution vector
                            zeros( [Shim.Parameters.nActiveChannels 1] ), ... % initial model (currents) guess
                            Shim.Parameters.CG ) ;

% ------- 
% TODO: CHECK SOL. VECTOR FOR CONDITIONS 
isCurrentSolutionOk = false ;

if ~isCurrentSolutionOk
    
    [X0,X1,X2,X3] = ShimCom.getchanneltobankmatrices( ) ;

    Options = optimset(...
        'DerivativeCheck','off',...
        'GradObj','on',...
        'Display', 'off',... %'iter-detailed',...
        'MaxFunEvals',36000,...
        'TolX',1e-11,...
        'TolCon',1E-8);

    A=M*A;
    tic
    [currents] = fmincon( ...
        @shimcost,...
        ones(Specs.Amp.nActiveChannels,1),...
        [],...
        [],...
        [],...
        [],...
        -Params.maxCurrentPerChannel*ones(Specs.Amp.nActiveChannels,1),...
        Params.maxCurrentPerChannel*ones(Specs.Amp.nActiveChannels,1),...
        @first_order_norm,...
        Options);
    toc

    Shim.Model.currents = currents ;
    Shim = Shim.forwardmodelfield ;
end

function [f, df] = shimcost( currents )
    
     y=A*currents - b;
     f=y'*y;
     df=2*A'*y;
    
end

function [C, Ceq] = first_order_norm( currents )
% C(x) <= 0
% (e.g. x = currents)
%
    C   = 0; 
    Ceq = [];
    waterLevel = 1E-8;
    % split currents up into banks
    % -------
    % bank 0
    i0 = X0*currents;

    % -------
    % bank 1
    i1 = X1*currents;

    % -------
    % bank 2
    i2 = X2*currents;

    % -------
    % bank 3
    i3 = X3*currents;

    % Overall abs current cannot exceed Specs.Amp.maxCurrentPerBank (e.g. 20 A)
    C(1) = sum( abs(i0) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 
    C(2) = sum( abs(i1) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 
    C(3) = sum( abs(i2) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 
    C(4) = sum( abs(i3) + waterLevel ) - Specs.Amp.maxCurrentPerBank ; 

    % pos. current cannot exceed Specs.Amp.maxCurrentPerRail (e.g. + 10 A) 
    C(5) = abs(sum( ((i0>0) .* i0) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(6) = abs(sum( ((i1>0) .* i1) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(7) = abs(sum( ((i2>0) .* i2) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(8) = abs(sum( ((i3>0) .* i3) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    
    % neg. current cannot be below Specs.Amp.maxCurrentPerRail (e.g. - 10 A) 
    C(9)  = abs(sum( ((i0<0) .* i0) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(10) = abs(sum( ((i1<0) .* i1) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(11) = abs(sum( ((i2<0) .* i2) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 
    C(12) = abs(sum( ((i3<0) .* i3) + waterLevel )) - Specs.Amp.maxCurrentPerRail ; 

end

end
% % =========================================================================
% function [] = assessshimvolume( Shim )
%     fprintf( ['\n Comparing Defined Vs. Desired (Image) shim volumes\n'...
%              '(In mm, relative to isocentre)\n...'] )
%
%     [X,Y,Z] = Shim.getvoxelpositions( ) ;
%     X = X( logical(Shim.Hdr.MaskingImage) ) ;
%     Y = Y( logical(Shim.Hdr.MaskingImage) ) ;
%     Z = Z( logical(Shim.Hdr.MaskingImage) ) ;
%     
%     fprintf( ['\n Shim volume defined over ' ...
%                  '\n in Z (read [rows]):    ' num2str( [ min( X(:) ) max( X(:) ) ] ) ...
%                  '\n in X (p.e. [columns]): ' num2str( [ min( Y(:) ) max( Y(:) ) ] ) ...
%                  '\n in Y (slice):          ' num2str( [ min( Z(:) ) max( Z(:) ) ] ) ...
%                  '\n']) ;
%
%     [X,Y,Z] = Shim.Field.getvoxelpositions( ) ;
%     X = X( logical(Shim.Field.Hdr.MaskingImage) ) ;
%     Y = Y( logical(Shim.Field.Hdr.MaskingImage) ) ;
%     Z = Z( logical(Shim.Field.Hdr.MaskingImage) ) ;
%
%     fprintf( ['\n Image volume defined over ' ...
%                  '\n in Z (read [rows]):    ' num2str( [ min( X(:) ) max( X(:) ) ] ) ...
%                  '\n in X (p.e. [columns]): ' num2str( [ min( Y(:) ) max( Y(:) ) ] ) ...
%                  '\n in Y (slice):          ' num2str( [ min( Z(:) ) max( Z(:) ) ] ) ...
%                  '\n']) ;
%
% end
% =========================================================================
function Shim = forwardmodelfield( Shim )
% FORWARDMODELFIELD
%
% Shim = FORWARDMODELFIELD( Shim ) ;
%
%   Predicts shim field (output: Shim.Model.field) for given set of currents 
%   (input: Shim.Model.currents)
    
A = Shim.getshimoperator() ;

Shim.Model.field = reshape( A*Shim.Model.currents, size( Shim.Field.img ) ) ;

end
% =========================================================================
function M = gettruncationoperator( Shim )
% GETTRUNCATIONOPERATOR
%
% M = GETTRUNCATIONOPERATOR( Shim ) ;
%
% Truncation .Hdr.MaskingImageing) operator (e.g. M*b, 'picks out' the voi
% voxels from vector b)

nVoxelsImg = numel( Shim.Field.Hdr.MaskingImage ) ;
nVoxelsVoi = nnz( Shim.Field.Hdr.MaskingImage ) ;

indicesVoi = find( Shim.Field.Hdr.MaskingImage(:) ) ;

M = sparse( [1:nVoxelsVoi], indicesVoi, ones([nVoxelsVoi 1]), nVoxelsVoi, nVoxelsImg ) ;

end
% =========================================================================
function A = getshimoperator( Shim )
% GETSHIMOPERATOR
%
% A = GETSHIMOPERATOR( Shim ) ;
%
%   where A * vectorOfShimCurrents = shimField

nVoxelsImg = Shim.getnumberofvoxels() ;

A = zeros( nVoxelsImg, Shim.Parameters.nActiveChannels ) ; 

for channel = 1 : Shim.Parameters.nActiveChannels
    A(:, channel) = reshape( Shim.img(:,:,:, channel), [nVoxelsImg 1] ) ;
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)

function [Field, Extras] = mapfield( Params )
% MAPFIELD
%
% Field = MAPFIELD( Params ) 
%
% Params must contain the following fields
%   .Path.Mag.echo1
%   .Path.Mag.echo2
%   .Path.Phase.echo1
%   .Path.Phase.echo2
%   .threshold     
%
%       .threshold (as a fraction (<1) of max measured magnitude intensity)
%           determines the phase region to be unwrapped (i.e. areas of 
%           low signal are ignored)
%       
% =========================================================================

% =========================================================================

PhaseEcho1 = MaRdI(  Params.Path.Phase.echo1  ) ;
MagEcho1   = MaRdI(  Params.Path.Mag.echo1  ) ;

PhaseEcho2 = MaRdI(  Params.Path.Phase.echo2  ) ;
MagEcho2   = MaRdI(  Params.Path.Mag.echo2  ) ;

%-------
% Mask 
PhaseEcho1.Hdr.MaskingImage = ( MagEcho1.img ) > Params.threshold ;
PhaseEcho2.Hdr.MaskingImage = ( MagEcho2.img ) > Params.threshold ;

% -------------------------------------------------------------------------
% Unwrap phase 
% (Using Abdul-Rahman method)
PhaseEcho1 = PhaseEcho1.unwrapphase( ) ;
PhaseEcho2 = PhaseEcho2.unwrapphase( ) ;

%-------
% Compute field map
Field = PhaseEcho1.mapfrequencydifference( PhaseEcho2 ) ;
Field.Hdr.MaskingImage = PhaseEcho1.Hdr.MaskingImage .* PhaseEcho2.Hdr.MaskingImage ; 

Extras.PhaseEcho1 = PhaseEcho1 ;
Extras.PhaseEcho2 = PhaseEcho2 ;



end
% =========================================================================

end
% =========================================================================
% =========================================================================

end

