classdef ShimDes < matlab.mixin.SetGet
%SHIMDES - Shim Design 
%
% .......
% 
% Usage
%   
%   TODO:
%       documentation
% .......
%
% Note 
%   
% Contributions:
%
% Kai-Ming Lo, Resonance Research Inc:
%   Original plot coil script
%
% The ___ function is modified from a previous function ('b1sim_dc_core()' from
% Jason P. Stockmann) which was itself adapted from Fa-Hsuan Lin's Biot-Savart
% solver software.
%
% Please cite: 
%
% Fa-Hsuan Lin, "Magnetic field by Biot-Savart's Law"
% http://maki.bme.ntu.edu.tw/?page_id=333
%
%
% =========================================================================
% Updated::20180927::ryan.topfer@polymtl.ca
% =========================================================================

properties
    pumcin ;
    Wires ;
    Params ;
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimDes( Params )
%SHIMDES - Shim Design 

DEFAULT_PUMCINFILENAME             = '/Users/ryan/Projects/Shimming/Acdc7t/scripts/acdc7t.pmcn' ;
DEFAULT_PATHTOSAGITTALREFERENCEIMG = '/Users/ryan/Projects/Shimming/Acdc7t/scripts/SagittalReferenceImg34/';

if nargin == 0
    Params.dummy = [] ;
end

if ~myisfield( Params, 'pumcinFilename' ) || isempty( Params.pumcinFilename ) 
    Params.pumcinFilename = DEFAULT_PUMCINFILENAME  ;
end

Shim.pumcin = load( Params.pumcinFilename ) ;

if ~myisfield( Params, 'pathToSagittalReferenceImg' ) || isempty( Params.pathToReferenceImg ) 
    Params.pathToSagittalReferenceImg = DEFAULT_PATHTOSAGITTALREFERENCEIMG ;
end

Shim.Params = Params ;

Shim.matchcoordinatesystems( Params ) ;

% TODO: 
%   place the following into its own function, e.g. Shim.populatewirepattern() ?
nPoints = size( Shim.pumcin, 1 ) ; 

wireStartPoint = Shim.pumcin(1, 2:5) ;
iChannel            = 0 ;
nSegmentsPerChannel = 0 ;

% If 2 points are neared than 'tolerance' they're considered identical 
% (i.e. loop beginning == loop ending)
tolerance      = 0.001 ; % [units: mm]

for iPoint = 1 : nPoints
        
    if Shim.pumcin( iPoint, 5 ) == 0
        % wire point is a start point
        wireStartPoint = Shim.pumcin( iPoint, 2:4 ) ;
        iChannel = iChannel + 1 ;
       
        iSegment = 1 ; 
        Shim.Wires{iChannel, iSegment}.start = wireStartPoint ;
    
    else
       
        iSegment = iSegment + 1 ;

        Shim.Wires{iChannel, iSegment-1}.stop = Shim.pumcin( iPoint, 2:4 ) ;
            
        if ( norm( Shim.pumcin( iPoint, 2:4 ) - wireStartPoint ) < tolerance )
            % wire point is an end point (same position as start point but non-zero weight (4th column in Shim.pumcin))
            nSegmentsPerChannel(iChannel) = iSegment - 1 ;
        else
            Shim.Wires{iChannel, iSegment}.start  = Shim.pumcin( iPoint, 2:4 ) ;
        end
     end
end

Shim.plotcoilarray() ;

end
% =========================================================================
function nChannels = getnchannels( Shim )
%GETNCHANNELS - Return the # of shim channels
% 
% nChannels = GETNCHANNELS( Shim )

% a new coil element is demarcated by a zero-weight (5th column) in the .pumcin wire pattern:
nChannels = length( Shim.pumcin(:, 5) ) - nnz( Shim.pumcin(:, 5) ) ;

% Right-hand side is the simpler way to determine nChannels. The assertion is to enforce consistency:
assert( nChannels == size( Shim.Wires, 1 ), ...
    'Problem: Shim.pumcin and Shim.wires representations are inconsistent. Which one is correct?' ) ; 

end
% =========================================================================
function nWireSegments = getnwiresegmentsperchannel( Shim )
%GETNWIRESEGMENTSPERCHANNEL - Return the # of wire segments per coil element 
% 
% nWireSegments = GETNWIRESEGMENTSPERCHANNEL( Shim )

[iChannelStart, iChannelStop] = Shim.getchannelpumcinindices( Shim.pumcin ) ;

nWireSegments = [ iChannelStop - iChannelStart ] ;

end
% =========================================================================
function [] = plotcoilarray( Shim )
%PLOTCOILARRAY - Opens new figure to plot coil array (wires)
%
% [] = PLOTCOILARRAY( Shim ) 

nChannels           = Shim.getnchannels() ;
nSegmentsPerChannel = Shim.getnwiresegmentsperchannel() ;

% from http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html
%
% Long Name	Short Name	RGB Triplet
% blue	        b	[0,0,1]
% black	        k	[0,0,0]
% red	        r	[1,0,0]
% green	        g	[0,1,0]
% yellow	    y	[1,1,0]
% cyan	        c	[0,1,1]
% magenta	    m	[1,0,1]
% white	        w	[1,1,1]
%
% coilColors = [ 0 0 1 ; 1 0 0 ; 0 1 0 ; 1 1 0 ; 0 1 1 ; 1 0 1 ] ;
% if nChannels > size( coilColors, 1 )
% .... TODO ( if desired ) : assign other colours besides blue to the coils

figure

for iChannel = 1 : nChannels

    % the average position of a coil element: text label to be positioned here.
    centerOfMass = [0 0 0] ;

    for iSegment = 1 : nSegmentsPerChannel(iChannel)

        h = line( [Shim.Wires{iChannel,iSegment}.start(1); Shim.Wires{iChannel,iSegment}.stop(1)], ...
                  [Shim.Wires{iChannel,iSegment}.start(2); Shim.Wires{iChannel,iSegment}.stop(2)], ...
                  [Shim.Wires{iChannel,iSegment}.start(3); Shim.Wires{iChannel,iSegment}.stop(3)] ) ;
        
        set(h, 'color', [0 0 1]);
        set(h, 'linewidth', 3);
        hold on
        
        centerOfMass = centerOfMass + [Shim.Wires{iChannel,iSegment}.stop(1), Shim.Wires{iChannel,iSegment}.stop(2), Shim.Wires{iChannel,iSegment}.stop(3)] ;

    end

    centerOfMass = centerOfMass/nSegmentsPerChannel(iChannel) ;

    text( centerOfMass(1), centerOfMass(2), centerOfMass(3), num2str(iChannel) )

end

% dbstop in ShimDes at 185
% % TODO : add a slice or 2 of the localizer to the figure
% % hold on
% % slice(fov_x,fov_y,fov_z,b1_effect_total,[],[],30); %z=0 plane; can change arguments to plot different slices
% % caxis([0,1e-8]);
% % Shim.Params.pathToCoronalReferenceImg = '/Users/ryan/Projects/Shimming/Acdc7t/scripts/CoronalReferenceImg34/' ;
% % SagImg = MaRdI( Shim.Params.pathToSagittalReferenceImg ) ;
% % CorImg = MaRdI( Shim.Params.pathToCoronalReferenceImg ) ;
% % [X,Y,Z]= SagImg.getvoxelpositions() ; 
%
% Params.patientId    = '34' ;
% Params.dataLoadDir  = ['~/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_' Params.patientId '/'] ; 
%
% % Reload actual experimental Params struct (saved @time of acquisitions)
% load( [Params.dataLoadDir '20180607T180006-Params.mat'] ) ;
% ImgArray      = cell( 1, 2);
% ImgArray{1,1} = MaRdI('/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_34/04_gre_field_mapping_shim0_ins/echo_4.92/' ) ;
% ImgArray{1,2} = MaRdI('/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_34/06_gre_field_mapping_shim0_ins/echo_4.92/' ) ;
%
% FieldInspired = FieldEval( ImgArray, Params ) ;
%
% [X,Y,Z ] = FieldInspired.getvoxelpositions() ;
% gridSize = FieldInspired.getgridsize() ;
%
% gridSpacing  = ImgArray{1,1}.getfieldofview()./ImgArray{1,1}.getgridsize() ;
%
% startPos = ( -ImgArray{1,1}.getfieldofview() + gridSpacing )/2 ;
% endPos   = ( ImgArray{1,1}.getfieldofview() - gridSpacing )/2 ;
%
% [Y, X, Z] = ndgrid( [ startPos(2):gridSpacing(2):endPos(2)], ...
%                     [ startPos(1):gridSpacing(1):endPos(1)], ...
%                     [ startPos(3):gridSpacing(3):endPos(3)] );
%
% hold on 
% slice( X, Y, Z, ImgArray{1,1}.img, [], Y(1), [] ) ;

hold off

end
% =========================================================================
function [dr] = matchcoordinatesystems( Shim, Params )
%MATCHCOORDINATESYSTEM 
%

if myisfield( Params, 'translations' )
    assert( numel( Params.translations ) == 3, 'Expected 3-component vector [dx dy dz] for Params.translations' )
    % Apply the user-supplied coordinate translations and exit

    dr = Params.translations ;

    % Unused variable in this case but leaving as comment for explanatory purposes:
    % Params.isRepositioningCoilArray = false;     
    
    % apply the translations
    Shim.pumcin(:, 2) = Shim.pumcin(:, 2) + Params.translations(1) ;
    Shim.pumcin(:, 3) = Shim.pumcin(:, 3) + Params.translations(2) ;
    Shim.pumcin(:, 4) = Shim.pumcin(:, 4) + Params.translations(3) ;

    return;
else
    assert( myisfield( Params, 'pathToSagittalReferenceImg' ), 'Missing required argument: Params.pathToSagittalReferenceImg.' )

    Params.isRepositioningCoilArray = true;
end

    
% save dx dy dz for acdc_33 + 34, to use w/future data sets? (translating CAD coordinates automatically)
dx = 0 ;

Img = MaRdI( Params.pathToSagittalReferenceImg ) ;

[X, Y, Z] = Img.getvoxelpositions();

nSegmentsPerChannel = Shim.getnwiresegmentsperchannel() ;

[iChannelStart, iChannelStop] = Shim.getchannelpumcinindices( Shim.pumcin ) ;

% ------
% YZ coordinates for coils 2, 10, and 15 corresponding to the min/max Z point of each coil 
% (i.e. a saggital plane along X is imagined to intersect the 3 coils at these 6 points):
minYZ_CAD = zeros( 3, 2 ) ;
maxYZ_CAD = zeros( 3, 2 ) ;
cosYZ_CAD = zeros( 3, 2 ) ; % cos = Center of Mass

XYZ_Coil2 = Shim.pumcin( iChannelStart(2):iChannelStop(2), 2:4 ) ;
[~, iMinZ_Coil2] = min( XYZ_Coil2(:,3) ) ;
[~, iMaxZ_Coil2] = max( XYZ_Coil2(:,3) ) ;

XYZ_Coil10 = Shim.pumcin( iChannelStart(10):iChannelStop(10), 2:4 ) ;
[~, iMinZ_Coil10] = min( XYZ_Coil10(:,3) ) ;
[~, iMaxZ_Coil10] = max( XYZ_Coil10(:,3) ) ;

XYZ_Coil15 = Shim.pumcin( iChannelStart(15):iChannelStop(15), 2:4 ) ;
[~, iMinZ_Coil15] = min( XYZ_Coil15(:,3) ) ;
[~, iMaxZ_Coil15] = max( XYZ_Coil15(:,3) ) ;

% coil 2
minYZ_CAD( 1, : ) = [ XYZ_Coil2(iMinZ_Coil2, 2) XYZ_Coil2(iMinZ_Coil2, 3) ] ;
maxYZ_CAD( 1, : ) = [ XYZ_Coil2(iMaxZ_Coil2, 2) XYZ_Coil2(iMaxZ_Coil2, 3) ] ;
cosYZ_CAD( 1, : ) = [ mean( [ XYZ_Coil2(iMinZ_Coil2, 2)  XYZ_Coil2(iMaxZ_Coil2, 2) ] ), ...
                      mean( [ XYZ_Coil2(iMinZ_Coil2, 3)  XYZ_Coil2(iMaxZ_Coil2, 3) ] ) ] ;

% coil 10 
minYZ_CAD( 2, : ) = [ XYZ_Coil10(iMinZ_Coil10, 2) XYZ_Coil10(iMinZ_Coil10, 3) ] ;
maxYZ_CAD( 2, : ) = [ XYZ_Coil10(iMaxZ_Coil10, 2) XYZ_Coil10(iMaxZ_Coil10, 3) ] ;
cosYZ_CAD( 2, : ) = [ mean( [ XYZ_Coil10(iMinZ_Coil10, 2)  XYZ_Coil10(iMaxZ_Coil10, 2) ] ), ...
                      mean( [ XYZ_Coil10(iMinZ_Coil10, 3)  XYZ_Coil10(iMaxZ_Coil10, 3) ] ) ] ;

% coil 15
minYZ_CAD( 3, : ) = [ XYZ_Coil15(iMinZ_Coil15, 2) XYZ_Coil15(iMinZ_Coil15, 3) ] ;
maxYZ_CAD( 3, : ) = [ XYZ_Coil15(iMaxZ_Coil15, 2) XYZ_Coil15(iMaxZ_Coil15, 3) ] ;
cosYZ_CAD( 3, : ) = [ mean( [ XYZ_Coil15(iMinZ_Coil15, 2)  XYZ_Coil15(iMaxZ_Coil15, 2) ] ), ...
                      mean( [ XYZ_Coil15(iMinZ_Coil15, 3)  XYZ_Coil15(iMaxZ_Coil15, 3) ] ) ] ;

% approximate surface of the former:
zz_CAD = linspace( min( minYZ_CAD(:,2) ),  max( maxYZ_CAD(:,2) ), 50 ) ;

[z_surfaceOfFormer, iRows] = sort( [ maxYZ_CAD(:, 2 ) ; cosYZ_CAD(:, 2) ] ) ;
tmp = [ maxYZ_CAD(:, 1 ) ; cosYZ_CAD(:, 1) ] ; 
y_surfaceOfFormer = tmp( iRows ) ;

yy_CAD = spline( z_surfaceOfFormer, y_surfaceOfFormer, zz_CAD ) ;

coilLabels = {'2'; '10'; '15'} ;

% the shifted coil coordinates, to be updated by user selection:
minYZ_userFit = minYZ_CAD ;
maxYZ_userFit = maxYZ_CAD ;
cosYZ_userFit = cosYZ_CAD ;

zz_userFit    = zz_CAD ;
yy_userFit    = yy_CAD ;

% coil translations, to be updated by user selection
dy = 0 ;
dz = 0 ;

isUserSatisfied = false ;

if Params.isRepositioningCoilArray
    display( 'Positioning coil (see figure)' )
    display( 'Hint: When repositioning, try to position Coil 10 as near to the subject as possible, without touching.' )
end

while(~isUserSatisfied)
    
    close all
    imagesc( Y(1,:), Z(:,1), Img.img, [50 400] );
    colormap('gray')
    set(gca,'YDir','normal') ; % 'YDir' property is 'reverse' by default, causing the image to display upside down
    axis square
    hold on
    
    % ------ 
    % updated positions based on user selection
    minYZ_userFit(:, 1) = minYZ_userFit(:, 1) + dy ;
    maxYZ_userFit(:, 1) = maxYZ_userFit(:, 1) + dy ;
    
    minYZ_userFit(:, 2) = minYZ_userFit(:, 2) + dz ;
    maxYZ_userFit(:, 2) = maxYZ_userFit(:, 2) + dz ;
    
    cosYZ_userFit(:, 1) = cosYZ_userFit(:, 1) + dy ;
    cosYZ_userFit(:, 2) = cosYZ_userFit(:, 2) + dz ;

    yy_userFit = yy_userFit + dy ;
    zz_userFit = zz_userFit + dz ;

    hold on    
    for iCoil = 1 : 3
        plot( [ minYZ_userFit( iCoil, 1 ) maxYZ_userFit( iCoil, 1 ) ], [ minYZ_userFit( iCoil, 2 ) maxYZ_userFit( iCoil, 2 ) ], 'y' )
        hold on
        plot( cosYZ_userFit(iCoil, 1), cosYZ_userFit(iCoil, 2), 'or' ) ;
        text( cosYZ_userFit(iCoil, 1), cosYZ_userFit(iCoil, 2), ['Coil ' coilLabels{iCoil}], 'Color', [1 1 0] )
        hold on
    end
    
    plot( yy_userFit, zz_userFit, 'r' )
    title( 'Coil placement awaiting user approval via command line.' )

    response = input(['Is the coil position satisfactory? ' ...
        'Enter 0 to reposition the coil; 1 (or enter) to accept & continue: ']) ;

     if ~isempty(response)
        isUserSatisfied = logical(response) ;
        
        if ~isUserSatisfied
            display( 'Reposition the central coil (10) by using the cursor to select a new position (see figure)' )
            title( 'Reposition the central coil (10) by selecting a new position' )
            [y1, z1] = ginput(1) ; 
            
            dy = y1 - cosYZ_userFit(2, 1) ;
            dz = z1 - cosYZ_userFit(2, 2) ;
        end

     else
         isUserSatisfied = true ;

     end

end

title( 'Coil position' )

% the final user-selected coil translations, relative to its original position:
dy = y1 - cosYZ_CAD(2, 1) ;
dz = z1 - cosYZ_CAD(2, 2) ;

dr = [dx dy dz] ;
display( ['Translation vector (from PUMCIN to patient reference coordinates) : ' num2str(dr)] ) ;

% apply the translations
Shim.pumcin(:, 2) = Shim.pumcin(:, 2) + dx ;
Shim.pumcin(:, 3) = Shim.pumcin(:, 3) + dy ;
Shim.pumcin(:, 4) = Shim.pumcin(:, 4) + dz ;

end
% =========================================================================
function [ b1_z ] = generateshimreferencemaps( Shim, Img )
%GENERATESHIMREFERENCEMAPS
%
% dBz = GENERATESHIMREFERENCEMAPS( Shim, Img )
%
% Returns the dBz, the longitudinal (z) component of the magnetic induction field
% due to a unit current circulating within each Shim channel. Img is a MaRdI-type
% reference image used to define the spatial (x-y-z) positions at which to compute
% dBz.
%
% NOTE
%   Img voxels should be isotropically spaced.
assert( nargin == 2, 'Function requires 2 inputs' )

% -------
% Grid for calculation of B1 output field map 
[X ,Y, Z ] = Img.getvoxelpositions() ;
gridSize   = Img.getgridsize() ;

assert( numel( unique(  Img.getvoxelsize() ) ) == 1, 'Img voxels must be isotropic.' )

%%%%%%%%%%% Biot-Savart's law %%%%%%%%%%%%%%
fprintf('GENERATING B FIELD...\n') ;

nChannels = Shim.getnchannels() ;
nSegmentsPerChannel = Shim.getnwiresegmentsperchannel() ;

%% NOTE 
%   Components besides -z might be of interest later, so I'll leave this commented.
% b1_x      = zeros( [ gridSize nChannels ] ) ;
% b1_y      = zeros( [ gridSize nChannels ] ) ;
b1_z      = zeros( [ gridSize nChannels ] ) ;
% b1_effect = zeros( [ gridSize nChannels ] ) ;

for iChannel = 1 : nChannels 
    
    clear coilCurrents ;

    fprintf('Channel [%d]...\n', iChannel);

    for iSegment = 1 : nSegmentsPerChannel( iChannel ) 
        coilCurrents{iSegment} = Shim.Wires{iChannel, iSegment} ;
    end

    [~, ~, Bz] = ShimDes.generatecoilbfield( coilCurrents, X, Y, Z ) ;

    % b1_x(:,:,:,iChannel) = (42.576E6)*Bx ;
    % b1_y(:,:,:,iChannel) = (42.576E6)*By ;
    b1_z(:,:,:,iChannel) = (42.576E6)*Bz ;

    % b1_effect(:,:,:,iChannel) = b1_x(:,:,:,iChannel)+sqrt(-1.0).*b1_y(:,:,:,iChannel);

end

% b1_effect_total = squeeze( sqrt( sum( abs(b1_effect).^2, 4) ) ) ;

end
% =========================================================================

% =========================================================================

end
% =========================================================================
% =========================================================================    
methods( Static = true )
% =========================================================================    
function [ IXYZW ] = cadtopumcin( cadListFilename, dimsToFlip )
%CADTOPUMCIN
% 
% Routine to reformat a text file from AutoCAD describing a set of coil
% geometries and into PUMCIN format.
%
% Syntax
%
%   [ IXYZW ] = CADTOPUMCIN( cadListFilename )
%
%   Reads .txt file cadListFilename and outputs a PUMCIN-formatted version with the same name
%   but with suffix .pmcn rather than .txt
%   IXYZW is PUMCIN-formatted data: 5 columns:
%   
%   [ RowIndex, X-Coordinate, Y-Coordinate, Z-Coordinate, Weight ]
%   
% 
% Preliminary HOW-TO:
%
%   Getting XYZ coordinate tables of the original coil geometries from the original CAD model:
%
% Current method (could be improved):
%
% 1. Open the design in AutoCAD and turn off all layers except for the coils.
%
% 2. Change the QAFLAGS system parameter (this will change the output of the LIST command to exclude 
% output 'Press enter to continue' lines (i.e. the entire output of LIST is printed automatically).
%
% -Enter command: QAFLAGS
% -Set QAFLAGS to 2
% 
%
% 5. Use the LIST command to print the X-Y-Z coordinates of a coil's centerline
% 
% -Desired is the list of the coordinates making up the pre-fit centerline of each coil. 
% To be able to select the centerline specifically, the easiest way I found is just to click and delete the coil,
% which should reveal the centerline on its own.
% -With the centerlines exposed, click on one and type LIST
% -Copy the output (either from the screen output or the log file) following the
% line 'User Data: Fit Points' (e.g. should be approx ~20-30 lines like:
%                                    X = 0.0      , Y = 1.6      , Z = -78.0
%                                    X = 11.8     , Y = 1.6      , Z = -76.6
% -Repeat step 5 for ever coil, simply appending the copied output to the same text file
% (a single coil begins and ends with the same X-Y-Z coordinates, hence defining a completed 'loop')
%
% 
% Optional/Idea for future dev:
%
% Rather than explicitly copying-and-pasting the LIST output for each coil,
% One can turn the AutoCAD log file on with command: LOGFILEON
% +Get the name of the log file with command: LOGFILENAME

assert( nargin >= 1, 'Requires 1 input argument: name of the text file containing the output from AutoCAD LIST).' )
[ loadDir, name, ext ] = fileparts( cadListFilename ) ;
if isempty(loadDir)
    loadDir = './'
end
assert( strcmp( ext, '.txt' ) ) ;
fid = fopen( cadListFilename ) ;
tmp = textscan( fid, '%s %f %s %f %s %f', 'Delimiter', '=') ;
fclose(fid) ;

% 4th column is the weight of the wire segment (i.e = # wire turns, = 1 for single-turn loops)
XYZW = [ tmp{2} tmp{4} tmp{6} ones(size(tmp{2})) ] ;

% -------
% optional: flip axis orientation using flipDim argument
% (e.g. CAD orientation does not correspond to typical patient/dicom
% orientation (i.e. x increasing right to left, y increasing anterior to
% posterior, z increasing inferior to superior)
if nargin < 2 || isempty(dimsToFlip)
   dimsToFlip = [1 1 1] ;
end

XYZW(:,1:3) = [ dimsToFlip(1)*XYZW(:,1) dimsToFlip(2)*XYZW(:,2) dimsToFlip(3)*XYZW(:,3) ] ;

% -------
% Weight = 0 for the first point in every loop

% 1st loop is simple: Begins with row 1:
XYZW( 1, 4 ) = 0 ;

% Remaining loops require search:
nPoints = size( XYZW, 1 ) ; 

iCoil          = 1;
iCoilStart     = 1;
iCoilEnd       = [] ;

% If 2 points are nearer than 'tolerance' they're considered identical (i.e. loop beginning == loop ending)
tolerance = 0.001 ; % [units: mm]

iPoint = 1 ;
startPoint = XYZW( iPoint, 1:3 ) ;

while( iPoint < nPoints )
        
    iPoint = iPoint + 1 ;
    distanceToStartPoint = norm( XYZW( iPoint, 1:3 ) - startPoint(iCoil, :) ) ;

    if distanceToStartPoint < tolerance
       iCoilEnd(iCoil) = iPoint ; 
       iCoil  = iCoil + 1 ;
       iPoint = iPoint + 1 ;
       iCoilStart(iCoil) = iPoint ;
       if iPoint < nPoints
           startPoint( iCoil, : ) = XYZW( iPoint, 1:3 ) ; 
           % 1st point in a loop, indicated by weight = 0
           XYZW( iPoint, 4 ) = 0 ;
       end
    else
        % not a start point: receives weight = 1
           XYZW( iPoint, 4 ) = 1 ;
    end
end   

% add 1st column (row index)
IXYZW  = [ [1:nPoints]' XYZW ] ; 

% print file
fileId = fopen( [ loadDir '/' name '.pmcn'], 'w') ;
fprintf(fileId,'%03d %-12.5f %-12.5f %-12.5f %3.1f\n', IXYZW' );
fclose(fileId);

end
% =========================================================================
function [iChannelStart, iChannelStop] = getchannelpumcinindices( pumcin )
%GETCHANNELPUMCININDICES
%
% [iChannelStart, iChannelStop] = GEtCHANNELPUMCININDICES( pumcin ) 

% nPoints total (some are nonunique)
nPoints       = length( pumcin(:, 5 ) ) ;
iChannelStart = [ find( pumcin(:, 5) == 0 )  ]' ; 
iChannelStop  = [ [ iChannelStart(2:end)' - 1 ]' nPoints ] ;

end
% =========================================================================
function[ Bx, By, Bz ] = generatecoilbfield( current, X, Y, Z )
%GENERATECOILBFIELD - Biot-Savart model of the coil's magnetic induction
%
% [Bx, By, Bz] = GENERATECOILBFIELD(current, X, Y, Z)
%
% code adapted from Fa-Hsuan Lin's Biot-Savart solver software
%
% Please cite: 
%
% Fa-Hsuan Lin, "Magnetic field by Biot-Savart's Law"
% http://maki.bme.ntu.edu.tw/?page_id=333

assert( nargin == 4 )

gridSize = size( X ) ;

assert( all( gridSize == size( Y ) ) && all( gridSize == size( Z ) ) ) ;

XYZ = [X(:) Y(:) Z(:)] ;
nPositions = size(XYZ, 1) ;

fx = zeros(nPositions, 1) ;
fy = zeros(nPositions, 1) ;
fz = zeros(nPositions, 1) ;

nSegments = length( current ) ;

for iSegment = 1 : nSegments

    if(isfield(current{1},'weight'))
        w=current{iSegment}.weight;
    else
        w=1.0;
    end;
    
    a = repmat( norm(current{iSegment}.start - current{iSegment}.stop).^2, [ nPositions, 1 ] ) ;

    b = 2 * sum( repmat( current{iSegment}.stop - current{iSegment}.start, [nPositions, 1] ) ...
               .*( repmat( current{iSegment}.start, [nPositions,1] ) - XYZ ), ...
             2 ) ;

    c=sum((repmat(current{iSegment}.start,[nPositions,1])-XYZ).^2,2);
    
    s1=current{iSegment}.start;
    s1=repmat(s1,[nPositions,1]);
    s2=current{iSegment}.stop;
    s2=repmat(s2,[nPositions,1]);    
    
    px=(s2(:,2)-s1(:,2)).*(s2(:,3)-s1(:,3))-(s2(:,3)-s1(:,3)).*(s2(:,2)-s1(:,2));
    qx=(s2(:,2)-s1(:,2)).*(s1(:,3)-XYZ(:,3))-(s2(:,3)-s1(:,3)).*(s1(:,2)-XYZ(:,2));
    
    fx=fx+integral(px,qx,a,b,c).*w;
    
    py=(s2(:,3)-s1(:,3)).*(s2(:,1)-s1(:,1))-(s2(:,1)-s1(:,1)).*(s2(:,3)-s1(:,3));
    qy=(s2(:,3)-s1(:,3)).*(s1(:,1)-XYZ(:,1))-(s2(:,1)-s1(:,1)).*(s1(:,3)-XYZ(:,3));
    
    fy=fy+integral(py,qy,a,b,c).*w;
    
    pz=(s2(:,1)-s1(:,1)).*(s2(:,2)-s1(:,2))-(s2(:,2)-s1(:,2)).*(s2(:,1)-s1(:,1));
    qz=(s2(:,1)-s1(:,1)).*(s1(:,2)-XYZ(:,2))-(s2(:,2)-s1(:,2)).*(s1(:,1)-XYZ(:,1));
    
    fz=fz+integral(pz,qz,a,b,c).*w;

end;

Bx = reshape( fx, gridSize ) / 1e4 ;
By = reshape( fy, gridSize ) / 1e4 ;
Bz = reshape( fz, gridSize ) / 1e4 ;

function output=integral(p,q,a,b,c)
% integral{(q+p*t)/sqrt(c+b*t+a*t^2).^3}dt, t from 0 to 1

term1=q.*(2.*(2.*a+b)./(4.*a.*c-b.^2)./sqrt(a+b+c)-2.*b./(4.*a.*c-b.^2)./sqrt(c));
term2=p.*(2.*(b+2.*c)./(b.^2-4.*a.*c)./sqrt(a+b+c)-4.*c./(b.^2-4.*a.*c)./sqrt(c));

output=term1+term2;

end

end
% =========================================================================

end
% =========================================================================
% =========================================================================    

end
