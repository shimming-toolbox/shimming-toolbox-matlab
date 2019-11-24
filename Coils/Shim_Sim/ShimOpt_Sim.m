classdef ShimOpt_Sim < ShimOpt
%SHIMOPT_SIM - Shim Optimization for a simulated coil design 
%
% ShimOpt_Sim is a ShimOpt subclass. 
%
% NOTE: 
% 
%   ShimOpt_Sim was written for a specific design project: 
%
%   Lopez Rios et al. Integrated AC / DC coil and dipole Tx array for 7T MRI of
%   the spinal cord.  In: Proc.  27th Annu. Meet. ISMRM, Montreal, Canada,
%   2019.  Abstr. 0220.  
%
%   Until now, this has been its only application! Hence, some methods are likely
%   deprecated and others overly specific to that one project. Considerable
%   adaptation is probably in order to render the class up-to-date (i.e. compatible
%   with the rest of the realtime_shimming library) and properly general.
%
%   Nevertheless, a few methods may be of more or less immediate + general applicability, 
%   e.g. see
%
%   ShimOpt_Sim.cadtopumcin( ):
%       A routine to reformat a text file from AutoCAD describing a set of coil
%       geometries into PUMCIN format.
%       
%   ShimOpt_Sim.generatecoilbfield( ) 
%       Biot-Savart model of the coil's magnetic induction, adapted from software from Fa-Hsuan Lin
%
% .......
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
% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

properties  
    pumcin ;
    Params ;
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimOpt_Sim( Params, Field )
%SHIMOPT_SIM - Shim Optimization

if nargin < 1 || isempty( Params ) 
    Params.dummy = [] ;
end

Params     = ShimOpt_Sim.assigndefaultparameters( Params ) ;

Shim.Model  = [] ;
Shim.Aux    = [] ;
Shim.Params = [] ;

Shim.System.Specs    = Params.Specs ; 
Shim.System.currents = zeros( Shim.System.Specs.Amp.nActiveChannels, 1 ) ; 

Params.isGeneratingReferenceMaps = true ;

if Params.isGeneratingReferenceMaps

    DEFAULT_PUMCINFILENAME             = '/Users/ryan/Projects/Shimming/Acdc7t/ismrm2019/scripts/acdc7t.pmcn' ;
    DEFAULT_PATHTOSAGITTALREFERENCEIMG = '';

    if ~myisfieldfilled( Params, 'pumcinFilename' )  
        Params.pumcinFilename = DEFAULT_PUMCINFILENAME  ;
    end

    Shim.pumcin = load( Params.pumcinFilename ) ;

    if ~myisfieldfilled( Params, 'pathToSagittalReferenceImg' )  
        Params.pathToSagittalReferenceImg = DEFAULT_PATHTOSAGITTALREFERENCEIMG ;
    end

    Shim.Params = Params ;

    Shim.matchcoordinatesystems( Params ) ;
    Shim.plotcoilarray() ;

end

%-------
% associate host MRI 
Shim.Aux = ShimOpt_SphericalHarmonics( ) ;

if (nargin == 2) && ~isempty(Field)
    Shim.setoriginalfield( Field ) ;
else
    Shim.Field = [] ;
end


end
% =========================================================================
function [Corrections] = optimizeshimcurrents( Shim, Params )
%OPTIMIZESHIMCURRENTS 
%
% Corrections = OPTIMIZESHIMCURRENTS( Shim, Params )
%   
% Params can have the following fields 
%
%   .maxCorrectionPerChannel
%       [default: determined by ShimSpecs_Sim property: .Amp.maxCurrentPerChannel]
%
%   .minCorrectionPerChannel
%       [default: -.maxCorrectionPerChannel]

if nargin < 2 
    Params.dummy = [];
end

% if ~myisfield( Params, 'maxCorrectionPerChannel') || isempty( Params.maxCorrectionPerChannel ) 
%     Params.maxCorrectionPerChannel = Shim.System.Specs.Amp.maxCurrentPerChannel ; 
% end
%
% if ~myisfield( Params, 'minCorrectionPerChannel') || isempty( Params.minCorrectionPerChannel ) 
%     Params.minCorrectionPerChannel = -Params.maxCorrectionPerChannel ; 
% end

Corrections = optimizeshimcurrents@ShimOpt( Shim, Params ) ;

end
% =========================================================================
function nChannels = getnchannels( Shim )
%GETNCHANNELS - Return the # of shim channels
% 
% nChannels = GETNCHANNELS( Shim )

% a new coil element is demarcated by a zero-weight (5th column) in the .pumcin wire pattern:
nChannels = length( Shim.pumcin(:, 5) ) - nnz( Shim.pumcin(:, 5) ) ;

end
% =========================================================================
function nWireSegments = getnwiresegmentsperchannel( Shim, iChannel )
%GETNWIRESEGMENTSPERCHANNEL - Returns # of wire segments per coil element 
% 
% nWireSegments           = GETNWIRESEGMENTSPERCHANNEL( Shim )
% nWireSegmentsIthChannel = GETNWIRESEGMENTSPERCHANNEL( Shim, iChannel )
%
% If channel index iChannel input argument is provided, the returned scalar
% is the number of wire segments belonging to that channel.
%
% Otherwise (nargin == 1) the returned output is a vector of the number of
% wire segments for each channel in the Shim array.

[iChannelStart, iChannelStop] = Shim.getchannelpumcinindices( Shim.pumcin ) ;

nWireSegments = [ iChannelStop - iChannelStart ] ;

if nargin == 2 && ~isempty( iChannel )
    if isscalar( iChannel ) && iChannel <= Shim.getnchannels()
        nWireSegments = nWireSegments( iChannel ) ;
    else
        error('Invalid channel index.')
    end
end

end
% =========================================================================
function [] = plotcoilarray( Shim, isLabellingCoils )
%PLOTCOILARRAY - Opens new figure to plot coil array (wires)
%
% [] = PLOTCOILARRAY( Shim ) 

nChannels           = Shim.getnchannels() ;
nSegmentsPerChannel = Shim.getnwiresegmentsperchannel() ;

Wires = Shim.getwirepattern() ;

if nargin == 1
    isLabellingCoils = false ;
end
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
coilColors       = repmat([ 0 0 1], [nChannels 1] ) ;
coilColors(1,:)  = [ 1 0 0 ] ;
coilColors(5,:)  = [ 0 1 1 ] ;
coilColors(10,:) = [ 1 1 0 ] ;
coilColors(13,:) = [ 1 0 1 ] ;
coilColors(15,:) = [ 0 1 0 ] ;

% figure

for iChannel = 1 : nChannels


    for iSegment = 1 : nSegmentsPerChannel(iChannel)

        % h = line( [Wires{iChannel,iSegment}.start(1); Wires{iChannel,iSegment}.stop(1)], ...
        %           [Wires{iChannel,iSegment}.start(2); Wires{iChannel,iSegment}.stop(2)], ...
        %           [Wires{iChannel,iSegment}.start(3); Wires{iChannel,iSegment}.stop(3)] ) ;
        h = line( [Wires{iChannel,iSegment}.start(2); Wires{iChannel,iSegment}.stop(2)], ...
                  [Wires{iChannel,iSegment}.start(3); Wires{iChannel,iSegment}.stop(3)], ...
                  [Wires{iChannel,iSegment}.start(1); Wires{iChannel,iSegment}.stop(1)] ) ;
        
        set(h, 'color', coilColors(iChannel,:));
        set(h, 'linewidth', 5);
        hold on

    end

    if isLabellingCoils
        % the average position of a coil element: text label to be positioned here.
        centerOfMass = Shim.getcoilcenterofmass( iChannel ) ;

        text( centerOfMass(2), centerOfMass(3), centerOfMass(1), num2str(iChannel), 'color', coilColors(iChannel,:) )
    end
end

% % TODO : add a slice or 2 of the localizer to the figure
% % hold on
% % slice(fov_x,fov_y,fov_z,b1_effect_total,[],[],30); %z=0 plane; can change arguments to plot different slices
% % caxis([0,1e-8]);
% Shim.Params.pathToCoronalReferenceImg = '/Users/ryan/Projects/Shimming/Acdc7t/scripts/CoronalReferenceImg34/' ;
% SagImg = MaRdI( Shim.Params.pathToSagittalReferenceImg ) ;
% CorImg = MaRdI( Shim.Params.pathToCoronalReferenceImg ) ;
%  [X,Y,Z]= SagImg.getvoxelpositions() ; 

% Params.patientId    = '34' ;
% Params.dataLoadDir  = ['~/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_' Params.patientId '/'] ; 
%
% Reload actual experimental Params struct (saved @time of acquisitions)
% load( [Params.dataLoadDir '20180607T180006-Params.mat'] ) ;
% ImgArray      = cell( 1, 2);
% ImgArray{1,1} = MaRdI('/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_34/04_gre_field_mapping_shim0_ins/echo_4.92/' ) ;
% ImgArray{1,2} = MaRdI('/Users/ryan/Projects/Shimming/Acdc/20180616_Ismrm/data/acdc_34/06_gre_field_mapping_shim0_ins/echo_4.92/' ) ;
%
% FieldInspired = FieldEval( ImgArray, Params ) ;
%
%
%
% Img = MaRdI('/Users/ryan/Projects/Shimming/Acdc7t/data/acdc_spine_7T_06/04-gre_field_mapping_shim/echo_2.46/') ; 
%
% [X,Y,Z ] = Img.getvoxelpositions() ;
% gridSize = Img.getgridsize() ;
%
% gridSpacing  = Img.getfieldofview()./Img.getgridsize() ;
%
% startPos = ( -Img.getfieldofview() + gridSpacing )/2 ;
% endPos   = ( Img.getfieldofview() - gridSpacing )/2 ;
%
% % [X, Y, Z] = ndgrid( [ startPos(2):gridSpacing(2):endPos(2)], ...
% %                     [ startPos(1):gridSpacing(1):endPos(1)], ...
% %                     [ startPos(3):gridSpacing(3):endPos(3)] );
%
% % image is in sagittal orientation, so,
% % 1st (row) dim    = Z (in patient coordinate system)
% % 2nd (column) dim = Y (in PCS)
% % 3rd (slice) dim  = X (in PCS):
% [posR, posC, posS] = ndgrid( [ startPos(3):gridSpacing(3):endPos(3)], ...
%                              [ startPos(1):gridSpacing(1):endPos(1)], ...
%                              [ startPos(2):gridSpacing(2):endPos(2)] );
%
% hold on 
% slice( posS, posS, posC, Img.img, 0, [], [] ) ;
% colormap('gray')
% hold off

end
% =========================================================================
function centerOfMass = getcoilcenterofmass( Shim, iChannel )
%GETCOILCENTEROFMASS - Returns the average position of a coil element
%
% centerOfMass = GETCOILCENTEROFMASS( Shim, iChannel ) 

if ( nargin == 2 ) && isscalar( iChannel ) 
   
    assert( iChannel <= Shim.getnchannels(), 'Invalid channel index' )

    [iChannelStart, iChannelStop] = Shim.getchannelpumcinindices( Shim.pumcin ) ;

    % + 1 to iChannelStart since the corresponding XYZ position is repeated (with a nonzero weight) at iChannelStop:
    XYZW = Shim.pumcin( (iChannelStart(iChannel)+1):iChannelStop(iChannel), 2:5 ) ;
    
    centerOfMass = [0 0 0] ;
    % weighted average:
    centerOfMass(1) = sum( XYZW(:,4) .* XYZW(:,1) )./sum( XYZW(:,4) ) ; % X
    centerOfMass(2) = sum( XYZW(:,4) .* XYZW(:,2) )./sum( XYZW(:,4) ) ; % Y
    centerOfMass(3) = sum( XYZW(:,4) .* XYZW(:,3) )./sum( XYZW(:,4) ) ; % Z

else
    error( 'Invalid input')
end

end
% =========================================================================
function [dr] = matchcoordinatesystems( Shim, Params )
%MATCHCOORDINATESYSTEM 
%

DEFAULT_COILPLACEMENTFILENAME = [ './coilPlacement_' datestr(now,30) ] ;

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

if  ~myisfield( Params, 'coilPlacementFilename' ) || isempty(Params.coilPlacementFilename)
    Params.coilPlacementFilename = DEFAULT_COILPLACEMENTFILENAME ;
end



Img = MaRdI( Params.pathToSagittalReferenceImg ) ;

[X, Y, Z] = Img.getvoxelpositions();

nSegmentsPerChannel = Shim.getnwiresegmentsperchannel() ;

[iChannelStart, iChannelStop] = Shim.getchannelpumcinindices( Shim.pumcin ) ;

% ------
% YZ coordinates for coils 2, 10, and 15 corresponding to the min/max Z point of each coil 
% (i.e. a saggital plane along X is imagined to intersect the 3 coils at these 6 points):
minYZ_CAD = zeros( 3, 2 ) ;
maxYZ_CAD = zeros( 3, 2 ) ;
comYZ_CAD = zeros( 3, 2 ) ; % com = center of mass

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
com = Shim.getcoilcenterofmass( 2 ) ;
comYZ_CAD( 1, : ) = [ com(2) com(3) ] ;

% coil 10 
minYZ_CAD( 2, : ) = [ XYZ_Coil10(iMinZ_Coil10, 2) XYZ_Coil10(iMinZ_Coil10, 3) ] ;
maxYZ_CAD( 2, : ) = [ XYZ_Coil10(iMaxZ_Coil10, 2) XYZ_Coil10(iMaxZ_Coil10, 3) ] ;
com = Shim.getcoilcenterofmass( 10 ) ;
comYZ_CAD( 2, : ) = [ com(2) com(3) ] ;

% coil 15
minYZ_CAD( 3, : ) = [ XYZ_Coil15(iMinZ_Coil15, 2) XYZ_Coil15(iMinZ_Coil15, 3) ] ;
maxYZ_CAD( 3, : ) = [ XYZ_Coil15(iMaxZ_Coil15, 2) XYZ_Coil15(iMaxZ_Coil15, 3) ] ;
com = Shim.getcoilcenterofmass( 15 ) ;
comYZ_CAD( 3, : ) = [ com(2) com(3) ] ;

% approximate surface of the former:
zz_CAD = linspace( min( minYZ_CAD(:,2) ),  max( maxYZ_CAD(:,2) ), 50 ) ;

[z_surfaceOfFormer, iRows] = sort( [ minYZ_CAD(:, 2 ) ; maxYZ_CAD(:, 2 ) ; comYZ_CAD(:, 2) ] ) ;
tmp = [ minYZ_CAD(:, 1 ) ; maxYZ_CAD(:, 1 ) ; comYZ_CAD(:, 1) ] ; 
y_surfaceOfFormer = tmp( iRows ) ;

yy_CAD = spline( z_surfaceOfFormer, y_surfaceOfFormer, zz_CAD ) ;

coilLabels = {'2'; '10'; '15'} ;

% the shifted coil coordinates, to be updated by user selection:
minYZ_userFit = minYZ_CAD ;
maxYZ_userFit = maxYZ_CAD ;
comYZ_userFit = comYZ_CAD ;

zz_userFit    = zz_CAD ;
yy_userFit    = yy_CAD ;

% coil translations :
% dy,dz : to be updated by user selection
dy = 0 ;
dz = 0 ;
% dx : = 0, assumed? or, determine it automatically?
dx = 0 ; 

isUserSatisfied = false ;

if Params.isRepositioningCoilArray
    display( 'Positioning coil (see figure)' )
    display( 'Hint: When repositioning, try to position Coil 10 as near to the subject as possible, without touching.' )
end

while(~isUserSatisfied)
    
    close all
    figure('units','normalized','outerposition',[0 0 1 1])
    imagesc( Y(1,:), Z(:,1), Img.img, [50 200] );
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
    
    comYZ_userFit(:, 1) = comYZ_userFit(:, 1) + dy ;
    comYZ_userFit(:, 2) = comYZ_userFit(:, 2) + dz ;

    yy_userFit = yy_userFit + dy ;
    zz_userFit = zz_userFit + dz ;

    hold on    
    for iCoil = 1 : 3
        % plot( [ minYZ_userFit( iCoil, 1 ) maxYZ_userFit( iCoil, 1 ) ], [ minYZ_userFit( iCoil, 2 ) maxYZ_userFit( iCoil, 2 ) ], 'y' )
        % hold on
        plot( comYZ_userFit(iCoil, 1), comYZ_userFit(iCoil, 2), 'or' ) ;
        text( comYZ_userFit(iCoil, 1), comYZ_userFit(iCoil, 2), ['Coil ' coilLabels{iCoil}], 'Color', [1 1 0] )
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
            
            dy = y1 - comYZ_userFit(2, 1) ;
            dz = z1 - comYZ_userFit(2, 2) ;
        end

     else
         isUserSatisfied = true ;

     end

end

title( 'Coil position' )

% the final user-selected coil translations, relative to its original position:
dy = y1 - comYZ_CAD(2, 1) ;
dz = z1 - comYZ_CAD(2, 2) ;

dr = [dx dy dz] ;
display( ['Translation vector (from PUMCIN to patient reference coordinates) : ' num2str(dr)] ) ;

% apply the translations
Shim.pumcin(:, 2) = Shim.pumcin(:, 2) + dx ;
Shim.pumcin(:, 3) = Shim.pumcin(:, 3) + dy ;
Shim.pumcin(:, 4) = Shim.pumcin(:, 4) + dz ;
    
if ~isempty( Params.coilPlacementFilename )
    close all
    figure('units','normalized','outerposition',[0 0 1 1])
    imagesc( Y(1,:), Z(:,1), Img.img, [50 200] );
    colormap('gray')
    set(gca,'YDir','normal') ; % 'YDir' property is 'reverse' by default, causing the image to display upside down
    axis square
    title( ['CAD translations [X Y Z]: ' num2str(dr)] ) 
    hold on
    Shim.plotcoilarray() ;

    export_fig( Params.coilPlacementFilename, '-png') ;
    pause(1)
    close all
end


end
% =========================================================================
function [] = setoriginalfield( Shim, Field, currents )
%SETORIGINALFIELD 
%
% [] = SETORIGINALFIELD( Shim, Field )
% [] = SETORIGINALFIELD( Shim, Field, currents )
%
% Sets Shim.Field
%
% Field is a FieldEval type object with .img in Hz

if nargin < 2
    error('Not enough input arguments.') ;
elseif nargin == 2
    currents = 0;
    warning('Assuming field map was acquired with all shim channels at 0 A.');
end

Shim.Model.currents = currents ;
Shim.Field = Field.copy() ;

Shim.img = Shim.generateshimreferencemaps( Field ) ;
Shim.Hdr = Field.Hdr ;

Shim.setshimvolumeofinterest( Field.Hdr.MaskingImage ) ;

if ~isempty( Shim.Aux )  
    Shim.Aux.setoriginalfield( Shim.Field ) ;
end


end
% =========================================================================
function [ b1_z ] = generateshimreferencemaps( Shim, Img, Params )
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

% assert( numel( unique(  Img.getvoxelsize() ) ) == 1, 'Img voxels must be isotropic.' )

%%%%%%%%%%% Biot-Savart's law %%%%%%%%%%%%%%
fprintf('GENERATING B FIELD...\n') ;

nChannels = Shim.getnchannels() ;
nSegmentsPerChannel = Shim.getnwiresegmentsperchannel() ;

Wires = Shim.getwirepattern() ;

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
        coilCurrents{iSegment} = Wires{iChannel, iSegment} ;
    end

    [~, ~, Bz] = ShimOpt_Sim.generatecoilbfield( coilCurrents, X, Y, Z ) ;

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
methods(Access=protected)
% =========================================================================
function [ Wires ] = getwirepattern( Shim )
%GETWIREPATTERN
%
% Cycles through Shim.pumcin to define the Wires variable used
% by the Biot-Savart computation ShimOpt_Sim.generatecoilbfield()
%
% [ Wires ] = GETWIREPATTERN( Shim ) ;

nPoints = size( Shim.pumcin, 1 ) ; 

wireStartPoint = Shim.pumcin(1, 2:5) ;
iChannel            = 0 ;
nSegmentsPerChannel = 0 ;

% If 2 points are nearer than 'tolerance' they're considered identical 
% (i.e. loop beginning == loop ending)
tolerance      = 0.001 ; % [units: mm]

for iPoint = 1 : nPoints

    if Shim.pumcin( iPoint, 5 ) == 0
        % wire point is a start point
        wireStartPoint = Shim.pumcin( iPoint, 2:4 ) ;
        iChannel = iChannel + 1 ;

        iSegment = 1 ; 
        Wires{iChannel, iSegment}.start = wireStartPoint ;

    else

        iSegment = iSegment + 1 ;

        Wires{iChannel, iSegment-1}.stop = Shim.pumcin( iPoint, 2:4 ) ;

        if ( norm( Shim.pumcin( iPoint, 2:4 ) - wireStartPoint ) < tolerance )
            % wire point is an end point (same position as start point but non-zero weight (4th column in Shim.pumcin))
            nSegmentsPerChannel(iChannel) = iSegment - 1 ;
        else
            Wires{iChannel, iSegment}.start  = Shim.pumcin( iPoint, 2:4 ) ;
        end
     end
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static=true, Hidden=true)
% =========================================================================
function  [ Params ] = assigndefaultparameters( Params )
%ASSIGNDEFAULTPARAMETERS  
% 
% Params = ASSIGNDEFAULTPARAMETERS( Params )
% 
% Add default parameters fields to Params without replacing values (unless empty)
%
% DEFAULT_PATHTOSHIMREFERENCEMAPS = [] ;
%
% DEFAULT_ISINTERPOLATINGREFERENCEMAPS = true ;

% default: assumes host system is the Prisma @UNF
DEFAULT_INSTITUTIONNAME = 'IUGM' ;
DEFAULT_STATIONNAME     = 'MRC35049' ;

if ~myisfield( Params, 'InstitutionName' ) || isempty(Params.InstitutionName)
   Params.InstitutionName = DEFAULT_INSTITUTIONNAME ;
end

if ~myisfield( Params, 'StationName' ) || isempty(Params.StationName)
   Params.StationName = DEFAULT_STATIONNAME ;
end

end
% =========================================================================

% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
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
% .......
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
