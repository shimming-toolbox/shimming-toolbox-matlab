function [] = assertheadervalidity( Hdrs )
%ASSERTHEADERVALIDITY Issue error if input Hdrs is not as expected  
% 
% [] = ASSERTHEADERVALIDITY( Hdrs )
    
    if nargin ~= 1 || ~isstruct( Hdrs )
        error('Function requires single input argument: Hdrs struct array.' ) ;
    end
    
    %% Check fields exist
    REQUIRED_FIELDS = { 'SeriesInstanceUID'; 
                        'Rows'; 
                        'Columns'; 
                        'PixelSpacing' ; 
                        'ImageOrientationPatient' ; 
                        'ImagePositionPatient' ; 
                      } ;
    
    iMissingFields = ~isfield( Hdrs, REQUIRED_FIELDS ) ; 
    
    if any( iMissingFields )
        missingFields = { REQUIRED_FIELDS{ iMissingFields } } ;
        error( [ 'Required Hdr fields not defined:\n .' strjoin( missingFields, '\n .' ) ], '%s' ) ; 
    end
    
    assert( isfield( Hdrs, 'SliceThickness' ) || ...
            isfield( Hdrs, 'SpacingBetweenSlices' ), ...
        'Expected Hdr entry for at least one of the fields .SliceThickness or .SpacingBetweenSlices' ) ;

    %% Check field entries that should possess a unique value
    
    assert( length( unique( string( { Hdrs(:).SeriesInstanceUID } ) ) ) == 1, ...
        'All elements of the DICOM Headers struct array must belong to the same acquisition series' ) ;

    assert( length(unique( [ Hdrs(:).Rows ] )) == 1, ...
        'Expected unique Hdr entry for .Rows' ) ;
    
    assert( length( unique( [ Hdrs(:).Columns ] ) ) == 1, ...
        'Expected unique Hdr entry for .Columns' ) ;

    % NOTE/TODO
    % The following conditions need to be met for instantiating a MrdiGrid
    % object, however, they may not pass for acquisitions with multiple slice
    % groups (e.g. localizer)...
    
    % pixelSpacing = [ Hdrs(:).PixelSpacing ] ;
    %
    % assert( ( length( unique( pixelSpacing(1,:) ) ) == 1 ) || ...
    %         ( length( unique( pixelSpacing(2,:) ) ) == 1 ), ...
    %     'Expected unique Hdr entry for both .PixelSpacing(1) and .PixelSpacing(2)' ) ;
    %
    % if myisfield( Hdrs, 'SliceThickness' )
    %     assert( length( unique( [ Hdrs(:).SliceThickness ] ) ) == 1, ...
    %         'Expected unique Hdr entry for .SliceThickness' ) ;
    % end
    %
    % if myisfield( Hdrs, 'SpacingBetweenSlices' )
    %     assert( length( unique( [ Hdrs(:).SpacingBetweenSlices ] ) ) == 1, ...
    %         'Expected unique Hdr entry for .SpacingBetweenSlices' ) ;
    % end

end    
