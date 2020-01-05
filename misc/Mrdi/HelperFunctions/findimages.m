function list = findimages( imgDir )
%FINDIMAGES Return list (cell array) of .dcm OR .IMA files within a directory & 'echo_*' subdirectories
%
% list = FINDIMAGES( imageDirectory ) 

imgDir       = [ imgDir '/' ] ;
ListSubdirs  = dir( [ imgDir 'echo*'] );
nEchoSubdirs = length( ListSubdirs ) ;

if nEchoSubdirs > 0 
    
    List = finddcmorima( [imgDir ListSubdirs(1).name] ) ;
    nImgPerEcho = length(List) ;

    list = cell( nImgPerEcho*nEchoSubdirs, 1 ) ;

    for iImg = 1 : nImgPerEcho 
        list{iImg} = [ imgDir ListSubdirs(1).name '/' List(iImg).name ] ;
    end

    for iEcho = 2 : nEchoSubdirs
        ListiSubdir = finddcmorima( [imgDir ListSubdirs(iEcho).name] ) ;
        assert( length(ListiSubdir) == nImgPerEcho, 'Each echo subdirectory should contain the same number of images' ) ; 

        for iImg = 1 : nImgPerEcho
            list{ (iEcho-1)*nImgPerEcho + iImg} = [ imgDir ListSubdirs(iEcho).name '/' ListiSubdir(iImg).name ] ;
        end
    end

else
    List = finddcmorima( imgDir ) ;
    nImg = length(List) ;
    list = cell( nImg, 1 ) ;

    for iImg = 1 : nImg
        list{iImg} = [ imgDir List(iImg).name ] ;
    end
end


function List = finddcmorima( imgDir ) 
%Find .dcm OR .IMA files in imgDir
    List   = dir( [ imgDir '/*.dcm'] );

    if length(List) == 0 % try .IMA
        List = dir( [imgDir '/*.IMA'] ) ;
    end

    nImages = length( List ) ; 
    assert( nImages ~= 0, 'No .dcm or .IMA files found in given directory' ) ;

end

end
