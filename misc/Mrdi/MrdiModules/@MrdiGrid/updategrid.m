function [] = updategrid( Grid, X, Y, Z )
%UPDATEGRID Update image grid positions
% 
% UPDATEGRID( Grid, X, Y, Z ) 

if ( nargin ~= 4 ) || ~isequal( size(X), size(Y), size(Z) ) 
    error( ['Function requires 4 arguments: MrdiGrid object followed by ' ...
        '3 identically-sized position arrays X,Y,Z with regular spacing'] ) ;
end

%% -----
sizeNew = size(X) ;

%% ------
% update spacing

% Displacements [x; y; z;] between rows, columns, and slices, respectively: dR, dC, dS 
dR = [ X(2,1,1) - X(1,1,1) ; Y(2,1,1) - Y(1,1,1) ; Z(2,1,1) - Z(1,1,1) ] ;
dC = [ X(1,2,1) - X(1,1,1) ; Y(1,2,1) - Y(1,1,1) ; Z(1,2,1) - Z(1,1,1) ] ;
    
if sizeNew(3) > 1
    dS = [ X(1,1,2) - X(1,1,1) ; Y(1,1,1) - Y(1,1,2) ; Z(1,1,2) - Z(1,1,1) ] ;
else
    dS = 0 ;
end

spacingNew = sqrt( [ sum(dR.^2) sum(dC.^2) sum(dS.^2) ]  ) ; 

%% -----
% Direction cosines:

% column (expressing angle between column direction and X,Y,Z axes)
imageOrientationPatientNew(4:6) = dR/spacingNew(1) ;
    
% row (expressing angle btw row direction and X,Y,Z axes)
imageOrientationPatientNew(1:3) = dC/spacingNew(2) ;

%% -----
% update imagePositionPatient:
imagePositionPatientNew = [ X(1,1,:) ; Y(1,1,:) ; Z(1,1,:) ] ;

end

