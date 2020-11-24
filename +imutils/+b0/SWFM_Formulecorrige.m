function [b0] = SWFM_Formulecorrige(complVol, TE)
% SWFM computes c0 fieldmaps based on "Selective-Weighted-Field-Mapping"
% method.
%
% _SYNTAX_
% 
% [b0] = SWFM(complVol, TE)
%
% _DESCRIPTION_
%
% _INPUT ARGUMENTS_
%
%    complVol
%      complex 4D data set complVol(x,y,z,t)
%    TE
%      array of TEs in [s]
%
% _OUTPUTS_
%
%   b0
%     field map in units of Hz 
% 

gamma = 267.52218744 * 10^6; % rad*Hz/Tesla
factor= 1./(gamma.*pi.*TE(end));
            
% create magnitude and phase data volumes
mag_data = abs(complVol);
ph_data = angle(complVol);


X=zeros(size(complVol));
for t=1:length(TE)
    X(:,:,:,t)= mag_data(:,:,:,t).*exp(1i*ph_data(:,:,:,t));
end


% SDQWF correspond au "Signal Decay Quality Weighting Factor"
% FM représente la fraction de la magnitude de l'écho j par rapport à la
% somme des magnitudes.
% FB représente la fraction inverse du bruit de l'écho j par rapport à la
% somme des bruits du signal.

% Calcul du facteur FM 
magtotal=zeros( size(mag_data,1), size(mag_data,2), size(mag_data,3)  );
FM=zeros(size(mag_data));
for t=1:length(TE)
    magtotal(:,:,:)=magtotal(:,:,:)+mag_data(:,:,:,t);
end
for t=1:length(TE)
    FM(:,:,:,t)=mag_data(:,:,:,t)./magtotal;
end

    % Calcul du facteur FB présent dans la fonction arctangente à 4 cadrants

sigmatot=0;
FB=zeros(1,size(mag_data,4));
sigma = imutils.b0.Csigma(mag_data);  
for t=1:length(TE)
    sigmatot = sigmatot + sigma(t);
end
for t =1:length(TE)
    FB(t)=sigmatot./sigma(t);
end


% Calcul du facteur SDQWF
SDQWF=zeros(size(mag_data));

for t=1:length(TE)
    SDQWF(:,:,:,t)=FM(:,:,:,t).*FB(t);
end
    
% Calcul du facteur produit X1Xj présent dans la fonction arctangente à 4 cadrants
X1X_T=zeros(size(complVol));
deltaphi=zeros(size(complVol));
for t=1:length(TE)
    deltaphi(:,:,:,t)=ph_data(:,:,:,1)-ph_data(:,:,:,t);      
    X1X_T(:,:,:,t)= mag_data(:,:,:,1).*mag_data(:,:,:,t).*exp(1i.*(-deltaphi(:,:,:,t)));
end

% Calcul du facteur somme présent dans la fonction arctangente à 4 cadrants        

% cond1 = zeros( size(complVol) );
% cond2 = zeros( size(complVol) );
% cond3 = zeros( size(complVol) );
% cond4 = zeros( size(complVol) );
% cond5 = zeros( size(complVol) );
% 
% for t = 1:length(TE)
%     cond1(:,:,:,t) = real(X1X_T(:,:,:,t))>0 ;
%     cond2(:,:,:,t) = (real(X1X_T(:,:,:,t))<0 & imag(X1X_T(:,:,:,t))~=0) ;
%     cond3(:,:,:,t) = (real(X1X_T(:,:,:,t))>0 & imag(X1X_T(:,:,:,t))==0) ;
%     cond4(:,:,:,t) = (real(X1X_T(:,:,:,t))<0 & imag(X1X_T(:,:,:,t))==0) ;
%     cond5(:,:,:,t) = (real(X1X_T(:,:,:,t))==0 & imag(X1X_T(:,:,:,t))==0) ;
% end
% for t = 1 : length(TE)
%     for i= 1:size(complVol,1) 
%         for j= 1:size(complVol,2)
%             for k= 1: size(complVol,3)
%                 if cond1(i,j,k,t)==1
%                      sommable1(i,j,k,t) = (imag(X1X_T(i,j,k,t)).*SDQWF(i,j,k,t))./(sqrt( (imag(X1X_T(i,j,k,t))).^2 + (real(X1X_T(i,j,k,t))).^2 ) +real(X1X_T(i,j,k,t)) );
%                 elseif cond2(i,j,k,t)==1
%                     sommable2(i,j,k,t) = ((sqrt((imag( X1X_T(i,j,k,t) ) ).^2 + (real(X1X_T(i,j,k,t))).^2 ) - real(X1X_T(i,j,k,t))) .*SDQWF(i,j,k,t))./(imag(X1X_T(i,j,k,t)));
%                 end
%             end
%         end
%     end
% end

complexe_imag = zeros( size(complVol) );
complexe_reel = zeros( size(complVol) );

for t = 1:length(TE)
    complexe_reel(:,:,:,t) = real(X1X_T(:,:,:,t)) ;
    complexe_imag(:,:,:,t) = imag(X1X_T(:,:,:,t)) ;
end

b0=zeros( size(complVol,1) ,size(complVol,2) ,size(complVol,3));  
phase=zeros( size(complVol));  
indice=1;
for echo = 2 : length(TE)
    phase(:,:,:,indice)= +imutils.b0.arctan4quadrant(complexe_imag(:,:,:,echo),complexe_reel(:,:,:,echo),SDQWF(:,:,:,echo));
    b0(:,:,:) = b0(:,:,:) + phase(:,:,:,indice);
    indice = indice +1;
end
b0=factor.*b0;
b0=b0/(2*pi);
% Sans la fonction arctan4quadrant qui est la fonction arctangente à 4 quadrants 
% modifiée, on aurait obtenu l'algorithme suivant.

% b0=zeros( size(complVol,1) ,size(complVol,2) ,size(complVol,3));
% for t = 2:length(TE)
%     for i= 1:size(complVol,1) 
%        for j= 1:size(complVol,2)
%            for k= 1: size(complVol,3)
%                if cond1(i,j,k,t)==1
%                    b0(i,j,k)  = b0(i,j,k) + atan(sommable1(i,j,k,t));
%                elseif cond2(i,j,k,t)==1
%                    b0(i,j,k)  = b0(i,j,k) + atan(sommable2(i,j,k,t));
%                elseif cond3(i,j,k,t)==1
%                    b0(i,j,k)  = b0(i,j,k) + pi/2;
%                elseif cond4(i,j,k,t)==1
%                    b0(i,j,k)  = b0(i,j,k) + (-pi/2);
%                elseif cond5(i,j,k,t)==1
%                    b0(i,j,k)  = b0(i,j,k);
%                end
%            end
%        end
%     end
% end
% b0(:,:,:)= factor.*b0(:,:,:);    
% b0 = b0./(2*pi); % [Hz]


