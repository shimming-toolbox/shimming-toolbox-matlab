function [b0] = SWFM_FF(complVol, TE)
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
sigma=zeros( size(mag_data)  );
FB=zeros(size(mag_data));
sigmatot=zeros( size(mag_data,1), size(mag_data,2), size(mag_data,3)  );
for t=1:length(TE)
    sigma(:,:,:,t) = stdfilt(mag_data(:,:,:,t));
end
for t=1:length(TE)
    sigmatot(:,:,:)=sigmatot(:,:,:) + sigma(:,:,:,t);
end
for t =1:length(TE) 
    FB(:,:,:,t)=sigmatot./sigma(:,:,:,t);
end


% Calcul du facteur SDQWF
SDQWF=zeros(size(mag_data));
for t=1:length(TE)
    SDQWF(:,:,:,t)=FM(:,:,:,t).*FB(:,:,:,t);
end
    
% Calcul du facteur produit X1Xj présent dans la fonction arctangente à 4 cadrants
X1X_T=[];
for t=1:length(TE)
    X1X_T(:,:,:,t)= mag_data(:,:,:,1).*mag_data(:,:,:,t).*exp(1i*(ph_data(:,:,:,1)+ph_data(:,:,:,t)));
end

% Calcul du facteur somme présent dans la fonction arctangente à 4 cadrants        
sommable1=zeros(size(complVol));
sommable2=zeros(size(complVol));
somme1=zeros( size(complVol,1),size(complVol,2) ,size(complVol,3) );
somme2=zeros( size(complVol,1),size(complVol,2) ,size(complVol,3) );
cond1 = zeros( size(complVol) );
cond2 = zeros( size(complVol) );
cond3 = zeros( size(complVol) );
cond4 = zeros( size(complVol) );
cond5 = zeros( size(complVol) );

for t = 1:length(TE)
    cond1(:,:,:,t) = real(X1X_T(:,:,:,t))>0 ;
    cond2(:,:,:,t) = real(X1X_T(:,:,:,t))<0 & imag(X1X_T(:,:,:,t))~=0 ;
    cond3(:,:,:,t) = real(X1X_T(:,:,:,t))>0 & imag(X1X_T(:,:,:,t))==0 ;
    cond4(:,:,:,t) = real(X1X_T(:,:,:,t))<0 & imag(X1X_T(:,:,:,t))==0 ;
    cond5(:,:,:,t) = real(X1X_T(:,:,:,t))==0 & imag(X1X_T(:,:,:,t))==0 ;
end

cond1f = zeros( size(complVol) );
cond2f = zeros( size(complVol) );
cond3f = zeros( size(complVol) );
cond4f = zeros( size(complVol) );
cond5f = zeros( size(complVol) );

for t = 1:length(TE)
    cond1f(:,:,:) = cond1f(:,:,:)&cond1(:,:,:,t);
    cond2f(:,:,:) = cond2f(:,:,:)&cond2(:,:,:,t);
    cond3f(:,:,:) = cond3f(:,:,:)&cond3(:,:,:,t);
    cond4f(:,:,:) = cond4f(:,:,:)&cond4(:,:,:,t);
    cond5f(:,:,:) = cond5f(:,:,:)&cond5(:,:,:,t);
end


for t = 1 : length(TE)
    for i= 1:size(complVol,1) 
        for j= 1:size(complVol,2)
            for k= 1: size(complVol,3)
                if cond1f(i,j,k)==1
                     sommable1(i,j,k,t) = (imag(X1X_T(i,j,k,t)).*SDQWF(i,j,k,t))./(sqrt( (imag(X1X_T(i,j,k,t))).^2 + (real(X1X_T(i,j,k,t))).^2 ) +real(X1X_T(i,j,k,t)) );
                elseif cond2f(i,j,k)==1
                    sommable2(i,j,k,t) = ((sqrt((imag( X1X_T(i,j,k,t) ) ).^2 + (real(X1X_T(i,j,k,t))).^2 ) - real(X1X_T(i,j,k,t))) .*SDQWF(i,j,k,t))./(imag(X1X_T(i,j,k,t)));
                end
            end
        end
    end
end
                    

for t=1:length(TE)
    somme1(:,:,:) = somme1(:,:,:) + sommable1(:,:,:,t);
    somme2(:,:,:) = somme2(:,:,:) + sommable2(:,:,:,t);
end



% % Somme n'a pas que des valeurs comprises entre 180 et -180 degrés. Il faut
% % faire en sorte que ce soit le cas ?!   
for i = 1 : size( somme1, 1 )
      for j = 1 : size( somme1, 2 )
          for k = 1 : size( somme1, 3 )
              if somme1(i,j,k)>180 
                  while somme1(i,j,k)>180
                      somme1(i,j,k)=somme1(i,j,k)-180;
                  end
              elseif somme1(i,j,k)<-180 
                  while somme1(i,j,k)<-180
                      somme1(i,j,k)=somme1(i,j,k)+180;
                  end
              end
          end
      end
  end
  
   for i = 1 : size( somme2, 1 )
      for j = 1 : size( somme2, 2 )
          for k = 1 : size( somme2, 3 )
              if somme2(i,j,k)>180 
                  while somme2(i,j,k)>180
                      somme2(i,j,k)=somme2(i,j,k)-180;
                  end
              elseif somme2(i,j,k)<-180 
                  while somme2(i,j,k)<-180
                      somme2(i,j,k)=somme2(i,j,k)+180;
                  end
              end
          end
      end
   end
 
 % Je possède des valeurs d'angle donnée par "somme". Or, la fonction atan2
% à 4 quadrants prend en entrées la partie imaginaire Y et la partie 
% réelle du nombre complexe dont l'angle est somme. De plus, on va 
% imposer que la norme soit égale à 1.

% Attention !! Les angles donnés par somme sont des angles en degrés ?! 
Xco1=zeros( size(somme1,1) , size(somme1,2) , size(somme1,3)  );
Yco1=zeros( size(somme1,1) , size(somme1,2) , size(somme1,3)  );
Xco2=zeros( size(somme2,1) , size(somme2,2) , size(somme2,3)  );
Yco2=zeros( size(somme2,1) , size(somme2,2) , size(somme2,3)  );
% 1
for i=1:size(somme1,1)
    for j=1:size(somme1,2)
        for k=1:size(somme1,3)
            if somme1(i,j,k)< 90 && somme1(i,j,k)>0
                Xco1(i,j,k)=1./(sqrt(1+(tan(somme1(i,j,k))).^2));
                Yco1(i,j,k)=sqrt( (( tan( somme1(i,j,k) ) ).^2)./(1 + ( tan( somme1(i,j,k) ) ).^2));
            elseif somme1(i,j,k)> 90
                Xco1(i,j,k)=-1./(sqrt(1+(tan(somme1(i,j,k))).^2));
                Yco1(i,j,k)=sqrt( (( tan( somme1(i,j,k) ) ).^2)./(1 + ( tan( somme1(i,j,k) ) ).^2));
            elseif somme1(i,j,k)> -90 && somme1(i,j,k)<0
                Xco1(i,j,k)=1./(sqrt(1+(tan(somme1(i,j,k))).^2));
                Yco1(i,j,k)=-sqrt( (( tan( somme1(i,j,k) ) ).^2)./(1 + ( tan( somme1(i,j,k) ) ).^2));
            elseif somme1(i,j,k)< -90
                Xco1(i,j,k)=-1./(sqrt(1+(tan(somme1(i,j,k))).^2));
                Yco1(i,j,k)=-sqrt( (( tan( somme1(i,j,k) ) ).^2)./(1 + ( tan( somme1(i,j,k) ) ).^2));
            elseif somme1(i,j,k)==0
                Xco1 = 1 ;
                Yco1 = 0 ;
            elseif somme1(i,j,k)==90
                Xco1 = 0 ;
                Yco1 = 1 ;
            elseif somme1(i,j,k)==-90
                Xco1 = 0 ;
                Yco1 = -1 ;
            elseif somme1(i,j,k)==180 || somme1(i,j,k)==-180
                Xco1 = -1 ;
                Yco1 = 0 ;
            end
        end
    end
end
for i=1:size(somme2,1)
    for j=1:size(somme2,2)
        for k=1:size(somme2,3)
            if somme2(i,j,k)< 90 && somme2(i,j,k)>0
                Xco2(i,j,k)=1./(sqrt(1+(tan(somme2(i,j,k))).^2));
                Yco2(i,j,k)=sqrt( (( tan( somme2(i,j,k) ) ).^2)./(1 + ( tan( somme2(i,j,k) ) ).^2));
            elseif somme2(i,j,k)> 90
                Xco2(i,j,k)=-1./(sqrt(1+(tan(somme2(i,j,k))).^2));
                Yco2(i,j,k)=sqrt( (( tan( somme2(i,j,k) ) ).^2)./(1 + ( tan( somme2(i,j,k) ) ).^2));
            elseif somme2(i,j,k)> -90 && somme2(i,j,k)<0
                Xco2(i,j,k)=1./(sqrt(1+(tan(somme2(i,j,k))).^2));
                Yco2(i,j,k)=-sqrt( (( tan( somme2(i,j,k) ) ).^2)./(1 + ( tan( somme2(i,j,k) ) ).^2));
            elseif somme2(i,j,k)< -90
                Xco2(i,j,k)=-1./(sqrt(1+(tan(somme2(i,j,k))).^2));
                Yco2(i,j,k)=-sqrt( (( tan( somme2(i,j,k) ) ).^2)./(1 + ( tan( somme2(i,j,k) ) ).^2));
            elseif somme2(i,j,k)==0
                Xco2 = 1 ;
                Yco2 = 0 ;
            elseif somme2(i,j,k)==90
                Xco2 = 0 ;
                Yco2 = 1 ;
            elseif somme2(i,j,k)==-90
                Xco2 = 0 ;
                Yco2 = -1 ;
            elseif somme2(i,j,k)==180 || somme2(i,j,k)==-180
                Xco2 = -1 ;
                Yco2 = 0 ;
            end
        end
    end
end
b0=zeros( size(complVol,1),size(complVol,2),size(complVol,3)  )
for i= 1:size(complVol,1) 
   for j= 1:size(complVol,2)
       for k= 1: size(complVol,3)
           if cond1f(i,j,k)==1
               b0(i,j,k)  = factor.*atan2(Yco1(i,j,k),Xco1(i,j,k));
           elseif cond2f(i,j,k)==1
               b0(i,j,k)  = factor.*atan2(Yco2(i,j,k),Xco2(i,j,k));
           elseif cond3f(i,j,k)==1
               b0(i,j,k)  = factor.*(pi/2);
           elseif cond4f(i,j,k)==1
               b0(i,j,k)  = factor.*(-pi/2); 
           elseif cond5f(i,j,k)==1
               b0(i,j,k)  =  0;
           end
       end
   end
end   
b0 = b0./(2*pi); % [Hz]


