function theta = arctan4quadrant(y,x,SDQWF)
% Fonction arctangente modifiée et adaptée au calcul du champ magnétique 
% dans la méthode SWFM - Selective Weight Field Mapping.
% Theta correspond à la phase.
% y , x et complexdata représentent les données complexes obtenues par 
% simulation ou par IRM. SDQWF correspond au Signal Decay Quality Weight Factor.

for i = 1 : size(x,1)
    for j = 1 : size(x,2)
        for k = 1 : size(x,3)
            if x(i,j,k)> 0
                theta(i,j,k) = atan( (y(i,j,k).*SDQWF(i,j,k))./(sqrt( y(i,j,k)^2+x(i,j,k)^2 )+ x(i,j,k)));
            elseif x(i,j,k)<0 && y(i,j,k)~=0
                    theta(i,j,k)= atan( ((sqrt( y(i,j,k)^2+x(i,j,k)^2 )- x(i,j,k)).*SDQWF(i,j,k))./y(i,j,k));
            elseif x(i,j,k) > 0 && y(i,j,k)==0
                    theta(i,j,k) = pi/2;
            elseif x(i,j,k) < 0 && y(i,j,k)==0
                    theta(i,j,k) = -pi/2;
            elseif x(i,j,k) == 0 && y(i,j,k)==0
                    theta(i,j,k)=0;
            end
        end
    end
end


