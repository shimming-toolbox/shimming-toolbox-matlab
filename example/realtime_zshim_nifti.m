function realtime_zshim_nifti %(scan_obj, varargin)


[mag_img,mag_info,mag_json] = imutils.read_nii('nifti/echo_2.46_gre_field_mapping_PMUlog_20200313131814_3.nii');
[ph_img,ph_info,ph_json] = imutils.read_nii('nifti/echo_4.92_gre_field_mapping_PMUlog_20200313131814_4_ph.nii');
interpolated_pressure = load('nifti/interpolated_pressure.mat');

ph_img = 2*pi*ph_img/4095;
B0Fields = ph_img./(4.94e-3-2.46e-3); % [rad*Hz]
B0Fields = B0Fields/(2*pi); % [Hz]
  

%B0Fields = FieldEval( FM_mag_path, FM_phase_path, Params ); 

% Siemens PMU recording
%Pmu   = ProbeTracking(respTrace_path);

%% ------------------------------------------------------------------------
% link the two objects (interpolate PMU to the fieldmap time series)
%% ------------------------------------------------------------------------
%B0Fields.associateaux( Pmu );



%% ------------------------------------------------------------------------
% compute z-gradients  
%% ------------------------------------------------------------------------

%GzFields = B0Fields.copy();
GzFields = zeros(size(B0fields));

% scaling factor 
g = 1000/(42.576E6) ; % [units: mT/Hz] 

ImageRes = ph_info.PixelDimensions;
%ImageRes = B0Fields.getvoxelspacing() ; % [units: mm]

for measNo = 1:size(B0Fields,4) %B0Fields.getnumberofmeasurements
    [~,GzFields(:,:,1,measNo)] = gradient( ...
    squeeze(g*B0Fields(:,:,1,measNo)), ImageRes(1,2)/1000, ImageRes(1,3)/1000 ) ; % [units: mT/m]
end



%% ------------------------------------------------------------------------
% plot the Gz map time series
%% ------------------------------------------------------------------------

% figure
% 
% montage(squeeze(GzFields.img));
% caxis([-0.2 0.2])
% colorbar
% title('Gz map time series (mT/m)') ;
% 
% print('-djpeg','Gz_TimeSeries.jpeg');


%% ------------------------------------------------------------------------
% modeled static + respiratory fields (in Field.img and Field.Model.Riro.img respectively)
%% ------------------------------------------------------------------------

%GzField = FieldEval.modelfield( GzFields );


% reproduce Ryan's "modelfield")
Bt = zeros(size(GzFields,1),size(GzFields,2),2);
p_mean = mean( interpolated_pressure ) ;
pressure = interpolated_pressure;% - p_mean ;

% figure 
% 
% for indy = 70%:size(GzFields.img,1)
%     for indx = 1:5:size(GzFields.img,2)
%         if GzFields.Hdr.MaskingImage(indy,indx,1,1,1) == 1 
%             Bt(indy,indx,:) = polyfit(pressure',squeeze(GzFields.img(indy,indx,1,1,:)),1);
%             %Bt(indx,indy,1) = rms(pressure)* Bt(indx,indy,1);
% 
%             GzField.img(indy,indx) = Bt(indy,indx,2);
%             GzField.Model.Riro.img(indy,indx) = Bt(indy,indx,1); 
%             
% %             if (indy == 70) && (indx == 34)
%                 
%                 figure
%                 subplot(2,1,1);
%                 scatter(pressure,squeeze(GzFields.img(indy,indx,1,1,:)));
%                 hold on;
%                 plot(pressure,polyval(squeeze(Bt(indy,indx,:)),pressure),'k--');
%                 hold off;
%                 xl = xlim;
%                 yl = ylim;
%                 xt = 0.05 * (xl(2)-xl(1)) + xl(1);
%                 yt = 0.90 * (yl(2)-yl(1)) + yl(1);
%                 caption = sprintf('y = %f * x + %f', Bt(indy,indx,1), Bt(indy,indx,2));
%                 text(xt, yt, caption, 'FontSize', 16, 'Color', 'r', 'FontWeight', 'bold');
%                 hold off;
%                 
%                 subplot(2,1,2);
%                 plot(squeeze(GzFields.img(indy,indx,1,1,:)));
% %             end
%             
%         end
%     end
% end





for indy = 1:size(GzFields,1)
    for indx = 1:size(GzFields,2)
        if mask(indy,indx,1,1) == 1 
            Bt(indy,indx,:) = polyfit(pressure',squeeze(B0Fields(indy,indx,1,:)),1);
            %Bt(indx,indy,1) = rms(pressure)* Bt(indx,indy,1);

%             B0Field.img(indy,indx) = Bt(indy,indx,2);
%             B0Field.Model.Riro.img(indy,indx) = Bt(indy,indx,1);
                
        end
    end
end



%% ------------------------------------------------------------------------
% plot some results
%% ------------------------------------------------------------------------

figure

subplot(2,1,1);
imagesc( GzField.img ) ;
axis equal
caxis([-0.2 0.2])
colorbar
title('Static Gz [mT/m]') ;

subplot(2,1,2);
imagesc( GzField.Model.Riro.img/GzField.Model.Riro.Aux.Data.p ) ;
axis equal
caxis([-0.0001 0.0001])
colorbar
title('RIRO correction [mT/m/unit-PMU]') ;

%print('-djpeg','Gz_map.jpeg');


%% ------------------------------------------------------------------------
% interpolate field gradient images to target slices for shimming:
%% ------------------------------------------------------------------------
[X,Y,Z]  = Mag.getvoxelpositions() ;

% accelerate interpolation by restricting it to the region where signal exists: 
% (NOTE: mask here could instead be the shimVoi output from SCT)
mask = Mag.getreliabilitymask( 0.05 ) ; % returns mask for each echo (4th dim)
mask = logical( sum( mask, 4 ) ) ; % combine echo-specific masks

% 'interp/extrap' (nearest-neighbour substitution) :
GzField.resliceimg( X,Y,Z, mask ) ; % reslice static b0 image 
GzField.Model.Riro.resliceimg( X,Y,Z, mask ) ; % reslice RIRO image


%% ------------------------------------------------------------------------
% plot some results
%% ------------------------------------------------------------------------

figure 

for ind = 1:1:size(GzField.img,3)
    subplot_tight(2,size(GzField.img,3),ind)
    imagesc( GzField.img(:,:,ind) ) ;
    caxis([-0.2 0.2])
    %colorbar
    set(gca, 'XTickLabel', [],'XTick',[])
    set(gca, 'YTickLabel', [],'YTick',[])
    if ind == 1
        title('resampled static Gz (mT/m)') ;
        set(get(gca,'title'),'Position',[150 0.3 1])
        cb = colorbar('Location','northoutside');
        pos=get(cb,'Position');
        set(cb,'Position',pos+[0.4,0.08,0.1,0]);
    end
end


for ind = 1:1:size(GzField.img,3) 
    subplot_tight(2,size(GzField.img,3),size(GzField.img,3)+ind)
    imagesc( GzField.Model.Riro.img(:,:,ind)/GzField.Model.Riro.Aux.Data.p ) ;
    caxis([-0.0001 0.0001])
    %colorbar
    set(gca, 'XTickLabel', [],'XTick',[])
    set(gca, 'YTickLabel', [],'YTick',[])
    if ind == 1
        title('resampled RIRO correction (mT/m/unit-PMU)') ;
        set(get(gca,'title'),'Position',[150 0.3 1])
                cb = colorbar('Location','northoutside');
        pos=get(cb,'Position');
        set(cb,'Position',pos+[0.4,0.08,0.1,0]);
    end
end

print('-djpeg','resampled_Gz_map.jpeg');

%% ------------------------------------------------------------------------
% Simple (voxelwise) z-shim calculation:
%% ------------------------------------------------------------------------

% Flip static gradient field polarity (aim is to cancel it)
staticTarget = -GzField.img ;

% Flip RIRO gradient polarity + rescale to units of [mT/m/unit-PMU]
riroTarget   = -GzField.Model.Riro.img/GzField.Model.Riro.Aux.Data.p ;

% slicewise corrections within shimVoi (spinal cord volume)
nSlices = size( Mag.img, 3 ) ;

% static slicewise Gz correction [units: mT/m]
Corrections.static = zeros( nSlices, 1 ) ; 
% RIRO slicewise Gz correction [units: mT/m/unit-PMU]
Corrections.riro   = zeros( nSlices, 1 ) ;


for iSlice = 1 : nSlices

    sliceVoi = false( size( shimVoi ) ) ;
    sliceVoi( :,:,iSlice ) = shimVoi(:,:, iSlice ) ;

    Corrections.static( iSlice ) = median( staticTarget( sliceVoi ) ) ;
    Corrections.riro( iSlice )   = median( riroTarget( sliceVoi ) ) ;
    
end

%% ------------------------------------------------------------------------
% write to .txt file readable by sequence
%% ------------------------------------------------------------------------

fileID = fopen('zshim_gradients.txt','w');

for iSlice = 1:(nSlices)
    fprintf(fileID,'Vector_Gz[0][%i]= %.6f\n', iSlice-1, Corrections.static(iSlice)); 
    fprintf(fileID,'Vector_Gz[1][%i]= %.12f\n', iSlice-1, Corrections.riro(iSlice)); 
    fprintf(fileID,'Vector_Gz[2][%i]= %.3f\n', iSlice-1, GzField.Aux.Data.p); 
end

fclose(fileID);


figure

subplot(2,1,1)
hold on;
plot(Corrections.static(:),'b','LineWidth',2)
xlabel('slice number')
ylabel('[mT/m]')
title('static Gz corr');

subplot(2,1,2)
hold on;
plot((max(GzField.Aux.Data.pRaw(:))-GzField.Aux.Data.p)*Corrections.riro(:),'b','LineWidth',2)
xlabel('slice number')
ylabel('[mT/m]')
title('(max pressure - RMS pressure) * RIRO corr');

percent_of_Gzstatic = 100.*((max(GzField.Aux.Data.pRaw(:))-GzField.Aux.Data.p)*Corrections.riro(:))./Corrections.static(:);
for i = 1:size(Corrections.riro(:))
    caption = sprintf('%1.1f %%', percent_of_Gzstatic(i));
    text(i,((max(GzField.Aux.Data.pRaw(:))-GzField.Aux.Data.p)*Corrections.riro(i)),caption);
end

%                 xl = xlim;
%                 yl = ylim;
%                 xt = 0.05 * (xl(2)-xl(1)) + xl(1);
%                 yt = 0.90 * (yl(2)-yl(1)) + yl(1);
%                 caption = sprintf('y = %f * x + %f', Bt(indy,indx,1), Bt(indy,indx,2));
%                 text(xt, yt, caption, 'FontSize', 16, 'Color', 'r', 'FontWeight', 'bold');




% unless in offline processing mode (nargin = 1), copy Dynamic_Gradients.txt to mounted drive
if nargin > 1
   unix('cp zshim_gradients.txt /SYNGO_TRANSFER/SYNGO_TRANSFER/')
end

B0Fields.write('B0_tSeries/','nii',false)
GzFields.write('Gz_tSeries/','nii',false)


% note: Fields.write will convert Fields.Img from double to uint16, which is not
% compatible with various operations in the code


diary off
