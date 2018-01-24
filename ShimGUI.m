function varargout = ShimGUI(varargin)
% SHIMGUI MATLAB code for ShimGUI.fig

%      This function allows to check insp/exp mag and phase data, define
%      some Shim VOIs, optimizing shim currents and send these currents 
%      to the coil through
%      serial communication with an arduino uno board.
%
%      Usage:

%       You choose the directory with the dicomfiles from the scan with the
%       button browse your folder.
%      
%       The function SortData define by the script Sortdata.m is called in background and start to sort the
%       files in a new folder called sorted_data.
%
%       You also need a configuration file which defines the function: 
%       shimparameters(). Example: see: shimparameters.m

%
% Last Modified by GUIDE v2.5 22-Jan-2018 14:42:29


% Begin initialization code - DO NOT EDIT


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ShimGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ShimGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT

%==========================================================================
%   CONFIGURATION PARAMETERS
%==========================================================================

% --- Executes just before ShimGUI is made visible.
function ShimGUI_OpeningFcn(hObject, ~, handles, varargin)


%Call the Configuration file shimparameters.m to specify 
%every parameters needed for shimming simulations.
handles.Params = shimparameters();


% Choose default command line output for ShimGUI
handles.output = hObject;

%Default value for display parameters
handles.item_selected = 'Phase/Inspired';           %View selected
handles.lim=200;                                    %Constrast in Images
handles.limits=[-handles.lim handles.lim];
handles.lim_mag=0.5;
handles.limits_mag=[-handles.lim_mag handles.lim_mag]; 
handles.colormap='parula'; 
handles.n_region=0;                                 %Voi number of regions 
handles.n_region_exclude=0;                           
handles.Voi=[];
handles.Total_Voi=[]; %Voi mask     
handles.Voi_exclude=[];                                
handles.position=cell(1,1);                         %Voi regions positions
handles.position_exclude=cell(1,1);
handles.bound=cell(1,1);
handles.Sct_Voi=[];
handles.PredictedFieldInspired= [];                 %Predicted fields
handles.PredictedFieldExpired= [];
handles.FieldExpired=0;
handles.sliceDeleted=[];                            %Slice excluded from Voi                         
handles.sorted_data=' ';                            %Directory with sorted data
handles.cal_val=[];
handles.send_val=[];
display ('Load your directory with data from the scan');
  


guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = ShimGUI_OutputFcn(~, ~, handles) 
varargout{1} = handles.output;


%==========================================================================
%Options and Commands on the GUI Interface
%==========================================================================


% --- Executes on button press in Load Training Maps :
function Load_training_maps_Callback(hObject, ~, handles)

if ~myisfield(handles, 'Shims')
    handles.Shims = ShimUse(handles.Params) ;
end


%Reverse polarization for the shims =======================================

handles.Shims.Opt.img = -handles.Shims.Opt.img ;

% =========================================================================
% Loading of the Inspired maps
% =========================================================================

ImgArray_Inspired = cell(1 , 2);

%Loading of Magnitude Maps ================================================

display('Select folder with magnitude maps in inspired state')
if handles.sorted_data==' '
    %pathtomaginspired=uigetdir();
    ImgArray_Inspired{1,1} = MaRdI( [ '/Users/ancha_admin/data/sorted_data/08_gre_field_mapping_shim0_ins/ORIGINAL\PRIMARY\M\ND/1' ] ) ;  % ATTENTION HARDCODE FOR TEST
else
handles.pathtomaginspired=uigetdir(handles.sorted_data);
ImgArray_Inspired{1,1} = MaRdI( handles.pathtomaginspired ) ;
end

%Loading of Phase Maps ====================================================

display('Select folder with phase maps in inspired state')
if handles.sorted_data==' '
    %pathtophaseinspired=uigetdir();
    ImgArray_Inspired{1,2} = MaRdI( [ '/Users/ancha_admin/data/sorted_data/10_gre_field_mapping_shim0_ins/ORIGINAL\PRIMARY\P\ND/2' ] ) ;   % ATTENTION HARDCODE FOR TEST
else
handles.pathtophaseinspired = uigetdir(handles.sorted_data);
ImgArray_Inspired{1,2} = MaRdI( handles.pathtophaseinspired ) ;
end


handles.FieldInspired = ShimOpt.mapfield( ImgArray_Inspired, handles.Params ) ;
display('Loading Inspired Fieldmap--------------> Done');

handles.magInspired = ImgArray_Inspired{1,1};
display('Loading Inspired Magnitudes --------------> Done');

%Definition of new variables from the fieldmaps----------------------------

handles.Params.scaling = [min(handles.FieldInspired.img(:)) max(handles.FieldInspired.img(:))] ;
handles.dim = size(handles.FieldInspired.img);
handles.sliceSelected = round(handles.dim(3)*0.5);
handles.clear= zeros(size(handles.FieldInspired.img));

%Set the slice selected and the contrasts on the images--------------------

set(handles.sliceSelector,'value',0.5);
set(handles.Mag_contrast_selector,'value',0.5);
set(handles.Phase_contrast_selector,'value',0.3);
set(handles.commandLine,'string',num2str(handles.sliceSelected));

%Interpolation to the image grid-------------------------------------------

handles.FieldInspired.Hdr.Private_0019_1014 = [0 0 0] ;
handles.Shims.Opt.interpolatetoimggrid( handles.FieldInspired );

%Display Inspired Fieldmap ------------------------------------------------

ImageField(handles.FieldInspired,handles.Fieldmaps,handles);
handles.magInspired.img=handles.magInspired.img-0.5;

Field_parametersins=assessfielddistribution( handles.FieldInspired );
handles.meanins=strcat('Mean Abs=',num2str(Field_parametersins.meanAbs));
handles.medianins=strcat('Median =',num2str(Field_parametersins.median));
handles.stdins=strcat('Std dev =',num2str(Field_parametersins.std));
Setparameters(handles.Field_mean,handles.Field_median,handles.Field_std,handles.meanins,handles.medianins,handles.stdins);

%Display the Feedback from the background script---------------------------

File=fopen('background');
A=fscanf(File,'%c');
disp(A);
fclose(File);

guidata(hObject, handles) ;

% =========================================================================
% Loading of the Expired maps
% =========================================================================

function Load_expired_Callback(hObject, ~, handles)
ImgArray_Expired = cell( 1, 2 ) ;

%Loading of Magnitude Maps ================================================
display('Select folder with magnitude maps in expired state')
handles.pathtomagexpired=uigetdir(handles.sorted_data);
ImgArray_Expired{1,1} = MaRdI( handles.pathtomagexpired ) ;

%Loading of Phase Maps ====================================================
display('Select folder with phase maps in expired state')
handles.pathtophaseexpired=uigetdir(handles.sorted_data);
ImgArray_Expired{1,2} = MaRdI( handles.pathtophaseexpired ) ;

handles.FieldExpired = ShimOpt.mapfield( ImgArray_Expired, handles.Params ) ;
display('Loading Expired Fieldmap  --------------> Done');

handles.magExpired = ImgArray_Expired{1,1};
display('Loading Expired Magnitudes --------------> Done');

%Declaration of the field distribution parameters==========================
Field_parametersexp=assessfielddistribution( handles.FieldExpired );
handles.meanexp=strcat('Mean Abs=',num2str(Field_parametersexp.meanAbs));
handles.medianexp=strcat('Median =',num2str(Field_parametersexp.median));
handles.stdexp=strcat('Std dev =',num2str(Field_parametersexp.std));

handles.magExpired.img=handles.magExpired.img-0.5;

%Feedback from the background script=======================================

File=fopen('background');
A=fscanf(File,'%c');
disp(A);
fclose(File);
guidata(hObject, handles) ;



% --- Executes on button press in customVOI.
function customVOI_Callback(hObject,~, handles)

     if any(handles.sliceSelected == handles.sliceDeleted)>=1
      display('Change the slice selected, this one was removed from the Voi')
     else
            handles.n_region = handles.n_region+1;
% =========================================================================
% Selection of the Regions inside the VOI 
% =========================================================================                         
         if  isempty(handles.Voi)
                handles.rect = imrect(handles.Fieldmaps);
                setColor(handles.rect,'r');  
                handles.Voi =handles.rect.createMask;
                handles.position{1}=handles.rect.getPosition;        
         else        
                handles.rectn = imrect(handles.Fieldmaps);
                setColor(handles.rectn,'r');   
                handles.mask = handles.rectn.createMask;
                handles.Voi=or(handles.Voi,handles.mask);
                handles.position{end+1}=handles.rectn.getPosition;
         end
         
% =========================================================================
% DEFINE Validity mask for shim VOI 
% =========================================================================

            if (handles.FieldExpired ~=0)
            handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.FieldInspired, handles.FieldExpired ) ;
            else
            handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.FieldInspired) ;
            end
% =========================================================================
% Adjust shim VOI based on the rectangular selection on the image
% =========================================================================
            if ~isempty(handles.Voi_exclude);
                handles.Voi = handles.Voi.*handles.Voi_exclude;
            end          
            
            for j=1:handles.dim(3);
                handles.Params.shimVoi(:,:,j)=handles.Params.shimVoi(:,:,j).*handles.Voi;
            end
            
% =========================================================================
% Field distribution parameters inside the VOI
% =========================================================================           
            handles.VOI_parametersins=assessfielddistribution( handles.FieldInspired,handles.Params.shimVoi);
            handles.mean_roi=strcat('Mean Abs=',num2str(handles.VOI_parametersins.meanAbs));
            handles.median_roi=strcat('Median = ',num2str(handles.VOI_parametersins.median));
            handles.std_roi=strcat('Std dev = ',num2str(handles.VOI_parametersins.std));
            
            if (handles.FieldExpired ~=0)
            handles.VOI_parametersexp=assessfielddistribution( handles.FieldExpired,handles.Params.shimVoi);
            handles.mean_roiexp=strcat('Mean Abs=',num2str(handles.VOI_parametersexp.meanAbs));
            handles.median_roiexp=strcat('Median =',num2str(handles.VOI_parametersexp.median));
            handles.std_roiexp=strcat('Std dev =',num2str(handles.VOI_parametersexp.std));
            end
   
             handles.bound=bwboundaries(handles.Voi);
             
% =========================================================================
% Plot Images and field parameters inside the VOI 
% =========================================================================             
       switch handles.item_selected
            
          case 'Phase/Inspired'
              ImageFieldRoi(handles.FieldInspired,handles.Roi,handles.Voi,handles); 
              Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,handles.mean_roi,handles.median_roi,handles.std_roi)
          case 'Phase/Expired'
              ImageFieldRoi(handles.FieldExpired,handles.Roi,handles.Voi,handles);  
              Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,handles.mean_roiexp,handles.median_roiexp,handles.std_roiexp)
          case 'Mag/Inspired'             
              ImageMagRoi(handles.magInspired,handles.Roi,handles.Voi,handles);
              Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,'Mean : -','Median : -','Std dev :-')
          case 'Mag/Expired'
              ImageMagRoi(handles.magExpired,handles.Roi,handles.Voi,handles);
              Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,'Mean : -','Median : -','Std dev :-')
       end
     end
     
% =========================================================================
% Feedback from the backgroud script
% =========================================================================     
File=fopen('background');
A=fscanf(File,'%c');
disp(A);
fclose(File);
guidata(hObject, handles) ;



function Use_sct_Callback(hObject, eventdata, handles)
    
    
    
    
    handles.matlabpath=uigetdir();
    
    command =sprintf('%s', 'sct_propseg -i ',handles.matlabpath,'/gre_field_mapping_shim0_ins.nii',' -c t1 ','-ofolder ',handles.matlabpath,' -CSF');
    command2 =sprintf('%s', 'sct_propseg -i ',handles.matlabpath,'/gre_field_mapping_shim0_exp.nii',' -c t1 ','-ofolder ',handles.matlabpath,' -CSF');
    
    

            dicm2nii(handles.pathtomaginspired,handles.matlabpath,0);        
            
            [~,~] = unix(command);
            handles.sct_mask_ins=load_untouch_nii('/gre_field_mapping_shim0_ins_seg.nii');
            handles.sct_CSFmask_ins=load_untouch_nii('/gre_field_mapping_shim0_ins_CSF_seg.nii');
            
            handles.sctmask_ins = double(handles.sct_mask_ins.img);
            handles.sctCSFmask_ins=double(handles.sct_CSFmask_ins.img);
            
            
            for i=1:handles.dim(3)
            handles.Sct_Voi_ins(:,:,i)=rot90(handles.sctmask_ins(:,:,i));
            handles.Sct_CSF_Voi_ins(:,:,i)=rot90(handles.sctCSFmask_ins(:,:,i));
            end
                  
           
            
            if (handles.FieldExpired ~=0)
            handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.FieldInspired, handles.FieldExpired ) ;   
            dicm2nii(handles.pathtomagexpired,handles.matlabpath,0); 
            [~,~] = unix(command2);
            
            handles.sct_mask_exp=load_untouch_nii('/gre_field_mapping_shim0_exp_seg.nii');
            handles.sct_CSFmask_exp=load_untouch_nii('/gre_field_mapping_shim0_exp_CSF_seg.nii');
            
            handles.sctmask_exp = double(handles.sct_mask_exp.img);
            handles.sctCSFmask_exp=double(handles.sct_CSFmask_exp.img);
            
            
            for i=1:handles.dim(3)
            handles.Sct_Voi_exp(:,:,i)=rot90(handles.sctmask_exp(:,:,i));
            handles.Sct_CSF_Voi_exp(:,:,i)=rot90(handles.sctCSFmask_exp(:,:,i));
            end
            
            handles.Total_Voi=handles.Sct_Voi_ins+handles.Sct_CSF_Voi_ins+handles.Sct_Voi_exp+handles.Sct_CSF_Voi_exp;
            
            else
            handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.FieldInspired) ;
            handles.Total_Voi=handles.Sct_Voi_ins+handles.Sct_CSF_Voi_ins;
            end     
            
            for j=1:handles.dim(3);
                handles.Params.shimVoi(:,:,j)=handles.Params.shimVoi(:,:,j).*handles.Total_Voi(:,:,j);    
                handles.bound{j}=bwboundaries(handles.Total_Voi(:,:,j));        
            end
            
                               
           [handles.mean_roi,handles.median_roi,handles.std_roi,handles.mean_roiexp,handles.median_roiexp,handles.std_roiexp]=calculparameters(handles);

            Switch_plot(handles);         
            
            
guidata(hObject, handles) ;

%Executes on slice selector movement=======================================

function sliceSelector_Callback(hObject,~, handles)

dim = size(handles.FieldInspired.img);

%Slice number in the sagittal plan ========================================
maxslice = dim(3);                   

%Define the slice selected from the slider position
handles.sliceSelected = round(get(hObject,'Value')*(maxslice - 1) + 1);

Switch_plot(handles);

%Display the slice selected on the Shim GUI interface======================

set(handles.commandLine,'string',num2str(round(get(hObject,'Value')*(maxslice - 1) + 1)));
            
guidata(hObject, handles) ;



%View Selected Callback====================================================

function View_selected_Callback(hObject, ~, handles)
    
items = get(hObject,'String');
index_selected = get(hObject,'Value');
handles.item_selected = items{index_selected};

Switch_plot(handles);
            
guidata(hObject, handles);



%Generate Predicted maps callback  ========================================

function Prediction_Callback(hObject, ~, handles)
    
handles.Shims.Opt.setoriginalfield( handles.FieldInspired ) ;
handles.Shims.Opt.setshimvolumeofinterest( handles.Params.shimVoi) ;

%Shim currents optimization ===============================================

if (handles.FieldExpired ~=0)
handles.Params.isSolvingAugmentedSystem    = true ;
handles.Params.isPenalizingFieldDifference = true;
handles.Params.regularizationParameter     = 0 ;
[handles.Params.Inspired.currents, handles.Params.Expired.currents] = handles.Shims.Opt.optimizeshimcurrents( handles.Params, handles.FieldExpired ) ;
else
[handles.Params.Inspired.currents] = handles.Shims.Opt.optimizeshimcurrents(handles.Params) ;
end
    
handles.Shims.Opt.Model.currents   =  handles.Params.Inspired.currents ;

%==========================================================================
%Predicted Inspired field calculation
%==========================================================================

handles.PredictedFieldInspired =  handles.Shims.Opt.predictshimmedfield( ) ;   

%Mask for regions without any signal in the inspired Fieldmaps=============

mask = handles.FieldInspired.img;
mask(mask ~= 0) = 1;


handles.PredictedFieldInspired.img = handles.PredictedFieldInspired.img .* mask;

[handles.mean_pre,handles.median_pre,handles.std_pre,handles.mean_preexp,handles.median_preexp,handles.std_preexp]=calculparameters(handles);


%==========================================================================
%Predicted Expired field calculation
%==========================================================================
if (handles.FieldExpired ~=0)
handles.Shims.Opt.setoriginalfield( handles.FieldExpired ) ;
handles.Shims.Opt.Model.currents = handles.Params.Expired.currents ;
handles.PredictedFieldExpired =  handles.Shims.Opt.predictshimmedfield( ) ;

%Mask for regions without any signal in the expired Fieldmaps==============
mask2 = handles.FieldExpired.img;
mask2(mask2 ~= 0) = 1;
handles.PredictedFieldExpired.img = handles.PredictedFieldExpired.img .* mask2;

end  

%==========================================================================
%Plot Predicted field distribution parameters and Images
%==========================================================================                
switch handles.item_selected
                  
   case 'Phase/Inspired'         
         ImageField(handles.PredictedFieldInspired,handles.Predicted,handles); 
         Setparameters(handles.predictedmean,handles.predictedmedian,handles.predictedstd,handles.mean_pre,handles.median_pre,handles.std_pre)
         Setcurrents(handles.predicted1,handles.predicted2,handles.predicted3,handles.predicted4,handles.predicted5,handles.predicted6,handles.predicted7,handles.predicted8,handles.Params.Inspired.currents);

              
   case 'Phase/Expired'
         ImageField(handles.PredictedFieldExpired,handles.Predicted,handles);     
         Setparameters(handles.predictedmean,handles.predictedmedian,handles.predictedstd,handles.mean_preexp,handles.median_preexp,handles.std_preexp)
         Setcurrents(handles.predicted1,handles.predicted2,handles.predicted3,handles.predicted4,handles.predicted5,handles.predicted6,handles.predicted7,handles.predicted8,handles.Params.Expired.currents);

                
   case 'Mag/Inspired' 
         ImageField(handles.PredictedFieldInspired,handles.Predicted,handles);
         Setparameters(handles.predictedmean,handles.predictedmedian,handles.predictedstd,handles.mean_pre,handles.median_pre,handles.std_pre)
         Setcurrents(handles.predicted1,handles.predicted2,handles.predicted3,handles.predicted4,handles.predicted5,handles.predicted6,handles.predicted7,handles.predicted8,handles.Params.Inspired.currents);

                               
   case 'Mag/Expired'             
         ImageField(handles.PredictedFieldExpired,handles.Predicted,handles);  
         Setparameters(handles.predictedmean,handles.predictedmedian,handles.predictedstd,handles.mean_preexp,handles.median_preexp,handles.std_preexp)
         Setcurrents(handles.predicted1,handles.predicted2,handles.predicted3,handles.predicted4,handles.predicted5,handles.predicted6,handles.predicted7,handles.predicted8,handles.Params.Inspired.currents);
                                     
end


     if any(handles.sliceSelected == handles.sliceDeleted)==0
         axes(handles.Predicted);
            hold on
            if handles.n_region >= 1
             for k = 1:length(handles.bound)
               boundary = handles.bound{k};     
               plot(handles.Predicted,boundary(:,2),boundary(:,1), 'Color','black', 'LineWidth', 1); 
             end
            else
                handles.sct_seg = handles.bound{handles.sliceSelected};
              for k = 1:length(handles.sct_seg)                                   
               boundary = handles.sct_seg{k};     
               plot(handles.Predicted,boundary(:,2),boundary(:,1), 'Color','black', 'LineWidth', 0.8);
              end
            end
     end
guidata(hObject, handles);


% Clear ROI on the selected slice callback=================================

function delete_region_Callback(hObject, ~, handles)

%Modification du shim VOI--------------------------------------------------

handles.Params.shimVoi(:,:,handles.sliceSelected)=0;
handles.sliceDeleted(end+1)=handles.sliceSelected;

%Calculation of the field distribution parameters in the modified VOI------

      [handles.mean_roi,handles.median_roi,handles.std_roi,handles.mean_roiexp,handles.median_roiexp,handles.std_roiexp]=calculparameters(handles);


%Plot field distribution parameters and images-----------------------------

switch handles.item_selected
    
    case 'Phase/Inspired'
        ImageField(handles.FieldInspired,handles.Fieldmaps,handles);
        Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,handles.mean_roi,handles.median_roi,handles.std_roi);
        if handles.n_region >=1
            Imageclear(handles.limits,handles.Roi,handles.colormap,handles);
        end           
           
    case 'Phase/Expired'
         ImageField(handles.FieldExpired,handles.Fieldmaps,handles);
         Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,handles.mean_roiexp,handles.median_roiexp,handles.std_roiexp);
         if handles.n_region >=1
             Imageclear(handles.limits,handles.Roi,handles.colormap,handles);
         end           
    case 'Mag/Inspired'
         ImageMag(handles.magInspired,handles.Fieldmaps,handles);
         Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,'Mean Abs: -','Median : -','Std dev :-');
         if handles.n_region >=1
            Imageclear(handles.limits_mag,handles.Roi,gray,handles);
         end 
                 
    case 'Mag/Expired'
         ImageMag(handles.magExpired,handles.Fieldmaps,handles);
         Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,'Mean Abs: -','Median : -','Std dev :-');
         if handles.n_region >=1
            Imageclear(handles.limits_mag,handles.Roi,gray,handles);
         end
end

if handles.PredictedFieldExpired ~=0
   Imageclear(handles.limits,handles.Predicted,handles.colormap,handles);
end
guidata(hObject,handles);



% Clear VOI callback=======================================================

function Clear_Roi_Callback(hObject, ~, handles)

    switch handles.item_selected
            
         case 'Phase/Inspired'
             ImageField(handles.FieldInspired,handles.Fieldmaps,handles); 
             Imageclear(handles.limits,handles.Roi,handles.colormap,handles);

            
         case 'Phase/Expired'
             ImageField(handles.FieldExpired,handles.Fieldmaps,handles);
             Imageclear(handles.limits,handles.Roi,handles.colormap,handles);
             
                         
         case 'Mag/Inspired'
             ImageMag(handles.magInspired,handles.Fieldmaps,handles) ;
             Imageclear(handles.limits_mag,handles.Roi,gray,handles);
             
                 
         case 'Mag/Expired'
             ImageMag(handles.magExpired,handles.Fieldmaps,handles);
             Imageclear(handles.limits_mag,handles.Roi,gray,handles);
             
             
    end
     
     if handles.PredictedFieldInspired ~=0
        Imageclear(handles.limits,handles.Predicted,handles.colormap,handles);
     end
     
     
%==========================================================================
% reset_current VOI parameters
%==========================================================================

handles.Voi = [];
handles.Total_Voi=[];
handles.Sct_Voi=[];
handles.Voi_exclude=[];
handles.VOI_parametersins=0;
handles.VOI_parametersexp=0;
handles.position=cell(1,1);
handles.position_exclude=cell(1,1);
handles.n_region=0;
handles.n_region_exclude=0;
handles.Params.shimVoi(:,:,:)=0;
handles.PredictedFieldInspired=[];
handles.PredictedFieldExpired=[];
handles.sliceDeleted=[];
handles.bound=[];
Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,'Mean Abs: -','Median : -','Std dev :-');
Setparameters(handles.predictedmean,handles.predictedmedian,handles.predictedstd,'Mean Abs: -','Median : -','Std dev :-');

 guidata(hObject, handles);
 
 % Exclude_From_VOI callback===============================================
 
function Exclude_From_ROI_Callback(hObject, ~, handles)
    
    
     if ~isempty(handles.Voi)
            handles.n_region_exclude = handles.n_region_exclude+1;

% =========================================================================
% Selection of Regions to exclude from the VOI 
% =========================================================================

         if  isempty(handles.Voi_exclude)
            handles.rect = imrect(handles.Fieldmaps);
            setColor(handles.rect,'b');  
            handles.Voi_exclude =~handles.rect.createMask;
            handles.position_exclude{1}=handles.rect.getPosition;       
         else  
            handles.rectn = imrect(handles.Fieldmaps);
            setColor(handles.rectn,'b');
            
            handles.mask = ~handles.rectn.createMask;
            handles.Voi_exclude=handles.Voi_exclude.*handles.mask;
            handles.position_exclude{end+1}=handles.rectn.getPosition;
         end
    
% =========================================================================
% DEFINE SHIM VOI validity mask
% =========================================================================

    handles.Voi = handles.Voi.*handles.Voi_exclude;
    
    if (handles.FieldExpired ~=0)
    handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.FieldInspired, handles.FieldExpired ) ;
    else
    handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.FieldInspired) ;
    end
    handles.bound=bwboundaries(handles.Voi);

% Adjust shim VOI based on the rectangular selection on the image==========

            for j=1:handles.dim(3);
                handles.Params.shimVoi(:,:,j)=handles.Params.shimVoi(:,:,j).*handles.Voi;
            end
            
% =========================================================================
% Field distribution parameters inside the VOI
% =========================================================================
      [handles.mean_roi,handles.median_roi,handles.std_roi,handles.mean_roiexp,handles.median_roiexp,handles.std_roiexp]=calculparameters(handles);

      
% =========================================================================    
%Plot Images and field parameters inside the VOI
% =========================================================================
        switch handles.item_selected
            
                case 'Phase/Inspired'
                    ImageFieldRoi(handles.FieldInspired,handles.Roi,handles.Voi,handles);  
                    Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,handles.mean_roi,handles.median_roi,handles.std_roi)
                case 'Phase/Expired'
                    ImageFieldRoi(handles.FieldExpired,handles.Roi,handles.Voi,handles);   
                    Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,handles.mean_roiexp,handles.median_roiexp,handles.std_roiexp)

                case 'Mag/Inspired'             
                    ImageMagRoi(handles.magInspired,handles.Roi,handles.Voi,handles);
                    Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,'Mean Abs: -','Median : -','Std dev :-')
                case 'Mag/Expired'
                    ImageMagRoi(handles.magExpired,handles.Roi,handles.Voi,handles);
                    Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,'Mean Abs: -','Median : -','Std dev :-')
         end
     else
         display ('Warning : Create a VOI before using the button "Exclude"')
     end   
    
guidata (hObject,handles);

% --- Executes on button press in Extend_SCT.
function Extend_SCT_Callback(hObject, ~, handles)
    
    if ~isempty(handles.Total_Voi)
        handles.extension = imrect(handles.Fieldmaps);
        handles.Sct_extension =handles.extension.createMask;
        handles.Total_Voi(:,:,handles.sliceSelected)=or(handles.Total_Voi(:,:,handles.sliceSelected),handles.Sct_extension);
        
        for j=1:handles.dim(3);
         handles.Params.shimVoi(:,:,j)=handles.Params.shimVoi(:,:,j)+handles.Sct_extension;    
         handles.bound{j}=bwboundaries(handles.Total_Voi(:,:,j));        
        end
        
      [handles.mean_roi,handles.median_roi,handles.std_roi,handles.mean_roiexp,handles.median_roiexp,handles.std_roiexp]=calculparameters(handles);

      
% =========================================================================    
%Plot Images and field parameters inside the VOI
% =========================================================================
        switch handles.item_selected
            
                case 'Phase/Inspired'
                    ImageFieldRoi(handles.FieldInspired,handles.Roi,handles.Total_Voi(:,:,handles.sliceSelected),handles);  
                    Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,handles.mean_roi,handles.median_roi,handles.std_roi)
                case 'Phase/Expired'
                    ImageFieldRoi(handles.FieldExpired,handles.Roi,handles.Total_Voi(:,:,handles.sliceSelected),handles);   
                    Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,handles.mean_roiexp,handles.median_roiexp,handles.std_roiexp)

                case 'Mag/Inspired'             
                    ImageMagRoi(handles.magInspired,handles.Roi,handles.Total_Voi(:,:,handles.sliceSelected),handles);
                    Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,'Mean Abs: -','Median : -','Std dev :-')
                case 'Mag/Expired'
                    ImageMagRoi(handles.magExpired,handles.Roi,handles.Total_Voi(:,:,handles.sliceSelected),handles);
                    Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,'Mean Abs: -','Median : -','Std dev :-')
         end
    
    else
        display('Use Sct button to calculate a VOI before trying to extend it')
    end

guidata (hObject,handles);



% --- Executes on button press in Reduce_SCT.
function Reduce_SCT_Callback(hObject, ~, handles)
    
    if ~isempty(handles.Total_Voi)
        handles.reduction = imrect(handles.Fieldmaps);
        handles.Sct_reduction =~handles.reduction.createMask;
        handles.Total_Voi(:,:,handles.sliceSelected)=handles.Total_Voi(:,:,handles.sliceSelected).*handles.Sct_reduction;
        
        for j=1:handles.dim(3);
         handles.Params.shimVoi(:,:,j)=handles.Params.shimVoi(:,:,j).*handles.Sct_reduction;    
         handles.bound{j}=bwboundaries(handles.Total_Voi(:,:,j));        
        end
        

      [handles.mean_roi,handles.median_roi,handles.std_roi,handles.mean_roiexp,handles.median_roiexp,handles.std_roiexp]=calculparameters(handles);

      
% =========================================================================    
%Plot Images and field parameters inside the VOI
% =========================================================================
        switch handles.item_selected
            
                case 'Phase/Inspired'
                    ImageFieldRoi(handles.FieldInspired,handles.Roi,handles.Total_Voi(:,:,handles.sliceSelected),handles);  
                    Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,handles.mean_roi,handles.median_roi,handles.std_roi)
                case 'Phase/Expired'
                    ImageFieldRoi(handles.FieldExpired,handles.Roi,handles.Total_Voi(:,:,handles.sliceSelected),handles);   
                    Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,handles.mean_roiexp,handles.median_roiexp,handles.std_roiexp)

                case 'Mag/Inspired'             
                    ImageMagRoi(handles.magInspired,handles.Roi,handles.Total_Voi(:,:,handles.sliceSelected),handles);
                    Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,'Mean Abs: -','Median : -','Std dev :-')
                case 'Mag/Expired'
                    ImageMagRoi(handles.magExpired,handles.Roi,handles.Total_Voi(:,:,handles.sliceSelected),handles);
                    Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,'Mean Abs: -','Median : -','Std dev :-')
         end
    
    else
        display('Use Sct button to calculate a VOI before trying to Reduce it')
    end
    
    guidata (hObject,handles);



%Fieldmaps contrast selector callback======================================

function Phase_contrast_selector_Callback(hObject, ~, handles)

handles.min_contrast=max(handles.FieldInspired.img(:));
handles.lim = round(get(hObject,'Value')*(handles.min_contrast - 1) + 1);
handles.limits=[-handles.lim handles.lim];

Switch_plot(handles);

guidata (hObject,handles);


%Magntitude contrast selector callback=====================================

function Mag_contrast_selector_Callback(hObject, ~, handles)

handles.lim_mag = get(hObject,'Value');
if handles.lim_mag <=0.01
    handles.lim_mag=0.011;
end

handles.limits_mag=[-handles.lim_mag handles.lim_mag];

Switch_plot(handles);

guidata (hObject,handles);

% Browse your folder callback==============================================

function Browsebutton_Callback(hObject, ~, handles)
    
%Select the folder with the raw data from the scan-------------------------
    
handles.Raw_data= uigetdir();
[handles.pathtofolder,handles.foldername]=fileparts(handles.Raw_data);

%Define a new folder to store the sorted data------------------------------

handles.sorted_data=fullfile(handles.pathtofolder,'/sorted_data');

%Call the function to sort the data in background--------------------------

unix(sprintf('%s','/Applications/MATLAB_R2016a.app/bin/matlab -nodesktop -nojvm -r ", SortData(''',handles.Raw_data,''' ,''',handles.sorted_data,''');" &'));
display('Load your training maps')
guidata (hObject,handles);


% Send_current Callback----------------------------------------------------

function Send_current_Callback(hObject, ~, handles)
display(handles.Shims.Com.ComPort)

handles.Shims.Com.setandloadallshims(handles.send_val);

guidata (hObject,handles);


% Reset current Callback===================================================

function reset_current_Callback(hObject, ~, handles)
    
    handles.Shims.Com.resetallshims();

guidata (hObject,handles);

% Start Communication Callback=============================================

function start_comunication_Callback(hObject, ~, handles)

switch handles.item_selected
       case {'Phase/Inspired','Mag/Inspired'}
            handles.currents=handles.Params.Inspired.currents.*1000;
       case {'Phase/Expired','Mag/Expired'}
            handles.currents=handles.Params.Expired.currents.*1000;      
end
   
  for i=1:8 
  handles.cal_val(i) = (handles.currents(i) - handles.Params.coeffP2(i)) / handles.Params.coeffP1(i);
  end
   
  handles.send_val = ((1.25 - handles.cal_val * 0.001 * 0.22) * 26214);
  
 handles.Shims.Com.Open_ComPort();
  
  display('Ready to send current')
  display(handles.send_val)
  
guidata (hObject,handles);


%End Communication Callback================================================
    
function end_comunication_Callback(hObject, ~, handles)
handles.Shims.Com.close_ComPort();
display('Serial communication with the coil is close now')
guidata (hObject,handles);


% =========================================================================    
%Functions Image
% =========================================================================

function ImageField(Field,ax,handles)
        
             imageToPlot = Field.img(:,:,handles.sliceSelected);
             imagesc(imageToPlot,'parent',ax);
             colormap(ax,handles.colormap);
             caxis(ax,handles.limits);
             colorbar(ax);
             
function ImageFieldRoi(Field,ax,mask,handles)
        
             imageToPlot = Field.img(:,:,handles.sliceSelected).*mask;
             imagesc(imageToPlot,'parent',ax);
             colormap(ax,handles.colormap);
             caxis(ax,handles.limits);
             colorbar(ax);
             
function ImageMag(Field,ax,handles)
        
             imageToPlot = Field.img(:,:,handles.sliceSelected);
             imagesc(imageToPlot,'parent',ax);
             colormap(ax,gray);
             caxis(ax,handles.limits_mag);
             colorbar(ax);
             
function ImageMagRoi(Field,ax,mask,handles)
             imageToPlot = Field.img(:,:,handles.sliceSelected).*mask;
             imagesc(imageToPlot,'parent',ax);
             colormap(ax,gray);
             caxis(ax,handles.limits_mag);
             colorbar(ax);
   
function Imageclear(lim,ax,colormaps,handles)
             imagesc(handles.clear(:,:,handles.sliceSelected),'parent',ax);
             colormap(ax,colormaps);
             caxis(ax,lim);
             colorbar(ax);
                           
function[a,b,c,d,e,f] = calculparameters(handles) 
    
      handles.VOI_parametersins=assessfielddistribution( handles.FieldInspired,handles.Params.shimVoi);
      handles.mean_roi=strcat('Mean Abs=',num2str(handles.VOI_parametersins.meanAbs));
      a=handles.mean_roi;
      handles.median_roi=strcat('Median =',num2str(handles.VOI_parametersins.median));
      b=handles.median_roi;
      handles.std_roi=strcat('Std dev = ',num2str(handles.VOI_parametersins.std));
      c=handles.std_roi;
      d=0;
      e=0;
      f=0;
            
      if (handles.FieldExpired ~=0)
      handles.VOI_parametersexp=assessfielddistribution( handles.FieldExpired,handles.Params.shimVoi);
      a=0;
      b=0;
      c=0;
      handles.mean_roiexp=strcat('Mean Abs =',num2str(handles.VOI_parametersexp.meanAbs));
      d=handles.mean_roi;
      handles.median_roiexp=strcat('Median =',num2str(handles.VOI_parametersexp.median));
      e=handles.median_roi;
      handles.std_roiexp=strcat('Std dev =',num2str(handles.VOI_parametersexp.std));
      f=handles.std_roi;
      end
      
             
   function Switch_plot(handles)
           
     switch handles.item_selected
            
        case 'Phase/Inspired'
              ImageField(handles.FieldInspired,handles.Fieldmaps,handles);
              Setparameters(handles.Field_mean,handles.Field_median,handles.Field_std,handles.meanins,handles.medianins,handles.stdins);
            
              if ~isempty(handles.Voi)
                  ImageFieldRoi(handles.FieldInspired,handles.Roi,handles.Voi,handles);
                  Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,handles.mean_roi,handles.median_roi,handles.std_roi);
              end
              
              if ~isempty(handles.Total_Voi)
                  ImageFieldRoi(handles.FieldInspired,handles.Roi,handles.Total_Voi(:,:,handles.sliceSelected),handles);
                  Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,handles.mean_roi,handles.median_roi,handles.std_roi);
              end
               
              if ~isempty(handles.PredictedFieldInspired)
                  ImageField(handles.PredictedFieldInspired,handles.Predicted,handles);  
                  Setparameters(handles.predictedmean,handles.predictedmedian,handles.predictedstd,handles.mean_pre,handles.median_pre,handles.std_pre);
                  Setcurrents(handles.predicted1,handles.predicted2,handles.predicted3,handles.predicted4,handles.predicted5,handles.predicted6,handles.predicted7,handles.predicted8,handles.Params.Inspired.currents);

              end
           
              if any(handles.sliceSelected == handles.sliceDeleted)==1  
                   Imageclear(handles.limits,handles.Roi,handles.colormap,handles);
              end
             
              if handles.FieldInspired.img(:,:,handles.sliceSelected)==0
                  if handles.n_region >=1
                    Imageclear(handles.limits,handles.Roi,handles.colormap,handles);
                  end
              end 
            
         case 'Phase/Expired'
            if (handles.FieldExpired ~=0)
               ImageField(handles.FieldExpired,handles.Fieldmaps,handles);
               Setparameters(handles.Field_mean,handles.Field_median,handles.Field_std,handles.meanexp,handles.medianexp,handles.stdexp);
                         
            if ~isempty(handles.Voi)
                  ImageFieldRoi(handles.FieldExpired,handles.Roi,handles.Voi,handles); 
                  Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,handles.mean_roiexp,handles.median_roiexp,handles.std_roiexp);
            end  
            
             if ~isempty(handles.Total_Voi)
                  ImageFieldRoi(handles.FieldExpired,handles.Roi,handles.Total_Voi(:,:,handles.sliceSelected),handles); 
                  Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,handles.mean_roiexp,handles.median_roiexp,handles.std_roiexp);
            end  
            
            if ~isempty(handles.PredictedFieldExpired)
                  ImageField(handles.PredictedFieldExpired,handles.Predicted,handles); 
                  Setparameters(handles.predictedmean,handles.predictedmedian,handles.predictedstd,handles.mean_preexp,handles.median_preexp,handles.std_preexp);
                  Setcurrents(handles.predicted1,handles.predicted2,handles.predicted3,handles.predicted4,handles.predicted5,handles.predicted6,handles.predicted7,handles.predicted8,handles.Params.Expired.currents);

            end     
          
             if any(handles.sliceSelected == handles.sliceDeleted)==1  
                   Imageclear(handles.limits,handles.Roi,handles.colormap,handles);
             end
           
             if handles.FieldExpired.img(:,:,handles.sliceSelected)==0
                 if handles.n_region >=1
                   Imageclear(handles.limits,handles.Roi,handles.colormap,handles)
                 end
             end 
             else
                 display ('There is no Expired field maps');
             end    
             
        case 'Mag/Inspired'
             ImageMag(handles.magInspired,handles.Fieldmaps,handles);
             Setparameters(handles.Field_mean,handles.Field_median,handles.Field_std,'Mean Abs: -','Median : -','Std dev :-');
             
            if ~isempty(handles.Voi)
                ImageMagRoi(handles.magInspired,handles.Roi,handles.Voi,handles);
                Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,'Mean Abs: -','Median : -','Std dev :-')
            end
                
             if ~isempty(handles.Total_Voi)
                ImageMagRoi(handles.magInspired,handles.Roi,handles.Total_Voi(:,:,handles.sliceSelected),handles);
                Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,'Mean Abs: -','Median : -','Std dev :-')
             end
            
            if ~isempty(handles.PredictedFieldInspired)            
                ImageField(handles.PredictedFieldInspired,handles.Predicted,handles);
                Setparameters(handles.predictedmean,handles.predictedmedian,handles.predictedstd,handles.mean_pre,handles.median_pre,handles.std_pre);
                Setcurrents(handles.predicted1,handles.predicted2,handles.predicted3,handles.predicted4,handles.predicted5,handles.predicted6,handles.predicted7,handles.predicted8,handles.Params.Inspired.currents);

            end
            
            if any(handles.sliceSelected == handles.sliceDeleted)==1  
                 Imageclear(handles.limits_mag,handles.Roi,gray,handles)
            end
             
            if handles.FieldInspired.img(:,:,handles.sliceSelected)==0
                if handles.n_region >=1
                  Imageclear(handles.limits_mag,handles.Roi,gray,handles)
                end
            end 
                 
        case 'Mag/Expired'
            if (handles.FieldExpired ~=0)
             ImageMag(handles.magExpired,handles.Fieldmaps,handles);
             Setparameters(handles.Field_mean,handles.Field_median,handles.Field_std,'Mean Abs: -','Median : -','Std dev :-');
             
           if ~isempty(handles.Voi)  
               ImageMagRoi(handles.magExpired,handles.Roi,handles.Voi,handles);
               Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,'Mean Abs: -','Median : -','Std dev :-');
           end
           
           if ~isempty(handles.Total_Voi)  
               ImageMagRoi(handles.magExpired,handles.Roi,handles.Total_Voi(:,:,handles.sliceSelected),handles);
               Setparameters(handles.VOI_mean,handles.VOI_median,handles.VOI_std,'Mean Abs: -','Median : -','Std dev :-');
           end
            
           if ~isempty(handles.PredictedFieldExpired)
               ImageField(handles.PredictedFieldExpired,handles.Predicted,handles); 
               Setparameters(handles.predictedmean,handles.predictedmedian,handles.predictedstd,handles.mean_preexp,handles.median_preexp,handles.std_preexp);
               Setcurrents(handles.predicted1,handles.predicted2,handles.predicted3,handles.predicted4,handles.predicted5,handles.predicted6,handles.predicted7,handles.predicted8,handles.Params.Expired.currents);

           end
           
           if any(handles.sliceSelected == handles.sliceDeleted)==1  
                Imageclear(handles.limits_mag,handles.Roi,gray,handles)
           end
             
           if handles.FieldInspired.img(:,:,handles.sliceSelected)==0
               if handles.n_region >=1
                 Imageclear(handles.limits_mag,handles.Roi,gray,handles)
               end
           end
           else
                 display ('There is no Expired field maps');
           end    
             
             
     end
     
      if ~isempty(handles.Voi)
            if nnz(handles.FieldInspired.img(:,:,handles.sliceSelected))~= 0
                  if any(handles.sliceSelected == handles.sliceDeleted)==0
                        for i=1:handles.n_region
                              rectangle(handles.Fieldmaps,'Position',handles.position{i},'EdgeColor','r','LineWidth',1);
                        end
                        for i=1:handles.n_region_exclude
                              rectangle(handles.Fieldmaps,'Position',handles.position_exclude{i},'EdgeColor','b','LineWidth',1);
                        end
                  end
            end             
      end
         
      
    if ~isempty(handles.PredictedFieldInspired) 
         if any(handles.sliceSelected == handles.sliceDeleted)==0
             if nnz(handles.FieldInspired.img(:,:,handles.sliceSelected))~= 0
                axes(handles.Predicted);
                if handles.n_region >= 1
                    for k = 1:length(handles.bound)
                    boundary = handles.bound{k};     
                    plot(handles.Predicted,boundary(:,2),boundary(:,1), 'Color','black', 'LineWidth', 1); 
                    end
                else
                    handles.sct_seg = handles.bound{handles.sliceSelected};
                    for k = 1:length(handles.sct_seg)                                   
                    boundary = handles.sct_seg{k};
                    plot(handles.Predicted,boundary(:,2),boundary(:,1), 'Color','black', 'LineWidth', 0.8);
                    end
                end
             end
         end
    end
       

             
% =========================================================================    
%Functions Set parameters on the GUI interface
% =========================================================================  

function Setparameters(location1,location2,location3,parameter1,parameter2,parameter3)
         set(location1,'string',parameter1);
         set(location2,'string',parameter2);
         set(location3,'string',parameter3);
         
function Setcurrents(loc1,loc2,loc3,loc4,loc5,loc6,loc7,loc8, current)
         currents=current*1000;
         set(loc1,'string',strcat('Ch1 :',num2str(currents(1),4),'mA'));
         set(loc2,'string',strcat('Ch2 :',num2str(currents(2),4),'mA'));
         set(loc3,'string',strcat('Ch3 :',num2str(currents(3),4),'mA'));
         set(loc4,'string',strcat('Ch4 :',num2str(currents(4),4),'mA'));
         set(loc5,'string',strcat('Ch5 :',num2str(currents(5),4),'mA'));
         set(loc6,'string',strcat('Ch6 :',num2str(currents(6),4),'mA'));
         set(loc7,'string',strcat('Ch7 :',num2str(currents(7),4),'mA'));
         set(loc8,'string',strcat('Ch8 :',num2str(currents(8),4),'mA'));
         
             
    
% =========================================================================    
%Function Create Fcn for the view selected button
% ========================================================================= 

function View_selected_CreateFcn(hObject, ~, ~)
    
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',{'Phase/Inspired';'Phase/Expired';'Mag/Inspired';'Mag/Expired'});
