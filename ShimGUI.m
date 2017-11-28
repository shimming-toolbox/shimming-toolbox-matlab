function varargout = ShimGUI(varargin)
% SHIMGUI MATLAB code for ShimGUI.fig

%      This function allows to check insp/exp mag and phase data and define
%      ROIs for optimizing shimming.
%
%      Usage:

%      You need four training maps, in DICOM format, located in:
%        'ins/mag/'    Magnitudes Inspired
%        'ins/phase/'  Phases Inspired
%        'exp/mag/'    Magnitudes Expired
%        'exp/phase/'  Phases Expired
%      
%       You also need a configuration file which defines the function: 
%       shimparameters(). Example: see: shimparameters.m

%
% Last Modified by GUIDE v2.5 28-Nov-2017 13:56:55


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
handles.item_selected = 'Phase/Inspired';
handles.lim=200;
handles.limits=[-handles.lim handles.lim];
handles.lim_mag=0.5;
handles.limits_mag=[-handles.lim_mag handles.lim_mag];
handles.n_region=0;
handles.n_region_exclude=0;
handles.Voi=[];
handles.Voi_exclude=[];
handles.position=cell(1,1);
handles.position_exclude=cell(1,1);
handles.PredictedFieldInspired= [];
handles.PredictedFieldExpired= [];
handles.sliceDeleted=[];
handles.colormap='parula';
handles.bound=cell(1,1);
display ('Load your training maps');

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
    handles.Shims = ShimOptAcdc(handles.Params) ;
end

%Reverse polarization for the shims =======================================

handles.Shims.img = -handles.Shims.img ;

% =========================================================================
% Prepare the calibration maps
% =========================================================================

ImgArray_Inspired = cell(1 , 2);
ImgArray_Inspired{1,1} = MaRdI( [ handles.Params.dataLoadDir 'ins/mag/' ] ) ;
ImgArray_Inspired{1,2} = MaRdI( [ handles.Params.dataLoadDir 'ins/phase/' ] ) ;

handles.FieldInspired = ShimOpt.mapfield( ImgArray_Inspired, handles.Params ) ;
display('Loading Inspired Fieldmap--------------> Done');

handles.magInspired = ImgArray_Inspired{1,1};
display('Loading Inspired Magnitudes --------------> Done');


ImgArray_Expired = cell( 1, 2 ) ;
ImgArray_Expired{1,1} = MaRdI( [ handles.Params.dataLoadDir 'exp/mag/' ] ) ;
ImgArray_Expired{1,2} = MaRdI( [ handles.Params.dataLoadDir 'exp/phase/' ] ) ;

handles.FieldExpired = ShimOpt.mapfield( ImgArray_Expired, handles.Params ) ;
display('Loading Expired Fieldmap  --------------> Done');

handles.magExpired = ImgArray_Expired{1,1};
display('Loading Expired Magnitudes --------------> Done');


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
handles.Shims.interpolatetoimggrid( handles.FieldInspired );

%Display Inspired Fieldmap ------------------------------------------------

ImageField(handles.FieldInspired,handles.Fieldmaps,handles);
handles.magInspired.img=handles.magInspired.img-0.5;
handles.magExpired.img=handles.magExpired.img-0.5;


guidata(hObject, handles) ;



% --- Executes on button press in customROI.
function customROI_Callback(hObject,~, handles)

     if any(handles.sliceSelected == handles.sliceDeleted)>=1
                display('Change the slice selected, this one was removed from the Voi')
     else
            handles.n_region = handles.n_region+1;
                
            
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
% DEFINE SHIM VOI 
% =========================================================================

            handles.Params.shimVoi = handles.Shims.getvaliditymask( handles.Params, handles.FieldInspired, handles.FieldExpired ) ;

% Adjust shim VOI based on the rectangular selection on the image.

            if ~isempty(handles.Voi_exclude);
                handles.Voi = handles.Voi.*handles.Voi_exclude;
            end

            for j=1:handles.dim(3);
                handles.Params.shimVoi(:,:,j)=handles.Params.shimVoi(:,:,j).*handles.Voi;
            end
   

             handles.bound=bwboundaries(handles.Voi,'noholes');
             
       switch handles.item_selected
            
          case 'Phase/Inspired'
              ImageFieldRoi(handles.FieldInspired,handles.Roi,handles.Voi,handles);  
          case 'Phase/Expired'
              ImageFieldRoi(handles.FieldExpired,handles.Roi,handles.Voi,handles);              
          case 'Mag/Inspired'             
              ImageMagRoi(handles.magInspired,handles.Roi,handles.Voi,handles);
          case 'Mag/Expired'
              ImageMagRoi(handles.magExpired,handles.Roi,handles.Voi,handles);
       end
     end
guidata(hObject, handles) ;



% --- Executes on slider movement.
function sliceSelector_Callback(hObject,~, handles)

dim = size(handles.FieldInspired.img);
maxslice = dim(3);     %Nombre de slice dans le plan sagittal
handles.sliceSelected = round(get(hObject,'Value')*(maxslice - 1) + 1);

Switch_plot(handles);

set(handles.commandLine,'string',num2str(round(get(hObject,'Value')*(maxslice - 1) + 1)));
            
guidata(hObject, handles) ;



%View Selected Callback
function View_selected_Callback(hObject, ~, handles)
    
items = get(hObject,'String');
index_selected = get(hObject,'Value');
handles.item_selected = items{index_selected};

Switch_plot(handles);
            
guidata(hObject, handles);



% --- Executes on button press in Prediction.
function Prediction_Callback(hObject, ~, handles)

handles.Params.isSolvingAugmentedSystem    = true ;
handles.Params.isPenalizingFieldDifference = true;
handles.Params.regularizationParameter     = 0 ;
    
handles.Shims.setoriginalfield( handles.FieldInspired ) ;
handles.Shims.setshimvolumeofinterest( handles.Params.shimVoi) ;

[handles.Params.Inspired.currents, handles.Params.Expired.currents] = handles.Shims.optimizeshimcurrents( handles.Params, handles.FieldExpired ) ;
    
%Shimming predictions------------------------------------------------------

handles.Shims.Model.currents   =  handles.Params.Inspired.currents ;
handles.PredictedFieldInspired =  handles.Shims.predictshimmedfield( ) ;
                
handles.Shims.setoriginalfield( handles.FieldExpired ) ;
handles.Shims.Model.currents = handles.Params.Expired.currents ;
handles.PredictedFieldExpired =  handles.Shims.predictshimmedfield( ) ;
    
    
%Mask for regions without any signal in the Fieldmaps----------------------

mask = handles.FieldInspired.img;
mask(mask ~= 0) = 1;
    
handles.PredictedFieldInspired.img = handles.PredictedFieldInspired.img .* mask;
    
mask2 = handles.FieldExpired.img;
mask2(mask2 ~= 0) = 1;
handles.PredictedFieldExpired.img = handles.PredictedFieldExpired.img .* mask2;

            
                
switch handles.item_selected
                  
   case 'Phase/Inspired'         
         ImageField(handles.PredictedFieldInspired,handles.Predicted,handles); 
            
              
   case 'Phase/Expired'
         ImageField(handles.PredictedFieldExpired,handles.Predicted,handles);                                
                
   case 'Mag/Inspired' 

         ImageField(handles.PredictedFieldInspired,handles.Predicted,handles);
                               
   case 'Mag/Expired'             
         ImageField(handles.PredictedFieldExpired,handles.Predicted,handles);  
                                              
end

     if any(handles.sliceSelected == handles.sliceDeleted)==0
         axes(handles.Predicted);
            hold on
         for k = 1:length(handles.bound)
            boundary = handles.bound{k};     
            plot(handles.Predicted,boundary(:,2),boundary(:,1), 'Color','black', 'LineWidth', 2);           
         end
     end
guidata(hObject, handles);


% --- Executes on button press in delete_region.
function delete_region_Callback(hObject, ~, handles)

    
handles.Params.shimVoi(:,:,handles.sliceSelected)=0;
handles.sliceDeleted(end+1)=handles.sliceSelected;


switch handles.item_selected
    
    case 'Phase/Inspired'
        ImageField(handles.FieldInspired,handles.Fieldmaps,handles);
        if handles.n_region >=1
            Imageclear(handles.limits,handles.Roi,handles.colormap,handles);
        end           
           
    case 'Phase/Expired'
         ImageField(handles.FieldExpired,handles.Fieldmaps,handles);
         if handles.n_region >=1
             Imageclear(handles.limits,handles.Roi,handles.colormap,handles);
         end           
    case 'Mag/Inspired'
         ImageMag(handles.magInspired,handles.Fieldmaps,handles);  
         if handles.n_region >=1
            Imageclear(handles.limits_mag,handles.Roi,gray,handles);
         end 
                 
    case 'Mag/Expired'
         ImageMag(handles.magExpired,handles.Fieldmaps,handles);
         if handles.n_region >=1
            Imageclear(handles.limits_mag,handles.Roi,gray,handles);
         end
end

if handles.PredictedFieldExpired ~=0
   Imageclear(handles.limits,handles.Predicted,handles.colormap,handles);
end
guidata(hObject,handles);



% --- Executes on button press in Clear_Roi.
function Clear_Roi_Callback(hObject, ~, handles)

    
    switch handles.item_selected
            
         case 'Phase/Inspired'
             ImageField(handles.FieldInspired,handles.Fieldmaps,handles); 
             if handles.n_region >= 1
             Imageclear(handles.limits,handles.Roi,handles.colormap,handles);
             end
            
         case 'Phase/Expired'
             ImageField(handles.FieldExpired,handles.Fieldmaps,handles);
             if handles.n_region >= 1
             Imageclear(handles.limits,handles.Roi,handles.colormap,handles);
             end
                         
         case 'Mag/Inspired'
             ImageMag(handles.magInspired,handles.Fieldmaps,handles) ;
             if handles.n_region >= 1
             Imageclear(handles.limits_mag,handles.Roi,gray,handles);
             end
                 
         case 'Mag/Expired'
             ImageMag(handles.magExpired,handles.Fieldmaps,handles);
             if handles.n_region >= 1
             Imageclear(handles.limits_mag,handles.Roi,gray,handles);
             end
             
    end
     
     if handles.PredictedFieldExpired ~=0
        Imageclear(handles.limits,handles.Predicted,handles.colormap,handles);
     end
     
handles.Voi = [];
handles.Voi_exclude=[];
handles.position=cell(1,1);
handles.position_exclude=cell(1,1);
handles.n_region=0;
handles.n_region_exclude=0;
handles.Params.shimVoi(:,:,:)=0;
handles.PredictedFieldInspired=[];
handles.PredictedFieldExpired=[];
handles.sliceDeleted=[];
handles.bound=[];

 guidata(hObject, handles);
 
 % --- Executes on button press in Exclude_From_ROI.
function Exclude_From_ROI_Callback(hObject, ~, handles)
% hObject    handle to Exclude_From_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
     if ~isempty(handles.Voi)
            handles.n_region_exclude = handles.n_region_exclude+1;
            
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
% DEFINE SHIM VOI 
% =========================================================================

    handles.Voi = handles.Voi.*handles.Voi_exclude;
    handles.Params.shimVoi = handles.Shims.getvaliditymask( handles.Params, handles.FieldInspired, handles.FieldExpired ) ;
    handles.bound=bwboundaries(handles.Voi);

% Adjust shim VOI based on the rectangular selection on the image.

            for j=1:handles.dim(3);
                handles.Params.shimVoi(:,:,j)=handles.Params.shimVoi(:,:,j).*handles.Voi;
            end
    
        switch handles.item_selected
            
                case 'Phase/Inspired'
                    ImageFieldRoi(handles.FieldInspired,handles.Roi,handles.Voi,handles);  
                case 'Phase/Expired'
                    ImageFieldRoi(handles.FieldExpired,handles.Roi,handles.Voi,handles);              
                case 'Mag/Inspired'             
                 ImageMagRoi(handles.magInspired,handles.Roi,handles.Voi,handles);
                case 'Mag/Expired'
                    ImageMagRoi(handles.magExpired,handles.Roi,handles.Voi,handles);
         end
     else
         display ('Warning : Create a Voi before using the button "Exclude"')
     end   
    
    guidata (hObject,handles);




function Phase_contrast_selector_Callback(hObject, ~, handles)

handles.min_contrast=max(handles.FieldInspired.img(:));
handles.lim = round(get(hObject,'Value')*(handles.min_contrast - 1) + 1);
handles.limits=[-handles.lim handles.lim];

Switch_plot(handles);

guidata (hObject,handles);


function Mag_contrast_selector_Callback(hObject, ~, handles)

handles.lim_mag = get(hObject,'Value');
if handles.lim_mag <=0.01
    handles.lim_mag=0.011;
end

handles.limits_mag=[-handles.lim_mag handles.lim_mag];

Switch_plot(handles);

guidata (hObject,handles);

%Functions Images ---------------------------------------------------------

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
             
    
   function Switch_plot(handles)
           
     switch handles.item_selected
            
        case 'Phase/Inspired'
              ImageField(handles.FieldInspired,handles.Fieldmaps,handles);
            
              if ~isempty(handles.Voi)
                  ImageFieldRoi(handles.FieldInspired,handles.Roi,handles.Voi,handles); 
              end
              
              if ~isempty(handles.PredictedFieldInspired)
                  ImageField(handles.PredictedFieldInspired,handles.Predicted,handles);      
              end
           
              if any(handles.sliceSelected == handles.sliceDeleted)==1  
                   Imageclear(handles.limits,handles.Roi,handles.colormap,handles)
              end
             
              if handles.FieldInspired.img(:,:,handles.sliceSelected)==0
                  if handles.n_region >=1
                    Imageclear(handles.limits,handles.Roi,handles.colormap,handles)
                  end
              end 
            
         case 'Phase/Expired'
               ImageField(handles.FieldExpired,handles.Fieldmaps,handles);
                         
            if ~isempty(handles.Voi)
                  ImageFieldRoi(handles.FieldExpired,handles.Roi,handles.Voi,handles); 
            end  
          
          
            if ~isempty(handles.PredictedFieldExpired)
                  ImageField(handles.PredictedFieldExpired,handles.Predicted,handles);  
            end     
          
             if any(handles.sliceSelected == handles.sliceDeleted)==1  
                   Imageclear(handles.limits,handles.Roi,handles.colormap,handles)
             end
           
             if handles.FieldExpired.img(:,:,handles.sliceSelected)==0
                 if handles.n_region >=1
                   Imageclear(handles.limits,handles.Roi,handles.colormap,handles)
                 end
             end 
             
        case 'Mag/Inspired'
             ImageMag(handles.magInspired,handles.Fieldmaps,handles);
             
            if ~isempty(handles.Voi)
                ImageMagRoi(handles.magInspired,handles.Roi,handles.Voi,handles);
            end
            
            if ~isempty(handles.PredictedFieldInspired)            
                ImageField(handles.PredictedFieldInspired,handles.Predicted,handles);                         
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
             ImageMag(handles.magExpired,handles.Fieldmaps,handles);
             
           if ~isempty(handles.Voi)  
               ImageMagRoi(handles.magExpired,handles.Roi,handles.Voi,handles);
           end
            
           if ~isempty(handles.PredictedFieldExpired)
               ImageField(handles.PredictedFieldExpired,handles.Predicted,handles);             
           end
           
           if any(handles.sliceSelected == handles.sliceDeleted)==1  
                Imageclear(handles.limits_mag,handles.Roi,gray,handles)
           end
             
           if handles.FieldInspired.img(:,:,handles.sliceSelected)==0
               if handles.n_region >=1
                 Imageclear(handles.limits_mag,handles.Roi,gray,handles)
               end
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
                    hold on
                    for k = 1:length(handles.bound)
                       boundary = handles.bound{k};     
                       plot(handles.Predicted,boundary(:,2),boundary(:,1), 'Color','black', 'LineWidth', 2);           
                    end
              end
           end
       end

%Functions CreateFcn--------------------------------------------------------
function View_selected_CreateFcn(hObject, ~, ~)
    
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',{'Phase/Inspired';'Phase/Expired';'Mag/Inspired';'Mag/Expired'});
