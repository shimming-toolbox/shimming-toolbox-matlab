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
% Last Modified by GUIDE v2.5 26-Jan-2018 13:48:55


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
handles.itemSelected = 'Phase/Inspired';           %View selected
handles.lim=200;                                    %Constrast in Images
handles.limits=[-handles.lim handles.lim];
handles.limMag=0.5;
handles.limitsMag=[-handles.limMag handles.limMag]; 
handles.colormap='parula'; 
handles.nRegion=0;                                 %Voi number of regions 
handles.nRegionexclude=0;                           
handles.voi=[];
handles.totalVoi=[]; %Voi mask     
handles.voiExclude=[];                                
handles.position=cell(1,1);                         %Voi regions positions
handles.positionExclude=cell(1,1);
handles.bound=cell(1,1);
handles.sctVoi=[];
handles.predictedfieldInspired= [];                 %Predicted fields
handles.predictedfieldExpired= [];
handles.fieldExpired=0;
handles.sliceDeleted=[];                            %Slice excluded from Voi                         
handles.sortedData=' ';                            %Directory with sorted data
handles.calibrationValues=[];
handles.valuestoSend=[];
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

imgArrayInspired = cell(1 , 2);

%Loading of Magnitude Maps ================================================

display('Select folder with magnitude maps in inspired state')
if handles.sortedData==' '
    %pathtomagInspired=uigetdir();
    imgArrayInspired{1,1} = MaRdI( [ '/Users/ancha_admin/data/sorted_data/08_gre_field_mapping_shim0_ins/ORIGINAL\PRIMARY\M\ND/1' ] ) ;  % ATTENTION HARDCODE FOR TEST
else
handles.pathtomagInspired=uigetdir(handles.sortedData);
imgArrayInspired{1,1} = MaRdI( handles.pathtomagInspired ) ;
end

%Loading of Phase Maps ====================================================

display('Select folder with phase maps in inspired state')
if handles.sortedData==' '
    %pathtophaseInspired=uigetdir();
    imgArrayInspired{1,2} = MaRdI( [ '/Users/ancha_admin/data/sorted_data/10_gre_field_mapping_shim0_ins/ORIGINAL\PRIMARY\P\ND/2' ] ) ;   % ATTENTION HARDCODE FOR TEST
else
handles.pathtophaseInspired = uigetdir(handles.sortedData);
imgArrayInspired{1,2} = MaRdI(handles.pathtophaseInspired) ;
end


handles.fieldInspired = ShimOpt.mapfield( imgArrayInspired, handles.Params ) ;
display('Loading Inspired Fieldmap--------------> Done');

handles.magInspired = imgArrayInspired{1,1};
display('Loading Inspired Magnitudes --------------> Done');

%Definition of new variables from the fieldmaps----------------------------

handles.Params.scaling = [min(handles.fieldInspired.img(:)) max(handles.fieldInspired.img(:))] ;
handles.dim = size(handles.fieldInspired.img);
handles.sliceSelected = round(handles.dim(3)*0.5);
handles.clear= zeros(size(handles.fieldInspired.img));

%Set the slice selected and the contrasts on the images--------------------

set(handles.sliceSelector,'value',0.5);
set(handles.magcontrastSelector,'value',0.5);
set(handles.phasecontrastSelector,'value',0.3);
set(handles.commandLine,'string',num2str(handles.sliceSelected));

%Interpolation to the image grid-------------------------------------------

handles.fieldInspired.Hdr.Private_0019_1014 = [0 0 0] ;
handles.Shims.Opt.interpolatetoimggrid( handles.fieldInspired );

%Display Inspired Fieldmap ------------------------------------------------

imagefield(handles.fieldInspired,handles.Fieldmaps,handles);
handles.magInspired.img=handles.magInspired.img-0.5;

fieldParametersins=assessfielddistribution( handles.fieldInspired );
handles.meanIns=strcat('Mean Abs=',num2str(fieldParametersins.meanAbs));
handles.medianIns=strcat('Median =',num2str(fieldParametersins.median));
handles.stdIns=strcat('Std dev =',num2str(fieldParametersins.std));
setparameters(handles.fieldMean,handles.fieldMedian,handles.fieldStd,handles.meanIns,handles.medianIns,handles.stdIns);

%Display the Feedback from the background script---------------------------

File=fopen('background');
a=fscanf(File,'%c');
disp(a);
fclose(File);

guidata(hObject, handles) ;

% =========================================================================
% Loading of the Expired maps
% =========================================================================

function Load_expired_Callback(hObject, ~, handles)
imgArrayExpired = cell( 1, 2 ) ;

%Loading of Magnitude Maps ================================================
display('Select folder with magnitude maps in expired state')
handles.pathtomagExpired=uigetdir(handles.sortedData);
imgArrayExpired{1,1} = MaRdI( handles.pathtomagExpired ) ;

%Loading of Phase Maps ====================================================
display('Select folder with phase maps in expired state')
handles.pathtophaseExpired=uigetdir(handles.sortedData);
imgArrayExpired{1,2} = MaRdI( handles.pathtophaseExpired ) ;

handles.fieldExpired = ShimOpt.mapfield( imgArrayExpired, handles.Params ) ;
display('Loading Expired Fieldmap  --------------> Done');

handles.magExpired = imgArrayExpired{1,1};
display('Loading Expired Magnitudes --------------> Done');

%Declaration of the field distribution parameters==========================
fieldParametersexp=assessfielddistribution( handles.fieldExpired );
handles.meanExp=strcat('Mean Abs=',num2str(fieldParametersexp.meanAbs));
handles.medianExp=strcat('Median =',num2str(fieldParametersexp.median));
handles.stdExp=strcat('Std dev =',num2str(fieldParametersexp.std));

handles.magExpired.img=handles.magExpired.img-0.5;

%Feedback from the background script=======================================

feedback=fopen('background');
a=fscanf(feedback,'%c');
disp(a);
fclose(feedback);
guidata(hObject, handles) ;



% --- Executes on button press in customVOI.
function customVOI_Callback(hObject,~, handles)

     if any(handles.sliceSelected == handles.sliceDeleted)>=1
      display('Change the slice selected, this one was removed from the Voi')
     else
            handles.nRegion = handles.nRegion+1;
% =========================================================================
% Selection of the Regions inside the VOI 
% =========================================================================                         
         if  isempty(handles.voi)
                handles.rect = imrect(handles.Fieldmaps);
                setColor(handles.rect,'r');  
                handles.voi =handles.rect.createMask;
                handles.position{1}=handles.rect.getPosition;        
         else        
                handles.rectn = imrect(handles.Fieldmaps);
                setColor(handles.rectn,'r');   
                handles.mask = handles.rectn.createMask;
                handles.voi=or(handles.voi,handles.mask);
                handles.position{end+1}=handles.rectn.getPosition;
         end
         
% =========================================================================
% DEFINE Validity mask for shim VOI 
% =========================================================================

            if (handles.fieldExpired ~=0)
            handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.fieldInspired, handles.fieldExpired ) ;
            else
            handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.fieldInspired) ;
            end
% =========================================================================
% Adjust shim VOI based on the rectangular selection on the image
% =========================================================================
            if ~isempty(handles.voiExclude);
                handles.voi = handles.voi.*handles.voiExclude;
            end          
            
            for j=1:handles.dim(3);
                handles.Params.shimVoi(:,:,j)=handles.Params.shimVoi(:,:,j).*handles.voi;
            end
            
% =========================================================================
% Field distribution parameters inside the VOI
% =========================================================================           
            handles.voiparametersIns=assessfielddistribution( handles.fieldInspired,handles.Params.shimVoi);
            handles.meanRoi=strcat('Mean Abs=',num2str(handles.voiparametersIns.meanAbs));
            handles.medianRoi=strcat('Median = ',num2str(handles.voiparametersIns.median));
            handles.stdRoi=strcat('Std dev = ',num2str(handles.voiparametersIns.std));
            
            if (handles.fieldExpired ~=0)
            handles.voiparametersExp=assessfielddistribution( handles.fieldExpired,handles.Params.shimVoi);
            handles.meanRoiexp=strcat('Mean Abs=',num2str(handles.voiparametersExp.meanAbs));
            handles.medianRoiexp=strcat('Median =',num2str(handles.voiparametersExp.median));
            handles.stdRoiexp=strcat('Std dev =',num2str(handles.voiparametersExp.std));
            end
   
             handles.bound=bwboundaries(handles.voi);
             
% =========================================================================
% Plot Images and field parameters inside the VOI 
% =========================================================================             
       switch handles.itemSelected
            
          case 'Phase/Inspired'
              imagefieldroi(handles.fieldInspired,handles.Roi,handles.voi,handles); 
              setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,handles.meanRoi,handles.medianRoi,handles.stdRoi)
          case 'Phase/Expired'
              imagefieldroi(handles.fieldExpired,handles.Roi,handles.voi,handles);  
              setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,handles.meanRoiexp,handles.medianRoiexp,handles.stdRoiexp)
          case 'Mag/Inspired'             
              imagemagroi(handles.magInspired,handles.Roi,handles.voi,handles);
              setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,'Mean : -','Median : -','Std dev :-')
          case 'Mag/Expired'
              imagemagroi(handles.magExpired,handles.Roi,handles.voi,handles);
              setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,'Mean : -','Median : -','Std dev :-')
       end
     end
     
% =========================================================================
% Feedback from the backgroud script
% =========================================================================     
feedback=fopen('background');
a=fscanf(feedback,'%c');
disp(a);
fclose(feedback);
guidata(hObject, handles) ;



function Use_sct_Callback(hObject, eventdata, handles) 
    
            dicm2nii(handles.pathtomagInspired,handles.Params.matlabPath,0);        
            
            [~,~] = unix(handles.Params.command);
            handles.sctMaskins=load_untouch_nii('/gre_field_mapping_shim0_ins_seg.nii');
            handles.sctCSFMaskins=load_untouch_nii('/gre_field_mapping_shim0_ins_CSF_seg.nii');
            
            handles.sctMaskIns = double(handles.sctMaskins.img);
            handles.sctCSFMaskIns=double(handles.sctCSFMaskins.img);
            
            
            for i=1:handles.dim(3)
            handles.sctVoiIns(:,:,i)=rot90(handles.sctMaskIns(:,:,i));
            handles.sctCSFVoiIns(:,:,i)=rot90(handles.sctCSFMaskIns(:,:,i));
            end
                  
           
            
            if (handles.fieldExpired ~=0)
            handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.fieldInspired, handles.fieldExpired ) ;   
            dicm2nii(handles.pathtomagExpired,handles.Params.matlabPath,0); 
            [~,~] = unix(handles.Params.command2);
            
            handles.sctMaskexp=load_untouch_nii('/gre_field_mapping_shim0_exp_seg.nii');
            handles.sctCSFMaskexp=load_untouch_nii('/gre_field_mapping_shim0_exp_CSF_seg.nii');
            
            handles.sctMaskExp = double(handles.sctMaskexp.img);
            handles.sctCSFMaskExp=double(handles.sctCSFMaskexp.img);
            
            
            for i=1:handles.dim(3)
            handles.sctVoiExp(:,:,i)=rot90(handles.sctMaskExp(:,:,i));
            handles.sctCSFVoiExp(:,:,i)=rot90(handles.sctCSFMaskExp(:,:,i));
            end
            
            handles.totalVoi=handles.sctVoiIns+handles.sctCSFVoiIns+handles.sctVoiExp+handles.sctCSFVoiExp;
            
            else
            handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.fieldInspired) ;
            handles.totalVoi=handles.sctVoiIns+handles.sctCSFVoiIns;
            end     
            
            for j=1:handles.dim(3);
                handles.Params.shimVoi(:,:,j)=handles.Params.shimVoi(:,:,j).*handles.totalVoi(:,:,j);    
                handles.bound{j}=bwboundaries(handles.totalVoi(:,:,j));        
            end
            
                               
           [handles.meanRoi,handles.medianRoi,handles.stdRoi,handles.meanRoiexp,handles.medianRoiexp,handles.stdRoiexp]=calculparameters(handles);

            switchplot(handles);         
            
            
guidata(hObject, handles) ;

%Executes on slice selector movement=======================================

function sliceSelector_Callback(hObject,~, handles)

dim = size(handles.fieldInspired.img);

%Slice number in the sagittal plan ========================================
maxSlice = dim(3);                   

%Define the slice selected from the slider position
handles.sliceSelected = round(get(hObject,'Value')*(maxSlice - 1) + 1);

switchplot(handles);

%Display the slice selected on the Shim GUI interface======================

set(handles.commandLine,'string',num2str(round(get(hObject,'Value')*(maxSlice - 1) + 1)));
            
guidata(hObject, handles) ;



%View Selected Callback====================================================

function View_selected_Callback(hObject, ~, handles)
    
items = get(hObject,'String');
indexSelected = get(hObject,'Value');
handles.itemSelected = items{indexSelected};

switchplot(handles);
            
guidata(hObject, handles);



%Generate Predicted maps callback  ========================================

function Prediction_Callback(hObject, ~, handles)
    
handles.Shims.Opt.setoriginalfield( handles.fieldInspired ) ;
handles.Shims.Opt.setshimvolumeofinterest( handles.Params.shimVoi) ;

%Shim currents optimization ===============================================

if (handles.fieldExpired ~=0)
handles.Params.isSolvingAugmentedSystem    = true ;
handles.Params.isPenalizingFieldDifference = true;
handles.Params.regularizationParameter     = 0 ;
[handles.Params.Inspired.currents, handles.Params.Expired.currents] = handles.Shims.Opt.optimizeshimcurrents( handles.Params, handles.fieldExpired ) ;
else
[handles.Params.Inspired.currents] = handles.Shims.Opt.optimizeshimcurrents(handles.Params) ;
end
    
handles.Shims.Opt.Model.currents   =  handles.Params.Inspired.currents ;

%==========================================================================
%Predicted Inspired field calculation
%==========================================================================

handles.predictedfieldInspired =  handles.Shims.Opt.predictshimmedfield( ) ;   

%Mask for regions without any signal in the inspired Fieldmaps=============

mask = handles.fieldInspired.img;
mask(mask ~= 0) = 1;


handles.predictedfieldInspired.img = handles.predictedfieldInspired.img .* mask;

handles.predictedParametersins=assessfielddistribution(handles.predictedfieldInspired,handles.Params.shimVoi);
handles.meanPre=strcat('Mean Abs =',num2str(handles.predictedParametersins.meanAbs));
handles.medianPre=strcat('Median =',num2str(handles.predictedParametersins.median));
handles.stdPre=strcat('Std dev =',num2str(handles.predictedParametersins.std));


%==========================================================================
%Predicted Expired field calculation
%==========================================================================
if (handles.fieldExpired ~=0)
handles.Shims.Opt.setoriginalfield( handles.fieldExpired ) ;
handles.Shims.Opt.Model.currents = handles.Params.Expired.currents ;
handles.predictedfieldExpired =  handles.Shims.Opt.predictshimmedfield( ) ;

%Mask for regions without any signal in the expired Fieldmaps==============
mask2 = handles.fieldExpired.img;
mask2(mask2 ~= 0) = 1;
handles.predictedfieldExpired.img = handles.predictedfieldExpired.img .* mask2;

handles.predictedParametersexp=assessfielddistribution(handles.predictedfieldExpired,handles.Params.shimVoi);
handles.meanPreexp=strcat('Mean Abs=',num2str(handles.predictedParametersexp.meanAbs));
handles.medianPreexp=strcat('Median =',num2str(handles.predictedParametersexp.median));
handles.stdPreexp=strcat('Std dev =',num2str(handles.predictedParametersexp.std));

end  

%==========================================================================
%Plot Predicted field distribution parameters and Images
%==========================================================================                
switch handles.itemSelected
                  
   case 'Phase/Inspired'         
         imagefield(handles.predictedfieldInspired,handles.Predicted,handles); 
         setparameters(handles.predictedMean,handles.predictedMedian,handles.predictedStd,handles.meanPre,handles.medianPre,handles.stdPre)
         setcurrents(handles.predicted1,handles.predicted2,handles.predicted3,handles.predicted4,handles.predicted5,handles.predicted6,handles.predicted7,handles.predicted8,handles.Params.Inspired.currents);

              
   case 'Phase/Expired'
         imagefield(handles.predictedfieldExpired,handles.Predicted,handles);     
         setparameters(handles.predictedMean,handles.predictedMedian,handles.predictedStd,handles.meanPreexp,handles.medianPreexp,handles.stdPreexp)
         setcurrents(handles.predicted1,handles.predicted2,handles.predicted3,handles.predicted4,handles.predicted5,handles.predicted6,handles.predicted7,handles.predicted8,handles.Params.Expired.currents);

                
   case 'Mag/Inspired' 
         imagefield(handles.predictedfieldInspired,handles.Predicted,handles);
         setparameters(handles.predictedMean,handles.predictedMedian,handles.predictedStd,handles.meanPre,handles.medianPre,handles.stdPre)
         setcurrents(handles.predicted1,handles.predicted2,handles.predicted3,handles.predicted4,handles.predicted5,handles.predicted6,handles.predicted7,handles.predicted8,handles.Params.Inspired.currents);

                               
   case 'Mag/Expired'             
         imagefield(handles.predictedfieldExpired,handles.Predicted,handles);  
         setparameters(handles.predictedMean,handles.predictedMedian,handles.predictedStd,handles.meanPreexp,handles.medianPreexp,handles.stdPreexp)
         setcurrents(handles.predicted1,handles.predicted2,handles.predicted3,handles.predicted4,handles.predicted5,handles.predicted6,handles.predicted7,handles.predicted8,handles.Params.Inspired.currents);
                                     
end


     if any(handles.sliceSelected == handles.sliceDeleted)==0
         axes(handles.Predicted);
            hold on
            if handles.nRegion >= 1
             for k = 1:length(handles.bound)
               boundary = handles.bound{k};     
               plot(handles.Predicted,boundary(:,2),boundary(:,1), 'Color','black', 'LineWidth', 1); 
             end
            else
                handles.sctSeg = handles.bound{handles.sliceSelected};
              for k = 1:length(handles.sctSeg)                                   
               boundary = handles.sctSeg{k};     
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

      [handles.meanRoi,handles.medianRoi,handles.stdRoi,handles.meanRoiexp,handles.medianRoiexp,handles.stdRoiexp]=calculparameters(handles);


%Plot field distribution parameters and images-----------------------------

switch handles.itemSelected
    
    case 'Phase/Inspired'
        imagefield(handles.fieldInspired,handles.Fieldmaps,handles);
        setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,handles.meanRoi,handles.medianRoi,handles.stdRoi);
        if handles.nRegion >=1
            imageclear(handles.limits,handles.Roi,handles.colormap,handles);
        end           
           
    case 'Phase/Expired'
         imagefield(handles.fieldExpired,handles.Fieldmaps,handles);
         setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,handles.meanRoiexp,handles.medianRoiexp,handles.stdRoiexp);
         if handles.nRegion >=1
             imageclear(handles.limits,handles.Roi,handles.colormap,handles);
         end           
    case 'Mag/Inspired'
         imagemag(handles.magInspired,handles.Fieldmaps,handles);
         setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,'Mean Abs: -','Median : -','Std dev :-');
         if handles.nRegion >=1
            imageclear(handles.limitsMag,handles.Roi,gray,handles);
         end 
                 
    case 'Mag/Expired'
         imagemag(handles.magExpired,handles.Fieldmaps,handles);
         setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,'Mean Abs: -','Median : -','Std dev :-');
         if handles.nRegion >=1
            imageclear(handles.limitsMag,handles.Roi,gray,handles);
         end
end

if handles.predictedfieldExpired ~=0
   imageclear(handles.limits,handles.Predicted,handles.colormap,handles);
end
guidata(hObject,handles);



% Clear VOI callback=======================================================

function Clear_Roi_Callback(hObject, ~, handles)

    switch handles.itemSelected
            
         case 'Phase/Inspired'
             imagefield(handles.fieldInspired,handles.Fieldmaps,handles); 
             imageclear(handles.limits,handles.Roi,handles.colormap,handles);

            
         case 'Phase/Expired'
             imagefield(handles.fieldExpired,handles.Fieldmaps,handles);
             imageclear(handles.limits,handles.Roi,handles.colormap,handles);
             
                         
         case 'Mag/Inspired'
             imagemag(handles.magInspired,handles.Fieldmaps,handles) ;
             imageclear(handles.limitsMag,handles.Roi,gray,handles);
             
                 
         case 'Mag/Expired'
             imagemag(handles.magExpired,handles.Fieldmaps,handles);
             imageclear(handles.limitsMag,handles.Roi,gray,handles);
             
             
    end
     
     if handles.predictedfieldInspired ~=0
        imageclear(handles.limits,handles.Predicted,handles.colormap,handles);
     end
     
     
%==========================================================================
% reset_current VOI parameters
%==========================================================================

handles.voi = [];
handles.totalVoi=[];
handles.sctVoi=[];
handles.voiExclude=[];
handles.voiparametersIns=0;
handles.voiparametersExp=0;
handles.position=cell(1,1);
handles.positionExclude=cell(1,1);
handles.nRegion=0;
handles.nRegionexclude=0;
handles.Params.shimVoi(:,:,:)=0;
handles.predictedfieldInspired=[];
handles.predictedfieldExpired=[];
handles.sliceDeleted=[];
handles.bound=[];
setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,'Mean Abs: -','Median : -','Std dev :-');
setparameters(handles.predictedMean,handles.predictedMedian,handles.predictedStd,'Mean Abs: -','Median : -','Std dev :-');

 guidata(hObject, handles);
 
 % Exclude_From_VOI callback===============================================
 
function Exclude_From_ROI_Callback(hObject, ~, handles)
    
    
     if ~isempty(handles.voi)
            handles.nRegionexclude = handles.nRegionexclude+1;

% =========================================================================
% Selection of Regions to exclude from the VOI 
% =========================================================================

         if  isempty(handles.voiExclude)
            handles.rect = imrect(handles.Fieldmaps);
            setColor(handles.rect,'b');  
            handles.voiExclude =~handles.rect.createMask;
            handles.positionExclude{1}=handles.rect.getPosition;       
         else  
            handles.rectn = imrect(handles.Fieldmaps);
            setColor(handles.rectn,'b');
            
            handles.mask = ~handles.rectn.createMask;
            handles.voiExclude=handles.voiExclude.*handles.mask;
            handles.positionExclude{end+1}=handles.rectn.getPosition;
         end
    
% =========================================================================
% DEFINE SHIM VOI validity mask
% =========================================================================

    handles.voi = handles.voi.*handles.voiExclude;
    
    if (handles.fieldExpired ~=0)
    handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.fieldInspired, handles.fieldExpired ) ;
    else
    handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.fieldInspired) ;
    end
    handles.bound=bwboundaries(handles.voi);

% Adjust shim VOI based on the rectangular selection on the image==========

            for j=1:handles.dim(3);
                handles.Params.shimVoi(:,:,j)=handles.Params.shimVoi(:,:,j).*handles.voi;
            end
            
% =========================================================================
% Field distribution parameters inside the VOI
% =========================================================================
      [handles.meanRoi,handles.medianRoi,handles.stdRoi,handles.meanRoiexp,handles.medianRoiexp,handles.stdRoiexp]=calculparameters(handles);

      
% =========================================================================    
%Plot Images and field parameters inside the VOI
% =========================================================================
        switch handles.itemSelected
            
                case 'Phase/Inspired'
                    imagefieldroi(handles.fieldInspired,handles.Roi,handles.voi,handles);  
                    setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,handles.meanRoi,handles.medianRoi,handles.stdRoi)
                case 'Phase/Expired'
                    imagefieldroi(handles.fieldExpired,handles.Roi,handles.voi,handles);   
                    setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,handles.meanRoiexp,handles.medianRoiexp,handles.stdRoiexp)

                case 'Mag/Inspired'             
                    imagemagroi(handles.magInspired,handles.Roi,handles.voi,handles);
                    setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,'Mean Abs: -','Median : -','Std dev :-')
                case 'Mag/Expired'
                    imagemagroi(handles.magExpired,handles.Roi,handles.voi,handles);
                    setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,'Mean Abs: -','Median : -','Std dev :-')
         end
     else
         display ('Warning : Create a VOI before using the button "Exclude"')
     end   
    
guidata (hObject,handles);

% --- Executes on button press in Extend_SCT.
function Extend_SCT_Callback(hObject, ~, handles)
    
    if ~isempty(handles.totalVoi)
        handles.extension = imrect(handles.Fieldmaps);
        handles.sctExtension =handles.extension.createMask;
        handles.totalVoi(:,:,handles.sliceSelected)=or(handles.totalVoi(:,:,handles.sliceSelected),handles.sctExtension);
        
        for j=1:handles.dim(3);
         handles.Params.shimVoi(:,:,j)=handles.Params.shimVoi(:,:,j)+handles.sctExtension;    
         handles.bound{j}=bwboundaries(handles.totalVoi(:,:,j));        
        end
        
      [handles.meanRoi,handles.medianRoi,handles.stdRoi,handles.meanRoiexp,handles.medianRoiexp,handles.stdRoiexp]=calculparameters(handles);

      
% =========================================================================    
%Plot Images and field parameters inside the VOI
% =========================================================================
        switch handles.itemSelected
            
                case 'Phase/Inspired'
                    imagefieldroi(handles.fieldInspired,handles.Roi,handles.totalVoi(:,:,handles.sliceSelected),handles);  
                    setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,handles.meanRoi,handles.medianRoi,handles.stdRoi)
                case 'Phase/Expired'
                    imagefieldroi(handles.fieldExpired,handles.Roi,handles.totalVoi(:,:,handles.sliceSelected),handles);   
                    setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,handles.meanRoiexp,handles.medianRoiexp,handles.stdRoiexp)

                case 'Mag/Inspired'             
                    imagemagroi(handles.magInspired,handles.Roi,handles.totalVoi(:,:,handles.sliceSelected),handles);
                    setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,'Mean Abs: -','Median : -','Std dev :-')
                case 'Mag/Expired'
                    imagemagroi(handles.magExpired,handles.Roi,handles.totalVoi(:,:,handles.sliceSelected),handles);
                    setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,'Mean Abs: -','Median : -','Std dev :-')
         end
    
    else
        display('Use Sct button to calculate a VOI before trying to extend it')
    end

guidata (hObject,handles);



% --- Executes on button press in Reduce_SCT.
function Reduce_SCT_Callback(hObject, ~, handles)
    
    if ~isempty(handles.totalVoi)
        handles.reduction = imrect(handles.Fieldmaps);
        handles.sctReduction =~handles.reduction.createMask;
        handles.totalVoi(:,:,handles.sliceSelected)=handles.totalVoi(:,:,handles.sliceSelected).*handles.sctReduction;
        
        for j=1:handles.dim(3);
         handles.Params.shimVoi(:,:,j)=handles.Params.shimVoi(:,:,j).*handles.sctReduction;    
         handles.bound{j}=bwboundaries(handles.totalVoi(:,:,j));        
        end
        

      [handles.meanRoi,handles.medianRoi,handles.stdRoi,handles.meanRoiexp,handles.medianRoiexp,handles.stdRoiexp]=calculparameters(handles);

      
% =========================================================================    
%Plot Images and field parameters inside the VOI
% =========================================================================
        switch handles.itemSelected
            
                case 'Phase/Inspired'
                    imagefieldroi(handles.fieldInspired,handles.Roi,handles.totalVoi(:,:,handles.sliceSelected),handles);  
                    setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,handles.meanRoi,handles.medianRoi,handles.stdRoi)
                case 'Phase/Expired'
                    imagefieldroi(handles.fieldExpired,handles.Roi,handles.totalVoi(:,:,handles.sliceSelected),handles);   
                    setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,handles.meanRoiexp,handles.medianRoiexp,handles.stdRoiexp)

                case 'Mag/Inspired'             
                    imagemagroi(handles.magInspired,handles.Roi,handles.totalVoi(:,:,handles.sliceSelected),handles);
                    setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,'Mean Abs: -','Median : -','Std dev :-')
                case 'Mag/Expired'
                    imagemagroi(handles.magExpired,handles.Roi,handles.totalVoi(:,:,handles.sliceSelected),handles);
                    setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,'Mean Abs: -','Median : -','Std dev :-')
         end
    
    else
        display('Use Sct button to calculate a VOI before trying to Reduce it')
    end
    
    guidata (hObject,handles);



%Fieldmaps contrast selector callback======================================

function phasecontrastSelector_Callback(hObject, ~, handles)

handles.minContrast=max(handles.fieldInspired.img(:));
handles.lim = round(get(hObject,'Value')*(handles.minContrast - 1) + 1);
handles.limits=[-handles.lim handles.lim];

switchplot(handles);

guidata (hObject,handles);


%Magntitude contrast selector callback=====================================

function magcontrastSelector_Callback(hObject, ~, handles)

handles.limMag = get(hObject,'Value');
if handles.limMag <=0.01
    handles.limMag=0.011;
end

handles.limitsMag=[-handles.limMag handles.limMag];

switchplot(handles);

guidata (hObject,handles);

% Browse your folder callback==============================================

function Browsebutton_Callback(hObject, ~, handles)
    
%Select the folder with the raw data from the scan-------------------------
    
handles.rawData= uigetdir();
[handles.pathtofolder,handles.foldername]=fileparts(handles.rawData);

%Define a new folder to store the sorted data------------------------------

handles.sortedData=fullfile(handles.pathtofolder,'/sorted_data');

%Call the function to sort the data in background--------------------------

unix(sprintf('%s','/Applications/MATLAB_R2016a.app/bin/matlab -nodesktop -nojvm -r ", SortData(''',handles.rawData,''' ,''',handles.sortedData,''');" &'));
display('Load your training maps')
guidata (hObject,handles);


% Send_current Callback----------------------------------------------------

function Send_current_Callback(hObject, ~, handles)

handles.Shims.Com.setandloadallshims(handles.valuestoSend);

guidata (hObject,handles);


% Reset current Callback===================================================

function reset_current_Callback(hObject, ~, handles)
    
    handles.Shims.Com.resetallshims();

guidata (hObject,handles);

% Start Communication Callback=============================================

function start_comunication_Callback(hObject, ~, handles)

switch handles.itemSelected
       case {'Phase/Inspired','Mag/Inspired'}
            handles.currents=handles.Params.Inspired.currents.*1000;
       case {'Phase/Expired','Mag/Expired'}
            handles.currents=handles.Params.Expired.currents.*1000;      
end
   
 handles.valuestoSend = ShimComAcdc.convert_values(handles.currents,handles.Params.coeffP1,handles.Params.coeffP2);
  
 handles.Shims.Com.Open_ComPort();
  
  display('Ready to send current')
  display(handles.valuestoSend)
  
guidata (hObject,handles);


%End Communication Callback================================================
    
function end_comunication_Callback(hObject, ~, handles)
handles.Shims.Com.close_ComPort();
display('Serial communication with the coil is close now')
guidata (hObject,handles);


% =========================================================================    
%Functions Image
% =========================================================================

function imagefield(Field,ax,handles)
        
             imageToPlot = Field.img(:,:,handles.sliceSelected);
             imagesc(imageToPlot,'parent',ax);
             colormap(ax,handles.colormap);
             caxis(ax,handles.limits);
             colorbar(ax);
             
function imagefieldroi(Field,ax,mask,handles)
        
             imageToPlot = Field.img(:,:,handles.sliceSelected).*mask;
             imagesc(imageToPlot,'parent',ax);
             colormap(ax,handles.colormap);
             caxis(ax,handles.limits);
             colorbar(ax);
             
function imagemag(Field,ax,handles)
        
             imageToPlot = Field.img(:,:,handles.sliceSelected);
             imagesc(imageToPlot,'parent',ax);
             colormap(ax,gray);
             caxis(ax,handles.limitsMag);
             colorbar(ax);
             
function imagemagroi(Field,ax,mask,handles)
             imageToPlot = Field.img(:,:,handles.sliceSelected).*mask;
             imagesc(imageToPlot,'parent',ax);
             colormap(ax,gray);
             caxis(ax,handles.limitsMag);
             colorbar(ax);
   
function imageclear(lim,ax,colormaps,handles)
             imagesc(handles.clear(:,:,handles.sliceSelected),'parent',ax);
             colormap(ax,colormaps);
             caxis(ax,lim);
             colorbar(ax);
                           
function[a,b,c,d,e,f] = calculparameters(handles) 
    
      handles.voiparametersIns=assessfielddistribution( handles.fieldInspired,handles.Params.shimVoi);
      handles.meanRoi=strcat('Mean Abs=',num2str(handles.voiparametersIns.meanAbs));
      a=handles.meanRoi;
      handles.medianRoi=strcat('Median =',num2str(handles.voiparametersIns.median));
      b=handles.medianRoi;
      handles.stdRoi=strcat('Std dev = ',num2str(handles.voiparametersIns.std));
      c=handles.stdRoi;
      d=0;
      e=0;
      f=0;
            
      if (handles.fieldExpired ~=0)
      handles.voiparametersExp=assessfielddistribution( handles.fieldExpired,handles.Params.shimVoi);
      handles.voiparametersIns=assessfielddistribution( handles.fieldInspired,handles.Params.shimVoi);
      handles.meanRoi=strcat('Mean Abs=',num2str(handles.voiparametersIns.meanAbs));
      a=handles.meanRoi;
      handles.medianRoi=strcat('Median =',num2str(handles.voiparametersIns.median));
      b=handles.medianRoi;
      handles.stdRoi=strcat('Std dev = ',num2str(handles.voiparametersIns.std));
      c=handles.stdRoi;
      handles.meanRoiexp=strcat('Mean Abs =',num2str(handles.voiparametersExp.meanAbs));
      d=handles.meanRoiexp;
      handles.medianRoiexp=strcat('Median =',num2str(handles.voiparametersExp.median));
      e=handles.medianRoiexp;
      handles.stdRoiexp=strcat('Std dev =',num2str(handles.voiparametersExp.std));
      f=handles.stdRoiexp;
      end
      
             
   function switchplot(handles)
           
     switch handles.itemSelected
            
        case 'Phase/Inspired'
              imagefield(handles.fieldInspired,handles.Fieldmaps,handles);
              setparameters(handles.fieldMean,handles.fieldMedian,handles.fieldStd,handles.meanIns,handles.medianIns,handles.stdIns);
            
              if ~isempty(handles.voi)
                  imagefieldroi(handles.fieldInspired,handles.Roi,handles.voi,handles);
                  setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,handles.meanRoi,handles.medianRoi,handles.stdRoi);
              end
              
              if ~isempty(handles.totalVoi)
                  imagefieldroi(handles.fieldInspired,handles.Roi,handles.totalVoi(:,:,handles.sliceSelected),handles);
                  setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,handles.meanRoi,handles.medianRoi,handles.stdRoi);
              end
               
              if ~isempty(handles.predictedfieldInspired)
                  imagefield(handles.predictedfieldInspired,handles.Predicted,handles);  
                  setparameters(handles.predictedMean,handles.predictedMedian,handles.predictedStd,handles.meanPre,handles.medianPre,handles.stdPre);
                  setcurrents(handles.predicted1,handles.predicted2,handles.predicted3,handles.predicted4,handles.predicted5,handles.predicted6,handles.predicted7,handles.predicted8,handles.Params.Inspired.currents);

              end
           
              if any(handles.sliceSelected == handles.sliceDeleted)==1  
                   imageclear(handles.limits,handles.Roi,handles.colormap,handles);
              end
             
              if handles.fieldInspired.img(:,:,handles.sliceSelected)==0
                  if handles.nRegion >=1
                    imageclear(handles.limits,handles.Roi,handles.colormap,handles);
                  end
              end 
            
         case 'Phase/Expired'
            if (handles.fieldExpired ~=0)
               imagefield(handles.fieldExpired,handles.Fieldmaps,handles);
               setparameters(handles.fieldMean,handles.fieldMedian,handles.fieldStd,handles.meanExp,handles.medianExp,handles.stdExp);
                         
            if ~isempty(handles.voi)
                  imagefieldroi(handles.fieldExpired,handles.Roi,handles.voi,handles); 
                  setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,handles.meanRoiexp,handles.medianRoiexp,handles.stdRoiexp);
            end  
            
             if ~isempty(handles.totalVoi)
                  imagefieldroi(handles.fieldExpired,handles.Roi,handles.totalVoi(:,:,handles.sliceSelected),handles); 
                  setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,handles.meanRoiexp,handles.medianRoiexp,handles.stdRoiexp);
            end  
            
            if ~isempty(handles.predictedfieldExpired)
                  imagefield(handles.predictedfieldExpired,handles.Predicted,handles); 
                  setparameters(handles.predictedMean,handles.predictedMedian,handles.predictedStd,handles.meanPreexp,handles.medianPreexp,handles.stdPreexp);
                  setcurrents(handles.predicted1,handles.predicted2,handles.predicted3,handles.predicted4,handles.predicted5,handles.predicted6,handles.predicted7,handles.predicted8,handles.Params.Expired.currents);

            end     
          
             if any(handles.sliceSelected == handles.sliceDeleted)==1  
                   imageclear(handles.limits,handles.Roi,handles.colormap,handles);
             end
           
             if handles.fieldExpired.img(:,:,handles.sliceSelected)==0
                 if handles.nRegion >=1
                   imageclear(handles.limits,handles.Roi,handles.colormap,handles)
                 end
             end 
             else
                 display ('There is no Expired field maps');
             end    
             
        case 'Mag/Inspired'
             imagemag(handles.magInspired,handles.Fieldmaps,handles);
             setparameters(handles.fieldMean,handles.fieldMedian,handles.fieldStd,'Mean Abs: -','Median : -','Std dev :-');
             
            if ~isempty(handles.voi)
                imagemagroi(handles.magInspired,handles.Roi,handles.voi,handles);
                setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,'Mean Abs: -','Median : -','Std dev :-')
            end
                
             if ~isempty(handles.totalVoi)
                imagemagroi(handles.magInspired,handles.Roi,handles.totalVoi(:,:,handles.sliceSelected),handles);
                setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,'Mean Abs: -','Median : -','Std dev :-')
             end
            
            if ~isempty(handles.predictedfieldInspired)            
                imagefield(handles.predictedfieldInspired,handles.Predicted,handles);
                setparameters(handles.predictedMean,handles.predictedMedian,handles.predictedStd,handles.meanPre,handles.medianPre,handles.stdPre);
                setcurrents(handles.predicted1,handles.predicted2,handles.predicted3,handles.predicted4,handles.predicted5,handles.predicted6,handles.predicted7,handles.predicted8,handles.Params.Inspired.currents);

            end
            
            if any(handles.sliceSelected == handles.sliceDeleted)==1  
                 imageclear(handles.limitsMag,handles.Roi,gray,handles)
            end
             
            if handles.fieldInspired.img(:,:,handles.sliceSelected)==0
                if handles.nRegion >=1
                  imageclear(handles.limitsMag,handles.Roi,gray,handles)
                end
            end 
                 
        case 'Mag/Expired'
            if (handles.fieldExpired ~=0)
             imagemag(handles.magExpired,handles.Fieldmaps,handles);
             setparameters(handles.fieldMean,handles.fieldMedian,handles.fieldStd,'Mean Abs: -','Median : -','Std dev :-');
             
           if ~isempty(handles.voi)  
               imagemagroi(handles.magExpired,handles.Roi,handles.voi,handles);
               setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,'Mean Abs: -','Median : -','Std dev :-');
           end
           
           if ~isempty(handles.totalVoi)  
               imagemagroi(handles.magExpired,handles.Roi,handles.totalVoi(:,:,handles.sliceSelected),handles);
               setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,'Mean Abs: -','Median : -','Std dev :-');
           end
            
           if ~isempty(handles.predictedfieldExpired)
               imagefield(handles.predictedfieldExpired,handles.Predicted,handles); 
               setparameters(handles.predictedMean,handles.predictedMedian,handles.predictedStd,handles.meanPreexp,handles.medianPreexp,handles.stdPreexp);
               setcurrents(handles.predicted1,handles.predicted2,handles.predicted3,handles.predicted4,handles.predicted5,handles.predicted6,handles.predicted7,handles.predicted8,handles.Params.Expired.currents);

           end
           
           if any(handles.sliceSelected == handles.sliceDeleted)==1  
                imageclear(handles.limitsMag,handles.Roi,gray,handles)
           end
             
           if handles.fieldInspired.img(:,:,handles.sliceSelected)==0
               if handles.nRegion >=1
                 imageclear(handles.limitsMag,handles.Roi,gray,handles)
               end
           end
           else
                 display ('There is no Expired field maps');
           end    
             
             
     end
     
      if ~isempty(handles.voi)
            if nnz(handles.fieldInspired.img(:,:,handles.sliceSelected))~= 0
                  if any(handles.sliceSelected == handles.sliceDeleted)==0
                        for i=1:handles.nRegion
                              rectangle(handles.Fieldmaps,'Position',handles.position{i},'EdgeColor','r','LineWidth',1);
                        end
                        for i=1:handles.nRegionexclude
                              rectangle(handles.Fieldmaps,'Position',handles.positionExclude{i},'EdgeColor','b','LineWidth',1);
                        end
                  end
            end             
      end
         
      
    if ~isempty(handles.predictedfieldInspired) 
         if any(handles.sliceSelected == handles.sliceDeleted)==0
             if nnz(handles.fieldInspired.img(:,:,handles.sliceSelected))~= 0
                axes(handles.Predicted);
                if handles.nRegion >= 1
                    for k = 1:length(handles.bound)
                    boundary = handles.bound{k};     
                    plot(handles.Predicted,boundary(:,2),boundary(:,1), 'Color','black', 'LineWidth', 1); 
                    end
                else
                    handles.sctSeg = handles.bound{handles.sliceSelected};
                    for k = 1:length(handles.sctSeg)                                   
                    boundary = handles.sctSeg{k};
                    plot(handles.Predicted,boundary(:,2),boundary(:,1), 'Color','black', 'LineWidth', 0.8);
                    end
                end
             end
         end
    end
       

             
% =========================================================================    
%Functions Set parameters on the GUI interface
% =========================================================================  

function setparameters(location1,location2,location3,parameter1,parameter2,parameter3)
         set(location1,'string',parameter1);
         set(location2,'string',parameter2);
         set(location3,'string',parameter3);
         
function setcurrents(loc1,loc2,loc3,loc4,loc5,loc6,loc7,loc8, current)
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
