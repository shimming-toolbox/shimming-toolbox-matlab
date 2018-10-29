function varargout = ShimGui(varargin)
%SHIMGUI 

%      This function allows to check insp/exp mag and phase data, define
%      Shim VOIs, optimizing shim currents and send these currents 
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

%       You need Spinal cord toolbox installed in your computer to use the
%       options with Sct. 


% Last Modified by GUIDE v2.5 28-Mar-2018 19:49:17


% Begin initialization code - DO NOT EDIT


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ShimGui_OpeningFcn, ...
                   'gui_OutputFcn',  @ShimGui_OutputFcn, ...
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
% PARAMETERS CONFIGURATION
%==========================================================================

% --- Executes just before ShimGui is made visible.------------------------

function varargout = ShimGui_OutputFcn(~, ~, handles) 
% ShimGui_OutputFcn :
% - Define varargout
%--------------------------------------------------------------------------
varargout{1} = handles.output;


function ShimGui_OpeningFcn(hObject, ~, handles,varargin)   
% ShimGui_OpeningFcn :
%
% - Declare parameters from shimparameters.m 
% - Declare handles variables use in ShimGui functions
%
%--------------------------------------------------------------------------


handles.Shims = varargin{1};  % handles.Shims refers now the input instance of ShimUse( )

handles.Params = handles.Shims.Params ; % copy


handles.output = hObject;                           % Default command line output for ShimGui


handles.itemSelected = 'Mag/Inspired';              %View selected
handles.lim=200;                                    %Constrast for phase images
handles.limits=[-handles.lim handles.lim];
handles.limMag=0.5;                                 %Constrast for magnitude images
handles.limitsMag=[-handles.limMag handles.limMag]; 
handles.colormap='parula';                          %Color for phase images
handles.nRegion=0;                                  %Voi number of regions 
handles.nRegionexclude=0;                           
handles.voi=[];
handles.totalVoi=[];                                %Voi mask     
handles.voiExclude=[]; 
handles.position=cell(1,1);                         %Voi regions positions
handles.positionExclude=cell(1,1);
handles.bound=cell(1,1);                            %Voi boundary positions
handles.predictedfieldInspired= [];                 %Predicted fields
handles.predictedfieldExpired= [];
handles.fieldExpired=0;
handles.sliceDeleted=[];                            %Slices where user remove manually Voi                      
handles.sortedData=' ';                             %Directory to store sorted data
handles.calibrationValues=[];                      
handles.valuestoSend=[];                            %Values for shim currents
handles.middleCursor = 0.5;                         %Variable to set cursors in middle postion
handles.phaseContrast = 0.3;                        %Variable to set initial phase contrast 


if exist('environmentPath.txt')
    
    environmentPATH = textread('environmentPath.txt','%s');
    setenv('PATH', environmentPATH{1});                 % Set enviroment path with Sct and Python
end
  
guidata(hObject, handles);



%==========================================================================
%Options and Commands on the GUI Interface
%==========================================================================

% ----------- Executes on button press in "Browse your folder" ------------

function Browsebutton_Callback(hObject, ~, handles)
%Browsebutton_Callback : 
%
% - Select and sort the folder with dicom files from the scan
%--------------------------------------------------------------------------

handles.rawData= uigetdir();                        %User choose the folder
[handles.pathtofolder,handles.foldername]=fileparts(handles.rawData);

%Define a new folder to store the sorted data------------------------------

display(handles.Shims.Params.dataLoadDir);
handles.sortedData=fullfile(handles.Shims.Params.dataLoadDir,'sorted_data');

%Call function to sort the data in background------------------------------
calltoSortdata=sprintf('%s','/Applications/MATLAB_R2016a.app/bin/matlab -nodesktop -nojvm -r ", SortData(''');

sortCommand = sprintf('%s',calltoSortdata,handles.rawData,''' ,''',handles.sortedData,''');" &');

unix(sortCommand);
display('Load your training maps')
guidata (hObject,handles);


% ---------- Executes on button press in "Load Inspired Maps" -------------

function Load_training_data_Callback(hObject, ~, handles)
%Load_inspired_maps_Callback : 
%
% - Construct instance of ShimUse
% - Convert dicom files into Matlab data to construct magnitude images
% - Display these images on the GUI interface
%--------------------------------------------------------------------------


handles.Shims.loadtrainingdata();
display(handles.Shims.Data.Img{1,1,1});

                

%Set the slice selected and the contrasts on the images--------------------

handles.dim = size(handles.Shims.Data.Img{1,1,1}.img);                             %Dimension of the fieldmaps
handles.sliceSelected=round(handles.dim(3)*0.5);                           %Slice selected : Middle of the third dimensions.                       
handles.clear=zeros(size(handles.Shims.Data.Img{1,1,1}.img)); 


set(handles.sliceSelector,'value',handles.middleCursor);                                    
set(handles.Magnitude_contrast,'value',handles.middleCursor);
set(handles.Fieldmaps_contrast,'value',handles.phaseContrast);
set(handles.commandLine,'string',num2str(handles.sliceSelected));

 
handles.Shims.Data.Img{1,1,1}.img=handles.Shims.Data.Img{1,1,1}.img-handles.middleCursor;      %Use middleCursor to shift magnitude image from [0,1] to [-0.5,0.5]


%Display Inspired Fieldmap ------------------------------------------------

imagemag(handles.Shims.Data.Img{1,1,1},handles.Fieldmaps,handles);





%Display the from the background script---------------------------

%File=fopen('background');
%a=fscanf(File,'%c');
%disp(a);
%fclose(File);

guidata(hObject, handles) ;


% ---------- Executes on button press in "Load Expired Maps" --------------


% --- Executes on button press in Process_training_data.
function Process_training_data_Callback(hObject, eventdata, handles)
%
% - Convert dicom files into Matlab data to construct phase and magnitude images
%--------------------------------------------------------------------------

handles.Shims.processtrainingdata();


%Definition of new variables from the fieldmaps----------------------------

handles.Params.scaling = [min(handles.Shims.Data.Img{1,3,1}.img) max(handles.Shims.Data.Img{1,3,1}.img)]; 
handles.totalVoi=handles.Shims.Params.cordVoi;

handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask({handles.Shims.Data.Img{1,3,1}, handles.Shims.Data.Img{1,3,2}});

 for j=1:handles.dim(3)    
     handles.bound{j}=bwboundaries(handles.totalVoi(:,:,j));        
 end
    

%Declaration of the field distribution parameters -------------------------

fieldParametersexp=handles.Shims.Data.Img{1,3,1}.assessfielddistribution(handles.Shims.Params.cordVoi);
handles.meanIns=strcat('Mean Abs=',num2str(fieldParametersexp.meanAbs));
handles.medianIns=strcat('Median =',num2str(fieldParametersexp.median));
handles.stdIns=strcat('Std dev =',num2str(fieldParametersexp.std));
setparameters(handles.fieldMean,handles.fieldMedian,handles.fieldStd,handles.meanIns,handles.medianIns,handles.stdIns);

imagefield(handles.Shims.Data.Img{1,3,1},handles.Fieldmaps,handles);
imagefieldroi(handles.Shims.Data.Img{1,3,1},handles.Roi,handles.totalVoi(:,:,handles.sliceSelected),handles);



%fieldParametersins=assessfielddistribution( handles.Shims.Data.Img{1,1,1}.img);       %Calculate statictical parameters from the field distribution
%handles.meanIns=strcat('Mean Abs=',num2str(fieldParametersins.meanAbs));
%handles.medianIns=strcat('Median =',num2str(fieldParametersins.median));
%handles.stdIns=strcat('Std dev =',num2str(fieldParametersins.std));


guidata(hObject, handles) ;



% -------------- Executes on button press in "Create VOI" -----------------

function CreateVOI_Callback(hObject,~, handles)
%CreateVOI_Callback :
%    
% - Selection of a rectangular region on the field or magnitude map to
%   define a VOI mask : handles.voi (First use).
%
% - Addition of rectangular regions to an existing VOI mask handles.voi (Next use)
%
% - Display ROI (region of interest) of the selected slice.
% - This function define VOI on all slices
%--------------------------------------------------------------------------

         if any(handles.sliceSelected == handles.sliceDeleted)>=1
            display('Change the slice selected, this one was removed from the Voi')
         else
            handles.nRegion = handles.nRegion+1;

% Selection of the Regions inside the VOI ---------------------------------
                      
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
         
% DEFINE Validity mask for shim VOI  --------------------------------------

            if (handles.fieldExpired ~=0)
            handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.fieldInspired, handles.fieldExpired ) ;
            else
            handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.fieldInspired) ;
            end

% Adjust shim VOI based on the rectangular selection on the image ---------

            if ~isempty(handles.voiExclude);
                handles.voi = handles.voi.*handles.voiExclude;
            end          
            
            for j=1:handles.dim(3);                                        %Loop on each slice
                handles.Params.shimVoi(:,:,j)=handles.Params.shimVoi(:,:,j).*handles.voi;
            end
            
% Field distribution parameters inside the VOI ----------------------------
            
            [handles.meanRoi,handles.medianRoi,handles.stdRoi,handles.meanRoiexp,handles.medianRoiexp,handles.stdRoiexp]=calculparameters(handles);

% Define boundaries of the VOI---------------------------------------------

             handles.bound=bwboundaries(handles.voi);                       
             
% Plot Images and field parameters inside the VOI -------------------------
             
             switchplot(handles); 
          end
     
% Feedback from the backgroud script --------------------------------------     
            %feedback=fopen('background');
            %a=fscanf(feedback,'%c');
            %disp(a);
            %fclose(feedback);
guidata(hObject, handles) ;


 % ---------- Executes on button press in "Exclude from ROI" --------------
 
function Exclude_From_ROI_Callback(hObject, ~, handles)
% Exclude_From_ROI_Callback 
%
% - Use this function after creation of manually adjusted VOI 
% - Selection of rectangular areas to exclude from the VOI 
% - Recalculation of the statistical parameters in the new VOI
%
%--------------------------------------------------------------------------
    
     if ~isempty(handles.voi)
            handles.nRegionexclude = handles.nRegionexclude+1;             % New excluded region

% Selection of Regions to exclude from the VOI ----------------------------

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
    
         
% DEFINE SHIM VOI validity mask -------------------------------------------
    
    if (handles.fieldExpired ~=0)
    handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.fieldInspired, handles.fieldExpired ) ;
    else
    handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.fieldInspired) ;
    end
    
% Adjust shim VOI based on the rectangular selection on the image----------
                
            handles.voi = handles.voi.*handles.voiExclude;
            handles.bound=bwboundaries(handles.voi);
            
            for j=1:handles.dim(3);
                handles.Params.shimVoi(:,:,j)=handles.Params.shimVoi(:,:,j).*handles.voi;
            end
            
% Field distribution parameters inside the VOI ----------------------------

[handles.meanRoi,handles.medianRoi,handles.stdRoi,handles.meanRoiexp,handles.medianRoiexp,handles.stdRoiexp]=calculparameters(handles);

         
%Plot Images and field parameters inside the VOI --------------------------

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


% -------------- Executes on button press in "Sct_deepseg" ----------------

%function Sct_deepseg_Callback(hObject, ~, handles) 
% Sct_deepseg_Callback :
%            
% - Automatic segmentation of the spinal cord using sct_deepseg_sc from spinalCord toolbox
% - Definition of the Voi mask from the segmentation
% - Display ROI (region of interest) of the selected slice.
%--------------------------------------------------------------------------


%Convert dicom files to .nii-----------------------------------------------
              
%            dicm2nii(handles.pathtomagInspired,handles.Params.matlabPath,0); 
            
%Call Spinal cord toolbox for segmentation---------------------------------
%            [~,~] = unix(handles.Params.command);
            
%Find the name of the segmentation files-----------------------------------
%            dicomfilesIns=dir( [ handles.pathtomagInspired '/*.dcm'] );
%            filesheaderIns = dicominfo( [handles.pathtomagInspired '/' dicomfilesIns(1).name]);   
%            seriesnameIns=filesheaderIns.SeriesDescription;
            
%            sctdeepsegFileIns = sprintf('%s','/',seriesnameIns,'_seg.nii');
            
            %sctCSFsefFileIns =sprintf('%s','/',seriesnameIns,'_CSF_seg.nii') (If SCT is used with CSF segmentation options)  
         
%Save the segmentation into a mask-----------------------------------------
%            handles.sctMaskins=load_untouch_nii(sctdeepsegFileIns);
            %handles.sctCSFMaskins=load_untouch_nii(sctCSFsefFileIns);
            
%Convert to double for next steps------------------------------------------
%            handles.sctMaskIns = double(handles.sctMaskins.img);
            %handles.sctCSFMaskIns=double(handles.sctCSFMaskins.img);
            
            
%Rotation of the mask to adjust it on the images---------------------------
%            for i=1:handles.dim(3)
%            handles.sctVoiIns(:,:,i)=rot90(handles.sctMaskIns(:,:,i));
            %handles.sctCSFVoiIns(:,:,i)=rot90(handles.sctCSFMaskIns(:,:,i));
%            end
            
%Use of Sct on expired maps (Same process) --------------------------------
            
%            if (handles.fieldExpired ~=0)
%            handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.fieldInspired, handles.fieldExpired ) ;   
%            dicm2nii(handles.pathtomagExpired,handles.Params.matlabPath,0); 
            
%            [~,~] = unix(handles.Params.command2);
            
%            dicomfilesExp=dir( [ handles.pathtomagInspired '/*.dcm'] );
%            filesheaderExp = dicominfo( [handles.pathtomagInspired '/' dicomfilesExp(1).name]);   
%            seriesnameExp=filesheaderExp.SeriesDescription;
            
%            sctdeepsegFileExp = sprintf('%s','/',seriesnameExp,'_seg.nii');
            %sctCSFsefFileExp = sprintf('%s','/',seriesnameExp,'_CSF_seg.nii');
            
%            handles.sctMaskexp=load_untouch_nii(sctdeepsegFileExp);
            %handles.sctCSFMaskexp=load_untouch_nii(sctCSFsefFileExp);
            
%            handles.sctMaskExp = double(handles.sctMaskexp.img);
            %handles.sctCSFMaskExp=double(handles.sctCSFMaskexp.img);
            
            
%            for i=1:handles.dim(3)
%            handles.sctVoiExp(:,:,i)=rot90(handles.sctMaskExp(:,:,i));
            %handles.sctCSFVoiExp(:,:,i)=rot90(handles.sctCSFMaskExp(:,:,i));
%            end
            
%Gathering of Sct segmentations on inspired and expired maps---------------
            
%            handles.totalVoi=handles.sctVoiIns+handles.sctVoiExp;          %VOI mask from Sct segmentation
            
%            else
%            handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.fieldInspired) ;
%            handles.totalVoi=handles.sctVoiIns ; %+handles.sctCSFVoiIns;
%            end     
            
%           for j=1:handles.dim(3);                                        %Adjuts handles.Params.shimVoi
%                handles.Params.shimVoi(:,:,j)=handles.Params.shimVoi(:,:,j).*handles.totalVoi(:,:,j);    
%                 handles.bound{j}=bwboundaries(handles.totalVoi(:,:,j));        
%           end
            
%Calcul of parameters  inside the Voi defined by Sct-----------------------              
%            [handles.meanRoi,handles.medianRoi,handles.stdRoi,handles.meanRoiexp,handles.medianRoiexp,handles.stdRoiexp]=calculparameters(handles);

% Plot Images and field parameters inside the VOI -------------------------
%           switchplot(handles);         
                       
%guidata(hObject, handles) ;


% ------------- Executes on button press in "Sct_Centerline" --------------

function Sct_Centerline_Callback(hObject, ~, handles)
% Sct_Centerline_Callback : 
%            
% - Automatic segmentation of the spinal cord centerline using sct_centerline from spinalCord toolbox
% - Definition of the Voi mask with an extension from the centerline
%   segmentation
% - Display ROI (region of interest) of the selected slice.
%-------------------------------------------------------------------------- 

%Convert dicom files to .nii-----------------------------------------------
            dicm2nii(handles.pathtomagInspired,handles.Params.matlabPath,0); 
            
%Call Sct for the spinal cord centerline detection on inspired maps--------
            [~,~] = unix(handles.Params.commandbis);

%Find the name of the sct_centerline files---------------------------------
            dicomfilesIns=dir( [ handles.pathtomagInspired '/*.dcm'] );
            filesheaderIns = dicominfo( [handles.pathtomagInspired '/' dicomfilesIns(1).name]);   
            seriesnameIns=filesheaderIns.SeriesDescription;
            
            sctcenterlineFileIns = sprintf('%s','/',seriesnameIns,'_centerline_optic.nii');
            
%Save the segmentation into a mask-----------------------------------------
            handles.sctcenterlineMaskIns=load_untouch_nii( sctcenterlineFileIns);
 
%Convert to double for next steps------------------------------------------
            handles.sctCenterlineMaskIns = double( handles.sctcenterlineMaskIns.img);
            
%Rotation of the mask to adjust it on the images--------------------------- 
            for i=1:handles.dim(3)
             handles.sctcenterlineVoiIns(:,:,i)=rot90(handles.sctCenterlineMaskIns(:,:,i));
            end
 
            if (handles.fieldExpired ~=0)
            handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.fieldInspired, handles.fieldExpired ) ;   
            dicm2nii(handles.pathtomagExpired,handles.Params.matlabPath,0); 

%Call Sct_deepseg_sc for the spinal cord centerline detection on expired maps (Same process)
            [~,~] = unix(handles.Params.commandbis2);
          
            dicomfilesExp=dir( [ handles.pathtomagInspired '/*.dcm'] );
            filesheaderExp = dicominfo( [handles.pathtomagInspired '/' dicomfilesExp(1).name]);   
            seriesnameExp=filesheaderExp.SeriesDescription;            
            sctcenterlineFileExp = sprintf('%s','/',seriesnameExp,'_centerline_optic.nii');
            
            handles.sctcenterlineMaskexp=load_untouch_nii(sctcenterlineFileExp);            
            handles.sctCenterlineMaskExp = double(handles.sctcenterlineMaskexp.img);
            
            for i=1:handles.dim(3)
            handles.sctcenterlineVoiExp(:,:,i)=rot90(handles.sctCenterlineMaskExp(:,:,i));
            end
            
%Gathering of Sct segmentations on inspired and expired maps---------------

            handles.centralVoi=handles.sctcenterlineVoiIns+handles.sctcenterlineVoiExp;
                        
            else
            handles.Params.shimVoi = handles.Shims.Opt.getvaliditymask( handles.Params, handles.fieldInspired) ;
            handles.centralVoi=handles.sctcenterlineVoiIns;
            end
           
%3D extension of the VOI mask defined by sct_get_centerline ---------------
            structuralElement = strel('cube',3);
            handles.totalVoi=imdilate(handles.centralVoi,structuralElement);

 
            for j=1:handles.dim(3)
            handles.Params.shimVoi(:,:,j)=handles.Params.shimVoi(:,:,j).*handles.totalVoi(:,:,j);    
            handles.bound{j}=bwboundaries(handles.totalVoi(:,:,j));        
            end

%Calcul of parameters inside the Voi defined with Sct ---------------------               
            [handles.meanRoi,handles.medianRoi,handles.stdRoi,handles.meanRoiexp,handles.medianRoiexp,handles.stdRoiexp]=calculparameters(handles);
            
% Plot Images and field parameters inside the VOI -------------------------
            switchplot(handles); 
            
 guidata(hObject, handles) ;


%----------------Executes on "slice selector" movement---------------------

function sliceSelector_Callback(hObject,~, handles)
%sliceSelector_Callback
%
% - Modify the selected slice 
% - Display maps (training maps, and predicted maps if already calculated) and ROI of the new slice
%--------------------------------------------------------------------------

%Number of slice in the third dimension (specific to the sequence)---------

nSlice = handles.dim(3);                   

%Define the slice selected from the slider position -----------------------
handles.sliceSelected = round(get(hObject,'Value')*(nSlice - 1) + 1);


%Display the slice selected on the Shim GUI interface ---------------------
set(handles.commandLine,'string',num2str(round(get(hObject,'Value')*(nSlice - 1) + 1)));

switchplot(handles);
            
guidata(hObject, handles) ;



%------------------------View Selected Callback----------------------------

function View_selected_Callback(hObject, ~, handles)
%View Selected Callback:
%
% - Allow user to choose one of the four maps that will be displayed in GUI
%   interface
%
%       - Phase/Inspired
%       - Phase/Expired
%       - Mag/Inspired
%       - Mag/Expired
%
% - Note that predicted maps from simulations are only phase maps
%--------------------------------------------------------------------------
items = get(hObject,'String');
indexSelected = get(hObject,'Value');
handles.itemSelected = items{indexSelected};

switchplot(handles);
            
guidata(hObject, handles);



%---------Executes on button press in "Generate Predicted maps"------------

function Generate_predicted_maps_Callback(hObject, ~, handles)
%Generate_predicted_maps_Callback :
%
%
% - Calculate shim currents for shimming inside the VOI.
% - Simulate predicted fieldmaps with shim current inside the antenna
% - Display Shim currents, predicted maps and statistical parameters inside
%   the VOI after shim simulations
%
%--------------------------------------------------------------------------

%Set original parameters for shim optimization-----------------------------
handles.Shims.Opt.setoriginalfield( handles.Shims.Data.Img{1,3,1}.img ) ;
handles.Shims.Opt.setshimvolumeofinterest( handles.Params.shimVoi) ;

%Shim currents optimization -----------------------------------------------
handles.Params.isSolvingAugmentedSystem    = true ;
handles.Params.isPenalizingFieldDifference = true;
handles.Params.regularizationParameter     = 0 ;


handles.Shims.Opt.optimizeshimcurrents( handles.Shims.Params) ;

    
handles.Shims.Opt.Model.currents   =  handles.Params.Inspired.currents ;


%Predicted Inspired field calculation -------------------------------------
handles.predictedfieldInspired =  handles.Shims.Opt.predictshimmedfield( ) ;   

%Mask for regions without any signal in the inspired Fieldmaps ------------
mask = handles.fieldInspired.img;
mask(mask ~= 0) = 1;

handles.predictedfieldInspired.img = handles.predictedfieldInspired.img .* mask;

%Calcul of parameters inside the Voi in the predicted maps ----------------              
handles.predictedParametersins=assessfielddistribution(handles.predictedfieldInspired,handles.Params.shimVoi);
handles.meanPre=strcat('Mean Abs =',num2str(handles.predictedParametersins.meanAbs));
handles.medianPre=strcat('Median =',num2str(handles.predictedParametersins.median));
handles.stdPre=strcat('Std dev =',num2str(handles.predictedParametersins.std));


%Predicted Expired field calculation (Same process) -----------------------

if (handles.fieldExpired ~=0)
handles.Shims.Opt.setoriginalfield( handles.fieldExpired ) ;
handles.Shims.Opt.Model.currents = handles.Params.Expired.currents ;
handles.predictedfieldExpired =  handles.Shims.Opt.predictshimmedfield( ) ;

%Mask for regions without any signal in the expired Fieldmaps -------------
mask2 = handles.fieldExpired.img;
mask2(mask2 ~= 0) = 1;
handles.predictedfieldExpired.img = handles.predictedfieldExpired.img .* mask2;

%Calcul of parameters inside the Voi in the predicted map -----------------
handles.predictedParametersexp=assessfielddistribution(handles.predictedfieldExpired,handles.Params.shimVoi);
handles.meanPreexp=strcat('Mean Abs=',num2str(handles.predictedParametersexp.meanAbs));
handles.medianPreexp=strcat('Median =',num2str(handles.predictedParametersexp.median));
handles.stdPreexp=strcat('Std dev =',num2str(handles.predictedParametersexp.std));
end  

%Plot Predicted field distribution parameters, Images and Shim currents ---             
switch handles.itemSelected
                  
   case {'Phase/Inspired','Mag/Inspired'}        
         imagefield(handles.predictedfieldInspired,handles.Predicted,handles); 
         setparameters(handles.predictedMean,handles.predictedMedian,handles.predictedStd,handles.meanPre,handles.medianPre,handles.stdPre)
         setcurrents(handles.predicted1,handles.predicted2,handles.predicted3,handles.predicted4,handles.predicted5,handles.predicted6,handles.predicted7,handles.predicted8,handles.Params.Inspired.currents);
         handles.currents=handles.Params.Inspired.currents.*1000;       
   case {'Phase/Expired','Mag/Expired'}
         imagefield(handles.predictedfieldExpired,handles.Predicted,handles);     
         setparameters(handles.predictedMean,handles.predictedMedian,handles.predictedStd,handles.meanPreexp,handles.medianPreexp,handles.stdPreexp)
         setcurrents(handles.predicted1,handles.predicted2,handles.predicted3,handles.predicted4,handles.predicted5,handles.predicted6,handles.predicted7,handles.predicted8,handles.Params.Expired.currents);
         handles.currents=handles.Params.Expired.currents.*1000;    
end

%Plot VOI boundaries on the Predicted images ------------------------------
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


%------------Executes on button press in "Clear ROI on slice"--------------

function delete_region_Callback(hObject, ~, handles)
% delete_region_Callback :
%
% - Remove the selected slice from the VOI manually
% - Calculate and display parameters in the new VOI
%
%--------------------------------------------------------------------------

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

if handles.predictedfieldInspired ~=0
   imageclear(handles.limits,handles.Predicted,handles.colormap,handles);
end
guidata(hObject,handles);



%----------------Executes on button press in "Clear VOI"-------------------

function Clear_VOI_Callback(hObject, ~, handles)
% Clear_VOI_Callback :
%
% - Delete VOI on the GUI interface
% - Clear figures with the ROI and with predicted maps
%
%--------------------------------------------------------------------------

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
     

% Reset all parameters from the previous VOI ------------------------------

handles.voi = [];                   %VOI manually selected
handles.totalVoi=[];                %Total VOI from Sct                
handles.voiExclude=[];              %Region excluded from VOI
handles.voiparametersIns=0;
handles.voiparametersExp=0;
handles.position=cell(1,1);         %VOI regions positions
handles.positionExclude=cell(1,1);  %VOI exclude regions positions
handles.nRegion=0;
handles.nRegionexclude=0;
handles.Params.shimVoi(:,:,:)=0;
handles.predictedfieldInspired=[];
handles.predictedfieldExpired=[];
handles.sliceDeleted=[];
handles.bound=[];                   % VOI boundaries
setparameters(handles.voiMean,handles.voiMedian,handles.voiStd,'Mean Abs: -','Median : -','Std dev :-');
setparameters(handles.predictedMean,handles.predictedMedian,handles.predictedStd,'Mean Abs: -','Median : -','Std dev :-');

 guidata(hObject, handles);
 

%--------------Executes on button press in "Extend_SCT_VOI"----------------

function Extend_SCT_VOI_Callback(hObject, ~, handles)
% Extend_SCT_VOI_Callback :
%
%
% - Add rectangular regions to the VOI from Sct segmentation
% - Modifications occurs only on the selected slice
%
%--------------------------------------------------------------------------
    
    if ~isempty(handles.totalVoi)
        handles.extension = imrect(handles.Fieldmaps);
        handles.sctExtension =handles.extension.createMask;
        handles.totalVoi(:,:,handles.sliceSelected)=or(handles.totalVoi(:,:,handles.sliceSelected),handles.sctExtension); %Modification of the VOI
        
        for j=1:handles.dim(3);
         handles.Params.shimVoi(:,:,j)=handles.Params.shimVoi(:,:,j)+handles.sctExtension;    
         handles.bound{j}=bwboundaries(handles.totalVoi(:,:,j));        
        end
        
      [handles.meanRoi,handles.medianRoi,handles.stdRoi,handles.meanRoiexp,handles.medianRoiexp,handles.stdRoiexp]=calculparameters(handles);

          
%Plot Images and field parameters inside the new VOI ----------------------

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
        display('Use Sct_deepseg or Sct_centerline to calculate a VOI before trying to extend it')
    end

guidata (hObject,handles);



%--------------Executes on button press in "Reduce_SCT_VOI"----------------

    function Reduce_SCT_VOI_Callback(hObject, ~, handles)
% Reduce_SCT_VOI_Callback :
%
%
% - Remove rectangular regions from the VOI from Sct segmentation
% - Modifications occurs only on the selected slice
%
%--------------------------------------------------------------------------
    
    if ~isempty(handles.totalVoi)
        handles.reduction = imrect(handles.Fieldmaps);
        handles.sctReduction =~handles.reduction.createMask;
        handles.totalVoi(:,:,handles.sliceSelected)=handles.totalVoi(:,:,handles.sliceSelected).*handles.sctReduction; %Modification of the VOI
        
        for j=1:handles.dim(3);
         handles.Params.shimVoi(:,:,j)=handles.Params.shimVoi(:,:,j).*handles.sctReduction;    
         handles.bound{j}=bwboundaries(handles.totalVoi(:,:,j));        
        end
        

 [handles.meanRoi,handles.medianRoi,handles.stdRoi,handles.meanRoiexp,handles.medianRoiexp,handles.stdRoiexp]=calculparameters(handles);

      
%Plot Images and field parameters inside the new VOI ----------------------

        switch handles.itemSelected
            
                case 'Phase/Inspired'
                    imagefieldroi(handles.Shims.Data.Img{1,3,1},handles.Roi,handles.totalVoi(:,:,handles.sliceSelected),handles);  
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
        display('Use Sct_deepseg or Sct_centerline to calculate a VOI before trying to Reduce it')
    end
    
    guidata (hObject,handles);



%-----------Executes on "Fieldmap contrast"  slider movement---------------

function Fieldmaps_contrast_Callback(hObject, ~, handles)
% Fieldmaps_contrast_Callback : 
%
% - Modify contrast on displayed fieldmaps (Training, ROI and predicted)
%
%--------------------------------------------------------------------------

handles.minContrast=max(handles.fieldInspired.img(:));
handles.lim = round(get(hObject,'Value')*(handles.minContrast - 1) + 1);

handles.limits=[-handles.lim handles.lim]; %Define contrast on images

switchplot(handles);

guidata (hObject,handles);


%-----------Executes on "Magnitude contrast"  slider movement--------------

function Magnitude_contrast_Callback(hObject, ~, handles)
% Magnitude_contrast_Callback :
%
% - Modify contrast on displayed magnitude maps (Training and ROI)
%
%--------------------------------------------------------------------------
handles.limMag = get(hObject,'Value');

if handles.limMag <=0.01
    handles.limMag=0.01; % Limitation on the maximum contrast 
end

handles.limitsMag=[-handles.limMag handles.limMag]; %Define contrast on images

switchplot(handles);

guidata (hObject,handles);


%-----------Executes on button press in "Start Communication"--------------

function Start_comunication_Callback(hObject, ~, handles)
% start_comunication_Callback : 
%
% - Open serial communication port 
% - Reset Arduino Board (Needed to send command after openning the port)
%
%--------------------------------------------------------------------------
 
%Open serial comport------------------------------------------------------
 handles.Shims.Com.opencomport();
  
 display('Ready to do calibration');
guidata (hObject,handles);

%-------------Executes on button press in "Send Currents"------------------

function Send_current_Callback(hObject, ~, handles)
% Send_current_Callback : 
%
% - Send the command to update shim currents in the coils
% - Convert Shim current into DAC value (Values to send to the board)
%
%--------------------------------------------------------------------------

%Conversion of the shim currents------------------------------------------

%handles.convertedCurrents = handles.Shims.Com.ampstodac(handles.currents,handles.Params.feedbackcalibrationcoeffx,handles.Params.feedbackcalibrationcoeffy); 
%handles.valuetoSend=round(handles.convertedCurrents);

handles.Shims.Com.setandloadallshims(handles.currents);


guidata (hObject,handles);


%-------------Executes on button press in "Reset Currents"-----------------

function Reset_current_Callback(hObject, ~, handles)
% reset_current_Callback :
%
% - Reset shim currents to zero 
%
%--------------------------------------------------------------------------
handles.Shims.Com.resetallshims();

guidata (hObject,handles);



%-----------Executes on button press in "End communication"----------------
    
function End_comunication_Callback(hObject, ~, handles)
% End_comunication_Callback : 
%
% - Close serial communication port
%
%--------------------------------------------------------------------------
handles.Shims.Com.closecomport();    
display('Serial communication with the coil is close now')
guidata (hObject,handles);




%-------------Executes on button press in "Calibrate Shim"-----------------

function Calibration_Callback(hObject, ~, handles)
 % Calibration_Callback : 
 %
 % - Automatic calculation of Adc feedback calibration coefficients
 % - Call getcalibrationcoefficient() from shimComAcdc
 % - Plot feedback values in function of reference values
 % - Linear fit of the plot and exctraction of the coefficients
 % - Convert Shim current into DAC value (Values to send to the board)
 %
 %-------------------------------------------------------------------------
 
 if strcmp( handles.Shims.Com.ComPort.Status, 'closed' ) 
     handles.Shims.Com.opencomport() ;
 end

 calibrationvalues = handles.Shims.Com.getcalibrationcoefficient();

 % TODO : convert units to A
calibrationref = [-200;-100;0;100;200];
 
 display( 'Calibration in process...(Wait 41 seconds)' )

 for i=1:8
     % f1=fit(calibrationref,calibrationvalues(:,i),'poly1');
     % feedbackcoeff= coeffvalues(f1);
     % 
     % 20180607: RT changed to polyfit (which doesn't require the MATLAB curvefitting toolbox)
     feedbackcoeff = polyfit( calibrationref, calibrationvalues(:,i), 1 );

     handles.Shims.Com.Specs.Com.feedbackcalibrationcoeffx(i)=feedbackcoeff(1);
     handles.Shims.Com.Specs.Com.feedbackcalibrationcoeffy(i)=feedbackcoeff(2);
 end
 %Conversion of the shim currents with new coefficients -------------------
 
 if ~isempty(handles.predictedfieldInspired)
    handles.convertedCurrents = handles.Shims.Com.ampstodac(handles.currents,handles.Shims.Com.feedbackcalibrationcoeffx,handles.Shims.Com.feedbackcalibrationcoeffy); 
    handles.valuetoSend=round(handles.convertedCurrents);
 end
%  dateandtime=sprintf('%s',datetime);
%  
%  diaryName=strcat(dateandtime,'_Correction_coefficient');
%  diary(diaryName);
  display('Feedback coeff X :')
  display(handles.Shims.Com.Specs.Com.feedbackcalibrationcoeffx);
  display('Feedback coeff Y :')
  display(handles.Shims.Com.Specs.Com.feedbackcalibrationcoeffy);
  display('Ready to send currents');

 guidata (hObject,handles);
 

 %---------------Executes on button press in "Feedback"--------------------

function Feedback_Callback(hObject, ~, handles)
% Feedback_Callback : 
%
% - Get currents feedback from all channels of the antenna
% - Set feedback currents values on Gui interface
%--------------------------------------------------------------------------
handles.feedbackCurrents = handles.Shims.Com.getallchanneloutputs();
display(handles.feedbackCurrents);
setcurrents(handles.feedback1,handles.feedback2,handles.feedback3,handles.feedback4,handles.feedback5,handles.feedback6,handles.feedback7,handles.feedback8,(handles.feedbackCurrents/1000));

guidata (hObject,handles);
 



% =========================================================================    
%Functions Image
% =========================================================================

function imagefield(Field,ax,handles)
% imagefield :       
%
% - Display 2D image from a fieldmap.
%
%--------------------------------------------------------------------------
             imageToPlot = Field.img(:,:,handles.sliceSelected);
             imagesc(imageToPlot,'parent',ax);
             colormap(ax,handles.colormap);
             caxis(ax,handles.limits);
             c=colorbar(ax);
             c.Label.String = 'Hz';
             c.Label.FontSize = 15;
             
function imagefieldroi(Field,ax,mask,handles)
% imagefieldroi :       
%
% - Display a 2D image ROI from a fieldmap and a mask (define ROI).
%
%--------------------------------------------------------------------------
        
             imageToPlot = Field.img(:,:,handles.sliceSelected).*mask;
             imagesc(imageToPlot,'parent',ax);
             colormap(ax,handles.colormap);
             caxis(ax,handles.limits);
             c=colorbar(ax); 
             c.Label.String = 'Hz';
             c.Label.FontSize = 15;
             
function imagemag(Field,ax,handles)
% imagemag :       
%
% - Display 2D image from a magnitude map.
%
%--------------------------------------------------------------------------        
             imageToPlot = Field.img(:,:,handles.sliceSelected);
             imagesc(imageToPlot,'parent',ax);
             colormap(ax,gray);
             %caxis(ax,handles.limitsMag);
             colorbar(ax);
             
function imagemagroi(Field,ax,mask,handles)
% imagemagroi :       
%
% - Display 2D image ROI from a magnitude map.
%
%--------------------------------------------------------------------------
             imageToPlot = Field.img(:,:,handles.sliceSelected).*mask;
             imagesc(imageToPlot,'parent',ax);
             colormap(ax,gray);
             caxis(ax,handles.limitsMag);
             colorbar(ax);
   
function imageclear(lim,ax,colormaps,handles)
% imageclear :       
%
% - Remove images from GUI interface
%
%--------------------------------------------------------------------------
             imagesc(handles.clear(:,:,handles.sliceSelected),'parent',ax);
             colormap(ax,colormaps);
             caxis(ax,lim);
             colorbar(ax);
                           
function[a,b,c,d,e,f] = calculparameters(handles)
% calculparameters :       
%
% - Calcul field distribution statistical parameters inside the VOI
%
%--------------------------------------------------------------------------
    
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
% switchplot :       
%
% - Display fieldmaps, magnitude maps, ROI maps, predicted maps, statistical parameters,
%   Roi boundaries depending on "Display parameters"  from GUI interface.
%   
%--------------------------------------------------------------------------           
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
                  handles.currents=handles.Params.Inspired.currents.*1000;


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
                  handles.currents=handles.Params.Expired.currents.*1000; 
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
                handles.currents=handles.Params.Inspired.currents.*1000;
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
               handles.currents=handles.Params.Expired.currents.*1000; 
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
                    plot(handles.Predicted,boundary(:,2),boundary(:,1), 'Color','black', 'LineWidth', 1);
                    end
                end
             end
         end
    end


      
      if ~isempty(handles.totalVoi)
        if any(handles.sliceSelected == handles.sliceDeleted)==0
             if nnz(handles.fieldInspired.img(:,:,handles.sliceSelected))~= 0
                    axes(handles.Fieldmaps);
                    hold on 
                    handles.sctSeg = handles.bound{handles.sliceSelected};                     
                    for k = 1:length(handles.sctSeg)                                   
                    boundary = handles.sctSeg{k};
                    plot(handles.Fieldmaps,boundary(:,2),boundary(:,1), 'Color','b', 'LineWidth', 1);
                    end
                    hold off
             end
        end
      end
       

             
% =========================================================================    
%Functions Set on the GUI interface
% =========================================================================  

function setparameters(location1,location2,location3,parameter1,parameter2,parameter3)
% setparameters : 
%
% Set the position to display statistical parameters on GUI interface 
%
%-------------------------------------------------------------------------- 
         set(location1,'string',parameter1);
         set(location2,'string',parameter2);
         set(location3,'string',parameter3);

         
function setcurrents(loc1,loc2,loc3,loc4,loc5,loc6,loc7,loc8, current)
% setcurrents : 
%
% - Converts currents from mA to A
% - Set the position to display shim currents on GUI interface
%
%--------------------------------------------------------------------------
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
% View_selected_CreateFcn : 
%
% - Define choices for the button View selected 
%
%--------------------------------------------------------------------------
    
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',{'Phase/Inspired';'Phase/Expired';'Mag/Inspired';'Mag/Expired'});


% ------------------Executes when user  close ShimGui----------------------

function figure1_CloseRequestFcn(hObject, eventdata, handles)
% figure1_CloseRequestFcn : 
%
% - Close matlab process in background (SortData.m) when user close Shim GUI
%
%--------------------------------------------------------------------------

[~, matlabProcess] = unix('ps -ef | grep matlab');   % Search for all matlab processes in the computer 

arrayofMatlabProcess = strsplit(matlabProcess,'\n');

for i=1:size(arrayofMatlabProcess,2)
processDescription = strsplit(arrayofMatlabProcess{i},' ');    
if (size(processDescription,2) == 11 && strcmp(processDescription{11},'no'))   % Check if the process is in background
unix( ['kill -9 ' processDescription{4}]);                 % Kill the process
end
clear pid
end
delete(hObject);



function setchannelNumber_Callback(hObject, eventdata, handles)
% hObject    handle to setchannelNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setchannelNumber as text
%        str2double(get(hObject,'String')) returns contents of setchannelNumber as a double
handles.channelNumber = str2double(get(hObject,'string'));
display(handles.channelNumber);

guidata(hObject, handles) ;



function manualCurrent_Callback(hObject, eventdata, handles)
handles.manualCurrent = str2double(get(hObject,'string'));
display(handles.manualCurrent);
guidata(hObject, handles) ;



% --- Executes on button press in Update.
function sendmanualCurrent_Callback(hObject, eventdata, handles)
    
% coeff1=handles.Params.feedbackcalibrationcoeffx(handles.channelNumber);
% coeff2=handles.Params.feedbackcalibrationcoeffy(handles.channelNumber);
% handles.Shims.Com.setandloadshim(handles.channelNumber,handles.manualCurrent,coeff1,coeff2); 
handles.Shims.Com.setandloadshim( handles.channelNumber, handles.manualCurrent ) ; 

display('Channel Updated');


