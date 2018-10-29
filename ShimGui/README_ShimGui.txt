Introduction:  
===============================================================================
Shim GUI is a graphic interface for shimming simulations and experiments. It allows users to display field and magnitude maps, select manually or automatically (Using Spinal Cord Toolbox) a VOI and simulate shimming in this region (Software). Then, users can open serial communication between the computer and the antenna to send shim currents (Hardware). 

Shin GUI was developed by NeuroPoly, the neuroimagery laboratory of Polytechnique Montréal.

Getting started
===============================================================================

The following lines will help you install all you need to ensure that ShimGUI is working. 

Software needed for Shim simulations
-------------------------------------------------------------------------------

First, make sure that Matlab is installed on your computer. 
Download “realtime_shiming” repository from Github, add all the scripts to your Matlab path.

After that, Install spinal cord toolbox latest version from this link : 

    https://sourceforge.net/projects/spinalcordtoolbox/

During the installation : Add SCT to environment path ?? —-> YES 

Shimming simulations
===============================================================================

Data needed for Shim simulations
-------------------------------------------------------------------------------
Some data has to be recorded in a scanner before running simulations :

  - B0 reference maps (Without any patient inside the scan)
  - B0 field map of the spinal cord with the patient inside the scan (Inspired and expired states)
	(3d maps, generally from GRE sequence)
  - Magnitude maps of the spinal cord.


Configuration of Shim simulations
-------------------------------------------------------------------------------
There are two configuration files to fulfill
with your own project parameters : shimparameters.m and ShimSpecs.m.


Serial Communication - Shimming
===============================================================================

Hardware needed for serial communication
-------------------------------------------------------------------------------

A shim micro controller to convert data from Matlab into electrical currents. 

A digital to analog converter (DAC) to transmit many shim currents at the same time
and an analog to digital converter (ADC) to receive feedback from the antenna. 

Serial communication connections : 

  - USB to serial RS232

  - RS232 to optic fibre media converters (Optic fibre is needed for real time shimming) 

  - Connections to allow communication between RS232 and your shim micro controller. 		
	(If you use arduino, see : https://www.arduino.cc/en/Tutorial/ArduinoSoftwareRS232)

Adjustment on software for Serial Communication
-------------------------------------------------------------------------------

The interface Shim GUI is designed for serial communication with a micro controller
which is already coded in its own language and which responding to commands with 
specific actions. 

.. note:: For the original project we used an Arduino Uno microcontroller

The class ShimCom set communication functions depending on the type of microcontroller 
and how it is programmed.

For a new project, you must define communication functions from ShimCom in a new script 
ShimComNameoftheproject (See ShimComAcdc, ShimComRri as examples)


How to use GUI interface
===============================================================================
Once, you’ve done these modifications you can start using the graphic interface
by following these first steps : 

FIRST Step (How to start Matlab - 2 options) 
===============================================================================
   - Start Matlab from there terminal to use Sct. You’ll have to start Matlab from terminal each time.
  
   - Run the script getenvironmentPath.sh in the terminal with the command : sh getenvironmentPath.sh 
     (it saves environment path in a textile). Now you can start Matlab from where you want. 


Loading Data (Interface use)
===============================================================================

   - Press the Button “Browse your folder” and select the folder with raw data acquisitions from the scan.  SortData.m is called in background and sort data in a 
      new folder called sorted_data. 

   - Load inspired maps and expired maps. (Select folder with Dicom files from one acquisition (ex : 08_gre_field_mapping_shim0_ins) —> MaRdi is called to extract maps from         
     Dicom files) 

   - You can change “display parameters” when maps are loaded (Point of view, slice selected and contrasts)


VOI Control(Interface use)
===============================================================================

   - Automatic VOI definition : Press Sct_deepseg or Sct_Centerline buttons to choose which segmentation you want.

   - If both Inspired and Expired maps are loaded, segmentations are combined which gives a bigger VOI.
   
   - With buttons “extend  Sct” and “reduce Sct”, you can custom the VOI define by SCT on each slice.  



   - Manual VOI definition : Press “create VOI” and select rectangular regions on the field/magnitude maps to add them to the VOI. Press “reduce VOI” to remove rectangular
     selections from VOI.You adjust the VOI on all the slices at the same time.
     
   - You can remove from the VOI specific slices with “Clear ROI on slice”.


   - Delete VOI : Press “Clear VOI” to delete the VOI and reset parameters. 
Shim Simulations (Interface use)
=============================================================================== 

   - Simulations : Press “Generate Predicted maps” to display predicted maps and shim currents from simulations.

Serial communication (Interface use)
=============================================================================== 

   - Press “Start Communication” to open communication canal between the computer and the microcontroller (ex : Arduino)

   - Press “Send Currents” to start shimming (static for the moment).

   - Press “Reset currents” to stop shimming.

   - Press “end” to close the serial port. 






