/*
  ShimControl.ino

// =========================================================================
// Updated::20181022::ryan.topfer@polymtl.ca
// =========================================================================
*/

#include "AD5668.h"
#include <SPI.h>
#include <Wire.h>
#include "Adafruit_ADS1015.h"
#include <stdlib.h>     /* atoi */
#include <math.h>

//Declaration of global variables 
Adafruit_ADS1015 adc1 (0x48);
Adafruit_ADS1015 adc2 (0x49); /* Use thi for the 12-bit version */

#define mosiPin 11
#define sclkPin 13
#define ssPin 10
#define clrPin 8
#define ldacPin 7

/*
    Software SPI instance, "AD5668(mosiPin, sclkPin ,ssPin, clrPin, ldacPin);"
    mosiPin is connected to AD5668 Din pin (15) and sclkPin to AD5668 SCK
    pin (16). Remaining connections as explained above.
*/

// global variables re: DAC
AD5668 DAC = AD5668(mosiPin, sclkPin , ssPin, clrPin, ldacPin);

// 3 terms describing hardware:
const uint8_t ADC_NCHANNELS       = 4 ; 
const uint8_t ADC_RESOLUTION      = 12 ; // 12-bit
const float   ADC_RANGE_VOUT      = 2.048 ; // for GAIN_TWO setting [units: volts]
// 1 derived term for convenience:
const float ADC_VOLTSPERBIT       = ADC_RANGE_VOUT/( pow( 2.0, float(ADC_RESOLUTION-1) ) - 1.0 ) ;

// 3 terms describing hardware:
const uint8_t DAC_RESOLUTION      = 16 ; // 16-bit
const float DAC_VREF              = 1.25 ; // [units: volts]
const float DAC_PREAMP_RESISTANCE = 0.22 ; // [units: Ohms]

// 3 derived terms for convenience:
const float DAC_RANGE_VOUT  = 2.0*DAC_VREF ; // [units: volts]
const float DAC_BITSPERVOLT = ( pow( 2.0, float(DAC_RESOLUTION) ) - 1.0 )/DAC_RANGE_VOUT ; // =26214.0 [units: bit-counts]
const float DAC_ZERO        = DAC_VREF*DAC_BITSPERVOLT ; // =32767.5 [units: bit-counts]

// variables re: shim board
const uint8_t SHIM_NCHANNELS = 8 ;

void setup() {
  Serial.begin(115200);   //Baudrate of the serial communication : Maximum
  delay(100);

  Serial.println("-------Shim control board initialization--------");

  DAC.init();

  DAC.enableInternalRef(); // Uncomment this line to turn on the internal reference.

  DAC.powerDAC_Normal(B11111111); // Power up all channels normal

  Serial.println("Setting ADC Range: +/- 2.048V  (1 bit = 1.0mV)");
  adc1.setGain(GAIN_TWO); //+/- 2.048V  1 bit = 1mV
  adc2.setGain(GAIN_TWO); //+/- 2.048V  1 bit = 1mV

  adc1.begin();
  adc2.begin();

  resetallshims() ;
    
  Serial.println("Ready to receive commands ...");
}

void(* resetFunc) (void) = 0;//declare reset function at address 0

void loop() {

  int element;
  int data;
  int i;
  int bytesread;
  unsigned int DACvaluetosend [8];                                                   //Array to store the 8 currents converted in DAC values
  int64_t a [41];                                                                    // Array to store 40 digits from serial communication

  String inString = "";
  char incomingByte;

  float val;
  
  float p1[] = {0.65909, 0.65, 0.64547, 0.6499, 0.6591, 0.654, 0.65, 0.65};            // Coefficients from Feedback calibration
  float p2[] = {19.09, -10.908, -23.636, 20.91, 32.728, 11.308, -42.728, 34.546};
  float D [41];                                                                        // Array to store 40 digits from serial communication
  float currents [8];                                                                  //Array to store the 8 currents

  if (Serial.available() > 0) {
    incomingByte = Serial.read();
  }
  switch (incomingByte) 
  {

    /* case 'a':              // Update one channel input current  */
    /*   uint8_t iCh ; */
    /*  float req_val ; */
    /*    */
    /*   iCh = (uint8_t)(Serial.parseInt()- 1) ; */
    /*   val = Serial.parseFloat(); */
    /*   Serial.println(val); */
    /*   req_val = (val - p2[element]) / p1[element];    */
    /*   DAC.writeUpdateCh( iCh, ampstodac(req_val) );   */
    /*   break; */

    case 'x':                 // Calibration of Adc current feedback
    calibrateadc() ;
      break;

    case 'w':                // Set all channels to 0V
        resetallshims( ) ;

      break;
      
  //  case 'f':                 // Display current feedback from one channel
  //    element=Serial.parseInt();
  //    feedback(element);
  //    break;

    case 'e':                //Display raw current feedbacks from all channels
      feedback();
      break;

    case 'q':                
      queryallchannelcurrents();
      break;

    case 'y':                
      queryallchannelvoltages();
      break;

    case 'r':                // Reset the arduino
      Serial.println("\n");
      resetFunc();
      break;


    case 'o':                // Update all channels with 8 currents

       for ( data= 0; data <= 40; data++){
        
           bytesread=Serial.read();       // Read digits receive by serial communication one by one
           if (isDigit(bytesread)) {
           inString += (char)bytesread;}  //Save these values in a string
           a[data]=inString.toInt();
           D[data]=double(a[data]);       //Conversion from string to double
           Serial.print(D[data]);
           
           
          // Clear the string for new input:
          inString = "";
      }
      Serial.println("Coeff saved");
      // Reconstruction of the DAC values to send
      
      DACvaluetosend[1]=(unsigned int)(D[1]*10000 + D[2]*1000 + D[3]*100 + D[4]*10 + D[5]);
      DACvaluetosend[2]=(unsigned int)(D[6]*10000 + D[7]*1000 + D[8]*100 + D[9]*10 + D[10]);
      DACvaluetosend[3]=(unsigned int)(D[11]*10000 + D[12]*1000 + D[13]*100 + D[14]*10 + D[15]);
      DACvaluetosend[4]=(unsigned int)(D[16]*10000 + D[17]*1000 + D[18]*100 + D[19]*10 + D[20]);
      DACvaluetosend[5]=(unsigned int)(D[21]*10000 + D[22]*1000 + D[23]*100 + D[24]*10 + D[25]);
      DACvaluetosend[6]=(unsigned int)(D[26]*10000 + D[27]*1000 + D[28]*100 + D[29]*10 + D[30]);
      DACvaluetosend[7]=(unsigned int)(D[31]*10000 + D[32]*1000 + D[33]*100 + D[34]*10 + D[35]);
      DACvaluetosend[8]=(unsigned int)(D[36]*10000 + D[37]*1000 + D[38]*100 + D[39]*10 + D[40]);

      // Update all the channels with the reconstructed DAC values
      
      DAC.writeUpdateCh(0, DACvaluetosend[1]);
      DAC.writeUpdateCh(1, DACvaluetosend[2]);
      DAC.writeUpdateCh(2, DACvaluetosend[3]);
      DAC.writeUpdateCh(3, DACvaluetosend[4]);
      DAC.writeUpdateCh(4, DACvaluetosend[5]);
      DAC.writeUpdateCh(5, DACvaluetosend[6]);
      DAC.writeUpdateCh(6, DACvaluetosend[7]);
      DAC.writeUpdateCh(7, DACvaluetosend[8]);

      // Feedback of all channels after updating the currents
      
      //query();

   break;


  }
}

void calibrateadc( void )
{
  // test currents [units: A]
  float calibrationCurrents [] = { -0.2, -0.1, 0.0, 0.1, 0.2 } ; 
  uint8_t nCalibrationCurrents = 5 ;

  for ( uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++)
  {
    for (uint8_t iCalibrationCurrent = 0; iCalibrationCurrent < nCalibrationCurrents; iCalibrationCurrent++)
    {
      DAC.writeUpdateCh( iCh, ampstodac( calibrationCurrents[iCalibrationCurrent] ) ) ;
      delay(1000) ;

      Serial.println( querychannelcurrent( iCh ), 2);
    }
  }
}

void queryallchannelcurrents() 
{
    for (uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++)
    { 
        Serial.print("CH ") ; Serial.print(iCh + 1) ; Serial.print(" set to : ") ;
        Serial.print( querychannelcurrent( iCh ) , 5); Serial.println(" A");
    } 
}

void queryallchannelvoltages( void ) 
{
    for ( uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++ )
    {
        Serial.print("CH ") ; Serial.print(iCh + 1) ; Serial.print(" set to : ") ;
        Serial.print( querychannelvoltage( iCh ), 4); Serial.println(" V");
    } 
}

float querychannelcurrent( uint8_t iChannel ) 
{
    return (querychannelvoltage( iChannel ) - DAC_VREF)/DAC_PREAMP_RESISTANCE ;
}

float querychannelvoltage( uint8_t iChannel ) 
{
    return ADC_VOLTSPERBIT * float( readchanneladc( iChannel ) ) ; 
}



uint16_t readchanneladc( uint8_t iChannel ) 
{
    if ( iChannel < ADC_NCHANNELS )                                  
        return adc1.readADC_SingleEnded( iChannel ) ;
    else if ( iChannel < 2*ADC_NCHANNELS )
        return adc2.readADC_SingleEnded( iChannel - ADC_NCHANNELS ) ;
    else
		Serial.print("Error: Invalid channel index. ") ;
}

void resetallshims( void ) 
{
    for ( uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++ ) 
    {
        DAC.writeUpdateCh( iCh, ampstodac( 0.0 ) ) ;
    }
}






void feedback(){
  for (int chAdr = 0; chAdr <= 3; chAdr++)
  {
    uint16_t adc0 = adc1.readADC_SingleEnded(chAdr);
    Serial.println(((adc0 * 0.001) - 1.25) / 0.22 * 1000, 6);
  }
  for (int chAdr = 0; chAdr <= 3; chAdr++)
  {
    uint16_t adc0 = adc2.readADC_SingleEnded(chAdr);
    Serial.println(((adc0 * 0.001) - 1.25) / 0.22 * 1000, 6);
  }
}

unsigned int ampstodac( float current ) 
{
    return round( current*DAC_PREAMP_RESISTANCE*DAC_BITSPERVOLT - DAC_ZERO ) ;
}

/* float adctoamps( uint16_t adcCounts ) { */
/*  */
/*     return round( adcCounts */
/*     return round( DAC_PREAMP_RESISTANCE*current*DAC_BITSPERVOLT - DAC_ZERO ) ; */
/*  */
/* } */

//void feedback(int element){
  
//    int adc0 = adc1.readADC_SingleEnded(element-1); // Read the feedback from the Adc
//   Serial.print(((adc0 * 0.001) - 1.25) / 0.22 * 1000, 6); //Display current feedback 
//}

