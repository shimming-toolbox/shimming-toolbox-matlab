/*
  ShimControl.ino

// =========================================================================
// Updated::20180727::ryan.topfer@polymtl.ca
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
const uint8_t DAC_RESOLUTION      = 16 ; // 16-bit
const float DAC_VREF              = 1.25 ; // [units: volts]
const float DAC_PREAMP_RESISTANCE = 0.22 ; // [units: Ohms]

// 3 derived terms (for convenience):
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

  Serial.println("Setting ADC Range: +/- 2.048V  (1 bit = 1mV)");
  adc1.setGain(GAIN_TWO); //+/- 2.048V  1 bit = 1mV
  adc2.setGain(GAIN_TWO); //+/- 2.048V  1 bit = 1mV

  adc1.begin();
  adc2.begin();
  
    Serial.println("Ready to receive commands ...");
}

void(* resetFunc) (void) = 0;//declare reset function at address 0

void loop() {
  //Declaration of global variables 

  int sch;
  uint8_t iCh ;
  int element;
  int data;
  int i;
  int bytesread;
  int n_cal;
  unsigned int DACvaluetosend [8];                                                   //Array to store the 8 currents converted in DAC values
  int16_t adc0;
  int64_t a [41];                                                                    // Array to store 40 digits from serial communication

  String inString = "";
  char incomingByte;

  float val;
  float req_val;
  float cal_val [] = { -200, -100, 0, 100, 200};                                    // Value for Feedback calibration
  float p1[] = {0.65909, 0.65, 0.64547, 0.6499, 0.6591, 0.654, 0.65, 0.65};            // Coefficients from Feedback calibration
  float p2[] = {19.09, -10.908, -23.636, 20.91, 32.728, 11.308, -42.728, 34.546};
  float D [41];                                                                        // Array to store 40 digits from serial communication
  float currents [8];                                                                  //Array to store the 8 currents
  float vOut;



  if (Serial.available() > 0) {
    incomingByte = Serial.read();
  }
  switch (incomingByte) {

    case 'a':              // Update one channel input current 
      iCh = (Serial.parseInt()-1);
      //Serial.println(element);
      val = Serial.parseFloat();
      //Serial.println(val);
      req_val = (val - p2[element]) / p1[element];   
      setandloadchannel( iCh, req_val);
      break;
    /* case 'x':                 // Calibration of Adc current feedback */
    /*  */
    /*   for (sch = 1; sch <= 8; sch++) */
    /*   { */
    /*     for (n_cal = 0; n_cal < 5; n_cal++) */
    /*     { */
    /*       vOut = ((1.25 - cal_val[n_cal] * 0.001 * 0.22) * 26214); */
    /*       DAC.writeUpdateCh(sch - 1, vOut); */
    /*       delay(1000); */
    /*       if (sch < 5) */
    /*       { */
    /*         int adc0 = adc1.readADC_SingleEnded(sch - 1); */
    /*         Serial.println(((adc0 * 0.001) - 1.25) / 0.22 * 1000, 2); */
    /*         //Serial.println('\n'); */
    /*       } */
    /*       else */
    /*       { */
    /*         int adc0 = adc2.readADC_SingleEnded(sch - 5); */
    /*         Serial.println(((adc0 * 0.001) - 1.25) / 0.22 * 1000, 2); */
    /*         //Serial.println('\n'); */
    /*       } */
    /*     } */
    /*  */
    /*   } */
    /*  */
//      //Serial.print('\n'); 
//      break;

    case 'w':                // Set all channels to 0V

    for ( uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++ ) 
    {
        DAC.writeUpdateCh( iCh, ampstodac( 0.0 ) ) ;
    }
      /* DAC.writeUpdateCh(0, 1.259 * 26214.0); */
      /* DAC.writeUpdateCh(1, 1.246 * 26214.0); */
      /* DAC.writeUpdateCh(2, 1.251 * 26214.0); */
      /* DAC.writeUpdateCh(3, 1.257 * 26214.0); */
      /*  */
      /* DAC.writeUpdateCh(4, 1.259 * 26214.0); */
      /* DAC.writeUpdateCh(5, 1.253 * 26214.0); */
      /* DAC.writeUpdateCh(6, 1.244 * 26214.0); */
      /* DAC.writeUpdateCh(7, 1.262 * 26214.0); */

      break;
      
  //  case 'f':                 // Display current feedback from one channel
  //    element=Serial.parseInt();
  //    feedback(element);
  //    break;

    case 'e':                //Display raw current feedbacks from all channels
      feedback();
      break;


    case 'q':                // Querry on currents feedback (Display current with channels info)
      query();
      break;

    case 'y':                // Querry on volatge feedback
      queryVoltage();
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

void query() {
  Serial.println("\n");
  for (int chAdr = 0; chAdr <= 3; chAdr++)
  {
    int adc0 = adc1.readADC_SingleEnded(chAdr);
    Serial.print("CH "); Serial.print(chAdr + 1); Serial.print(" set to : ");
    Serial.print(((adc0 * 0.001) - 1.25) / 0.22 * 1000, 6); Serial.println(" mA");
  }
  for (int chAdr = 0; chAdr <= 3; chAdr++)
  {
    int adc0 = adc2.readADC_SingleEnded(chAdr);
    Serial.print("CH "); Serial.print(chAdr + 5); Serial.print(" set to : ");
    //Serial.print(adc0, 4); Serial.println(" V");
    Serial.print(((adc0 * 0.001) - 1.25) / 0.22 * 1000, 6); Serial.println(" mA");
  }
}

void queryVoltage() {
  for (int chAdr = 0; chAdr <= 3; chAdr++)
  {
    int adc0 = adc1.readADC_SingleEnded(chAdr);
    Serial.print("CH "); Serial.print(chAdr + 1); Serial.print(" set to : ");
    Serial.print(adc0, 4); Serial.println(" V");
    //Serial.print((adc0 * 0.001) , 6); Serial.println(" V");
  }
  for (int chAdr = 0; chAdr <= 3; chAdr++)
  {
    int adc0 = adc2.readADC_SingleEnded(chAdr);
    Serial.print("CH "); Serial.print(chAdr + 5); Serial.print(" set to : ");
    Serial.print(adc0, 4); Serial.println(" V");
    //Serial.print((adc0 * 0.001), 6); Serial.println(" V");
  }
}

void setandloadchannel(uint8_t iCh, float current) {
  float vOut;
  vOut = ((1.25 - current * 0.001 * 0.22) * 26214); // Convert current to DAC currentue

  DAC.writeUpdateCh( iCh, vOut);             // Update channel 
}

unsigned int ampstodac( float current ) {

    return round( DAC_PREAMP_RESISTANCE*current*DAC_BITSPERVOLT - DAC_ZERO ) ;

}

//void feedback(int element){
  
//    int adc0 = adc1.readADC_SingleEnded(element-1); // Read the feedback from the Adc
//   Serial.print(((adc0 * 0.001) - 1.25) / 0.22 * 1000, 6); //Display current feedback 
//}

void feedback(){
  for (int chAdr = 0; chAdr <= 3; chAdr++)
  {
    int adc0 = adc1.readADC_SingleEnded(chAdr);
    Serial.println(((adc0 * 0.001) - 1.25) / 0.22 * 1000, 6);
  }
  for (int chAdr = 0; chAdr <= 3; chAdr++)
  {
    int adc0 = adc2.readADC_SingleEnded(chAdr);
    Serial.println(((adc0 * 0.001) - 1.25) / 0.22 * 1000, 6);
  }
}
