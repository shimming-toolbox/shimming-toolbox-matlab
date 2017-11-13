// ARDUINO FIRMWARE FOR SHIM CONTROL BOARD
//
// v1 created by Alexandru Foias 
//
// ========================================================================= 
// Updated::20171017::ryan.topfer@polymtl.ca 
// ========================================================================= 

#include "AD5668.h"
#include <SPI.h>
#include <Wire.h>
#include <Adafruit_ADS1015.h>

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
AD5668 DAC = AD5668(mosiPin, sclkPin , ssPin, clrPin, ldacPin);

float val;
int qch;

uint8_t iCh;
char inByte;
uint16_t adcReading;
float 

int16_t results;
float vOut;



void setup() 
{
  Serial.begin(9600);
  delay(100);

  // initialize the DAC
  DAC.init();
  /* Serial.println("Initializing DAC ..."); */

  DAC.enableInternalRef(); 

  // Power up all channels normal
  DAC.powerDAC_Normal(B11111111); 

  /* Serial.println("ADC Range: +/- 4.096V (1 bit = 2mV)"); */
  adc1.setGain(GAIN_ONE); //+/- 4.096V  1 bit = 2mV
  adc2.setGain(GAIN_ONE); //+/- 4.096V  1 bit = 2mV

  adc1.begin();
  adc2.begin();
  Serial.print(true);
}


void loop() 
{

  if (Serial.available() > 0) 
  {
    inByte = Serial.read();
  
      switch (inByte) 
      {
        /* case 's': */
        /*   iCh = Serial.parseInt(); */
        /*   val = Serial.parseFloat(); */
        /*   Serial.println(val, 4); */
        /*   float vOut; */
        /*   vOut = (val * 26214); */
        /*   DAC.writeUpdateCh(iCh - 1, vOut); */
        /*   break; */
        /*  */
        /* case 't': */
        /*   val = Serial.parseFloat(); */
        /*   //Serial.println(val, 4); */
        /*   vOut = (val * 26214); */
        /*   for (int chAdr = 0; chAdr <= 7; chAdr++) */
        /*   { */
        /*     DAC.writeUpdateCh(chAdr, vOut); */
        /*   } */
        /*   break; */
        
        case 'p': // query i-th channel
          iCh  = uint8_t( Serial.parseInt() ) ;
          adcReading = adc1.readADC_SingleEnded(iCh);
          Serial.print(adcReading);         
          
        break;

        case 'z': // reset all channels to 0.0 V
          for (uint8_t iCh = 0; iCh <= 7; iCh++)
          {
            DAC.writeUpdateCh(iCh, 0);
          }
          break;
        

        /* case 'q': // query i-th channel */
        /*   for (int iCh = 0; iCh <= 3; iCh++) */
        /*   { */
        /*     adc0 = adc1.readADC_SingleEnded(iCh); */
        /*     Serial.print("Query CH "); Serial.print(iCh + 1); Serial.print(" set to : "); */
        /*     Serial.print(adc0 * 0.002, 4); Serial.println(" V"); */
        /*   } */
        /*   for (int iCh = 0; iCh <= 3; iCh++) */
        /*   { */
        /*     adc0 = adc2.readADC_SingleEnded(iCh); */
        /*     Serial.print("Query CH "); Serial.print(iCh + 5); Serial.print(" set to : "); */
        /*     Serial.print(adc0 * 0.002, 4); Serial.println(" V"); */
        /*   } */
        /*   break; */
        }
   } 
}

