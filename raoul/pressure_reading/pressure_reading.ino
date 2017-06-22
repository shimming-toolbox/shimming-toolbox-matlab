/* #include "SSC.h" */

#include <Wire.h>

#include <math.h>       /* pow */

int pressurePin = A5; // pressure sensor
int pinBits;

float pinVoltage;

//  A reading of 1 bit for the ADC corresponds to 0.0048mV
const float mVPerAdcBit = (5.0 / 1024.0) ;  // milli volts per ADC bit

//  A reading of 1 mV corresponds to 20 mBar ?
const float mBarPermV = (100.0 / 5.0) ;   

void setup()
{

  Serial.begin(9600);

}

void loop() 
{

  Serial.println(checkPressure());
  delay(10); // in ms
    /* ssc.update();  */
    /* Serial.println( ssc.pressure() ); */
    /* delay(10); // in ms */

}
