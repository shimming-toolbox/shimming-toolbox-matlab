#include "SSC.h"

#include <Wire.h>

#include <math.h>       /* pow */

//  create an SSC sensor with I2C address 0x28 and power pin 8.

SSC ssc(0x28, 8);

void setup() 
{

// ------- 
// See Honeywell SSC series Datasheet for specifications
  const uint16_t AdcResolution = pow(2, 14) ;

  const float minPressure      = 0.0 ; // kPa 
  const float maxPressure      = 10.0 ; // kPa

// ------- 
  Serial.begin(9600); //

  Wire.begin();

  //  set min / max reading and pressure, see datasheet for the values for your sensor

  ssc.setMinRaw(0);
  ssc.setMaxRaw(AdcResolution - 1);

  ssc.setMinPressure(minPressure);
  ssc.setMaxPressure(maxPressure); 

  //  start the sensor
  ssc.start();
 
}



void loop() 
{

    ssc.update(); 
    Serial.println( ssc.pressure() );
    delay(10); // in ms
  
}
