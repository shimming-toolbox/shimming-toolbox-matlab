#include "SSC.h"

#include <Wire.h>

//  create an SSC sensor with I2C address 0x28 and power pin 8.

SSC ssc(0x28, 8);

void setup() 
{

  Serial.begin(9600); //

  Wire.begin();

  //  set min / max reading and pressure, see datasheet for the values for your sensor

  ssc.setMinRaw(0);
  ssc.setMaxRaw(16383);

  ssc.setMinPressure(0.0);
  ssc.setMaxPressure(996.0);

  //  start the sensor
  ssc.start();
 
}



void loop() 
{

    ssc.update(); 
    Serial.println( ssc.pressure() );
    delay(10); // in ms
  
}
