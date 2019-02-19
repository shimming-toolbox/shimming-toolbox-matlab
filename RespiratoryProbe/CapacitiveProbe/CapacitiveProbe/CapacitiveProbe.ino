/* FreqCount - Example with serial output
 * http://www.pjrc.com/teensy/td_libs_FreqCount.html
 *
 * Updated::20181101::ryan.topfer@polymtl.ca
 */
#include <FreqCount.h>


const unsigned long BAUDRATE = 115200 ;

// "Gate interval"/time interval between frequency measurements [units: ms]
const unsigned long SAMPLING_PERIOD = 50 ;

void setup() 
{
  Serial.begin( BAUDRATE );
  FreqCount.begin( SAMPLING_PERIOD );
}

void loop() 
{
  if ( FreqCount.available() ) 
  {
    Serial.println( FreqCount.read(), DEC ) ;
    Serial.flush( ) ; // wait for transmission before continuing 
  }
}
