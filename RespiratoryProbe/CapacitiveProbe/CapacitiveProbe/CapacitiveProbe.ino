/* FreqCount - Example with serial output
 * http://www.pjrc.com/teensy/td_libs_FreqCount.html
 *
 * Updated::20181029::ryan.topfer@polymtl.ca
 */
#include <FreqCount.h>


const unsigned long BAUDRATE = 57600 ;

// "Gate interval"/time interval between frequency measurements [units: ms]
const unsigned long SAMPLING_PERIOD = 100 ;
 
// delay before serial buffer is flushed for new measurement [units: ms]
unsigned long DWELL_TIME = 10 ; 

void setup() 
{
  Serial.begin(57600);
  FreqCount.begin( SAMPLING_PERIOD );
}

void loop() 
{
  if ( FreqCount.available() ) 
  {
    Serial.println( FreqCount.read() ) ;
    delay( DWELL_TIME ) ; // Matlab has dwellTime ms to read the serial port
    Serial.flush( ) ; // clear the (now outdated) measurement from serial buffer
  }
}
