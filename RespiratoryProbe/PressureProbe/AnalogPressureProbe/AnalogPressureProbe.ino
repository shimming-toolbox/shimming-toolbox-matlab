// AnalogPressureProbe
//   
//  Reads voltage from Arduino analog input (pin A5) and prints raw value to
//  serial. 
//
// ========================================================================= 
// Updated::20181101::ryan.topfer@polymtl.ca
// ========================================================================= 

int pressurePin         = A5 ; // pressure sensor
unsigned long dwellTime = 10 ; // delay between samples [units: ms]

void setup()
{
  Serial.begin(115200) ;
  analogReference(EXTERNAL) ;
}

void loop()
{
  unsigned int p = analogRead( pressurePin ) ; // read the input pin
  Serial.println( p, DEC ) ;
  delay( dwellTime ) ; // Matlab has dwellTime ms to read the serial port
  Serial.flush() ; // clear the (now outdated) measurement
}
