int pressurePin = A5; // pressure sensor
unsigned long dwellTime = 10 ; // delay between samples [units: ms]

void setup()
{
  Serial.begin(115200) ;
  analogReference(EXTERNAL) ;
}

void loop()
{
  unsigned int p = analogRead( pressurePin ) ;
  Serial.println( p, DEC ) ;
  delay( dwellTime ) ; // Matlab has dwellTime ms to read the serial port
  Serial.flush() ; // clear the (now outdated) measurement
}



