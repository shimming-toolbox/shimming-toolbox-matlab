// AnalogPressureProbe
//   
//  Program awaits user to input sampling period in milliseconds, then it
//  records the HIGH/LOW trigger input (pin 2) followed by the voltage
//  on analog input (pin A5). It prints values to serial in that order,
//  separated by a space and followed by a carriage-return.
//
//  Note: In between pressure samples, the digital input is read as frequently
//  as possible (depends on sampling period). If a HIGH is read at any point
//  over the interval, the returned value will be 1, and 0 otherwise.
//
// ========================================================================= 

const int iPinPressure        = A5 ; // pressure sensor
const int iPinTrigger         = 2 ; // TTL trigger input

unsigned long samplingPeriod  = 0 ; // delay between samples [units: ms]

void setup()
{
  Serial.begin(115200) ;
  analogReference(EXTERNAL) ;

  pinMode( iPinTrigger, INPUT ) ;

  while (samplingPeriod <= 0)
  {
    if (Serial.available() > 0) 
      {
        samplingPeriod = Serial.parseInt();
      }
  }
}

void loop()
{
  unsigned long timeElapsed = millis() ;    
 
  bool trigger = false ;
 
  while ( millis() - timeElapsed < samplingPeriod )
  {
    trigger = trigger || ( HIGH == digitalRead( iPinTrigger ) ) ;
  }
  
  Serial.print( trigger, DEC ) ; // print trigger (0 or 1)
  Serial.print( " " ) ; 
  Serial.println( analogRead( iPinPressure ), DEC ) ;  // print pressure sample
}
