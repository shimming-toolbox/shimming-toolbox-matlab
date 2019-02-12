// AnalogPressureProbe
//   
//  Reads voltage from Arduino analog input (pin A5) and prints raw value to
//  serial. 
//
// ========================================================================= 
// Updated::20181101::ryan.topfer@polymtl.ca
// ========================================================================= 

int pressurePin         = A5 ; // pressure sensor

void setup()
{
  Serial.begin(115200) ;
  analogReference(EXTERNAL) ;
  Serial.println(analogRead( pressurePin ));
}

void loop() 
{
  
  char incomingByte;

  if (Serial.available() > 0) 
  {
    incomingByte = Serial.read();
  
      switch (incomingByte) 
      {
        case 'a':
          Serial.println(analogRead( pressurePin ));
      }
  }
}
