#include <FreqCount.h>

void setup() {
  Serial.begin(115200);
  FreqCount.begin(100);
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
          Serial.println(FreqCount.read());
      }
  }
}
