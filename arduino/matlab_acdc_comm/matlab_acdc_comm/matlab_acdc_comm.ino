// send command on serial port in the form of: a XXX,
// where XXX is a float
// if XXX > 100 the led turns on, else the led turns off

int ledPin = 13;
char matlabData;
float val;

char states[8] = {'a','b','c','d','e','f','g','h'};
void setup()
{
  pinMode(ledPin, OUTPUT);
  Serial.begin(115200);
}

void loop()
{
  if (Serial.available() > 0) // if there is data to read
  {
    matlabData = Serial.read(); // read data
    if (matlabData == states[0])
    { val = Serial.parseFloat();
      Serial.println(val);
      if (val > 100) {
        digitalWrite(ledPin, HIGH); // turn light on}
      } else
      { digitalWrite(ledPin, LOW);
      }
    }
  }
}
