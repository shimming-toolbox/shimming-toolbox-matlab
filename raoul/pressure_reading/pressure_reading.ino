int pressurePin = A5; // pressure sensor
int skip_cmpt;
int sampling_per;
int valP = 0;
float pinVoltage;

void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
  sampling_per=-1;

}

void loop() {

  if (sampling_per <= 0)
    { sampling_per = Serial.parseInt();
      if (sampling_per > 0)
      { skip_cmpt = sampling_per / 10;
        Serial.println(skip_cmpt);
      }
    }
    else {
      Serial.println(checkPressure());
      delay(10);
      }
  // put your main code here, to run repeatedly:

}
