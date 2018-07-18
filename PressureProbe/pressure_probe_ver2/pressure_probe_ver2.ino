int pressurePin = A5; // pressure sensor

void setup() {
  Serial.begin(2000000);
  analogReference(EXTERNAL);

}

void loop() {
      Serial.println(analogRead(pressurePin));
      delay(1);
      }



