int pressurePin = A5; // pressure sensor

void setup() {
  Serial.begin(115200);
  analogReference(EXTERNAL);

}

void loop() {
      Serial.println(analogRead(pressurePin));
      delay(10);
      }



