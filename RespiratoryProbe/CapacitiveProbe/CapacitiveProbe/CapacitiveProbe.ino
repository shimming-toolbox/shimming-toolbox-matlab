#include <FreqCount.h>

void setup() {
  Serial.begin(115200);
  FreqCount.begin(100);
}


void loop() {
  if (FreqCount.available()) {
    unsigned long count = FreqCount.read();
    Serial.println(count);
  }
}
