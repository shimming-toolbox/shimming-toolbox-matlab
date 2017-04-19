
#include "AD5668.h"
#include <SPI.h>
#include <Wire.h>
#include <Adafruit_ADS1015.h>

Adafruit_ADS1015 adc1 (0x48);
Adafruit_ADS1015 adc2 (0x49); /* Use thi for the 12-bit version */

#define mosiPin 11
#define sclkPin 13
#define ssPin 10
#define clrPin 8
#define ldacPin 7


/*
    Software SPI instance, "AD5668(mosiPin, sclkPin ,ssPin, clrPin, ldacPin);"
    mosiPin is connected to AD5668 Din pin (15) and sclkPin to AD5668 SCK
    pin (16). Remaining connections as explained above.
*/
AD5668 DAC = AD5668(mosiPin, sclkPin , ssPin, clrPin, ldacPin);

void setup() {
  Serial.begin(9600);
  delay(100);

  Serial.println("-------Shim control board initilization--------");

  // initialize the DAC
  DAC.init();
  Serial.println("Initializing DAC ...");

  DAC.enableInternalRef(); // Uncomment this line to turn on the internal reference.

  Serial.println("Turning on DAC internal reference ...");

  DAC.powerDAC_Normal(B11111111); // Power up all channels normal
  Serial.println("Power up all DAC channels normal ...");

  //adc.setGain(GAIN_TWOTHIRDS); //+/- 6.144V  1 bit = 3mV
  adc1.setGain(GAIN_ONE); //+/- 4.096V  1 bit = 2mV
  adc2.setGain(GAIN_ONE); //+/- 4.096V  1 bit = 2mV
  Serial.println("Setting ADC Range: +/- 6.144V (1 bit = 3mV)");

  adc1.begin();
  adc2.begin();
  Serial.println("Starting ADC ...");
  Serial.println("Ready to receive commands ");
}
void loop() {

  float val;
  int qch;
  int sch;
  char incomingByte;
  int16_t adc0;
  int16_t results;
  float vOut;
  if (Serial.available() > 0) {
    incomingByte = Serial.read();
    Serial.println(incomingByte);
  }
  switch (incomingByte) {
    case 's':
      sch = Serial.parseInt();
      Serial.print("CH "); Serial.println(sch);
      val = Serial.parseFloat();
      Serial.println(val, 4);
      float vOut;
      vOut = (val * 26214);
      DAC.writeUpdateCh(sch - 1, vOut);
      break;

    case 't':
      val = Serial.parseFloat();
      //Serial.println(val, 4);
      vOut = (val * 26214);
      for (int chAdr = 0; chAdr <= 7; chAdr++)
      {
        DAC.writeUpdateCh(chAdr, vOut);
      }
      break;

    case 'z':
      for (int chAdr = 0; chAdr <= 7; chAdr++)
      {
        DAC.writeUpdateCh(chAdr, 0);
      }
      break;

    case 'q':
      for (int chAdr = 0; chAdr <= 3; chAdr++)
      {
        adc0 = adc1.readADC_SingleEnded(chAdr);
        Serial.print("Query CH "); Serial.print(chAdr + 1); Serial.print(" set to : ");
        Serial.print(adc0 * 0.002, 4); Serial.println(" V");
      }
      for (int chAdr = 0; chAdr <= 3; chAdr++)
      {
        adc0 = adc2.readADC_SingleEnded(chAdr);
        Serial.print("Query CH "); Serial.print(chAdr + 5); Serial.print(" set to : ");
        Serial.print(adc0 * 0.002, 4); Serial.println(" V");
      }
      break;
  }
}

