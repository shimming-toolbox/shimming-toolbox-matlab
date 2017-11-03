
#include "AD5668.h"
#include <SPI.h>
#include <Wire.h>
#include <Adafruit_ADS1015.h>

Adafruit_ADS1015 adc1 (0x49);
Adafruit_ADS1015 adc2 (0x48); /* Use thi for the 12-bit version */

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
  Serial.begin(115200);
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
  adc1.setGain(GAIN_TWO); //+/- 2.048V  1 bit = 1mV
  adc2.setGain(GAIN_TWO); //+/- 2.048V  1 bit = 1mV
  Serial.println("Setting ADC Range: +/- 2.048V  (1 bit = 1mV)");

  adc1.begin();
  adc2.begin();
  Serial.println("Starting ADC ...");
  Serial.println("Ready to receive commands ");
}

void(* resetFunc) (void) = 0;//declare reset function at address 0

void loop() {

  float val;
  int qch;
  int sch;
  int element;
  char incomingByte;
  int16_t adc0;
  int16_t measured;
  int16_t results;
  float vOut;
  char rc;
  float req_val;

  //float p1[] = {0.6498, 0.65, 0.65, 0.6394, 0.6545, 0.6591, 0.6607, 0.659};
//float p1[] = {0.6498, 0.65, 0.6454, 0.6394, 0.6409, 0.6591, 0.6607, 0.659};
float p1[] = {0.6498, 0.65, 0.6454, 0.6394, 0.65, 0.6591, 0.6454, 0.659};

  //float p2 [] = {25.44, -4.546, -24.54, 22.44, 41.82, 19.09, -20.75, 49.09};
  //float p2 [] = {25.44, -2.726, -23.64, 22.44, 46.36, 19.09, -20.75, 49.09};
float p2 [] = {25.44, -2.726, -32.73, 22.44, 56.36, 19.09, -28.18, 49.09};

  if (Serial.available() > 0) {
    incomingByte = Serial.read();
    Serial.println(incomingByte);
  }
  switch (incomingByte) {
          
    case 'u':
      sch = Serial.parseInt();
      Serial.print("CH "); Serial.println(sch);
      val = Serial.parseFloat();
      Serial.println(val, 4);
      //float vOut;
      vOut = ((1.25 - val * 0.001 * 0.22) * 26214);
      Serial.println(vOut);
      DAC.writeUpdateCh(sch - 1, vOut);
      break;
          
    case 'q':
      query();
      break;

    case 'r':
      Serial.println("\n");
      resetFunc();
      break;

    case 'a':
      element = 0;
      val = Serial.parseFloat();
      Serial.print("Channel 1 : "); Serial.print(val); Serial.println(" mA");
      req_val = (val - p2[element]) / p1[element];
      setCh(element, req_val);
      delay(5);
      query();
      break;

    case 'b':
      element = 1;
      val = Serial.parseFloat();
      Serial.print("Channel 2 : "); Serial.print(val); Serial.println(" mA");
      req_val = (val - p2[element]) / p1[element];
      setCh(element, req_val);
      delay(5);
      query();
      break;
    //req_val = (val - 2.72) / 0.65;

    case 'c':
      element = 2;
      val = Serial.parseFloat();
      Serial.print("Channel 3 : "); Serial.print(val); Serial.println(" mA");
      req_val = (val - p2[element]) / p1[element];
      setCh(element, req_val);
      delay(5);
      query();
      break;
    //req_val = (val + 32.72) / 0.6452;

    case 'd':
      element = 3;
      val = Serial.parseFloat();
      Serial.print("Channel 4 : "); Serial.print(val); Serial.println(" mA");
      req_val = (val - p2[element]) / p1[element];
      setCh(element, req_val);
      delay(5);
      query();
      break;

    case 'e':
      element = 4;
      val = Serial.parseFloat();
      Serial.print("Channel 5 : "); Serial.print(val); Serial.println(" mA");
      req_val = (val - p2[element]) / p1[element];
      setCh(element, req_val);
      delay(5);
      query();
      break;
    case 'f':
      element = 5;
      val = Serial.parseFloat();
      Serial.print("Channel 6 : "); Serial.print(val); Serial.println(" mA");
      req_val = (val - p2[element]) / p1[element];
      setCh(element, req_val);
      delay(5);
      query();
      break;
    case 'g':
      element = 6;
      val = Serial.parseFloat();
      Serial.print("Channel 7 : "); Serial.print(val); Serial.println(" mA");
      req_val = (val - p2[element]) / p1[element];
      setCh(element, req_val);
      delay(5);
      query();
      break;
    case 'h':
      element = 7;
      val = Serial.parseFloat();
      Serial.print("Channel 8 : "); Serial.print(val); Serial.println(" mA");
      req_val = (val - p2[element]) / p1[element];
      setCh(element, req_val);
      delay(5);
      query();
      break;

    case 'z':
      for (int chAdr = 0; chAdr <= 7; chAdr++)
      {
        req_val = (val - p2[chAdr]) / p1[chAdr];
        setCh(chAdr, req_val);
        delay(50);
      }
      delay (100);
      query();
      break;
  }
}

void query() {
  Serial.println("Current state");
  for (int chAdr = 0; chAdr <= 3; chAdr++)
  {
    int adc0 = adc1.readADC_SingleEnded(chAdr);
    Serial.print("CH "); Serial.print(chAdr + 1); Serial.print(" set to : ");
    //Serial.print(adc0, 4); Serial.println(" V");
    Serial.print(((adc0 * 0.001) - 1.25) / 0.22 * 1000, 6); Serial.println(" mA");
  }
  for (int chAdr = 0; chAdr <= 3; chAdr++)
  {
    int adc0 = adc2.readADC_SingleEnded(chAdr);
    Serial.print("CH "); Serial.print(chAdr + 5); Serial.print(" set to : ");
    //Serial.print(adc0, 4); Serial.println(" V");
    Serial.print(((adc0 * 0.001) - 1.25) / 0.22 * 1000, 6); Serial.println(" mA");
  }
}

void setCh(int element, float val) {
  //float val = Serial.parseFloat();
  float vOut;
  vOut = ((1.25 - val * 0.001 * 0.22) * 26214);
  //Serial.println(vOut);
  DAC.writeUpdateCh(element, vOut);
}


