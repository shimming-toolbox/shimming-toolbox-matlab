#include "AD5668.h"
#include <SPI.h>
#include <Wire.h>
#include "Adafruit_ADS1015.h"

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
  Serial.begin(115200);   //Baudrate of the serial communication : Maximum
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
  //declaration of 2 float for each channel val (current) and req_val (voltage) 
  float val0;
  float val;
  float val1;
  float val2;
  float val3;
  float val4;
  float val5;
  float val6;
  float val7;
  float val8;           
   
  float req_val;
  float req_val0;
  float req_val1;
  float req_val2;
  float req_val3;
  float req_val4;
  float req_val5;
  float req_val6;
  float req_val7;
  int sch;
  int element;
  char incomingByte;
  int16_t adc0;
  float vOut;
  int n_cal;
  float cal_val [] = { -200, -100, 0, 100, 200, 0};

  float p1[] = {0.65909, 0.65, 0.64547, 0.6499, 0.6591, 0.654, 0.65, 0.65};
  float p2[] = {19.09, -10.908, -23.636, 20.91, 32.728, 11.308, -42.728, 34.546};



  if (Serial.available() > 0) {
    incomingByte = Serial.read();
  }
  switch (incomingByte) {
    
    case 'u':                  // Set a current for one channel
      sch = Serial.parseInt();
      Serial.print("Calibrate CH "); Serial.println(sch);
      val = Serial.parseFloat();
      Serial.println(val, 4);
      //float vOut;
      vOut = ((1.25 - val * 0.001 * 0.22) * 26214);
      Serial.println(vOut);
      DAC.writeUpdateCh(sch - 1, vOut);
      break;

    case 'x':                 // Calibration of Adc current feedback

      for (sch = 1; sch <= 8; sch++)
      {
        Serial.print("Calibrate CH "); Serial.println(sch);

        for (n_cal = 0; n_cal <= 5; n_cal++)
        {
          vOut = ((1.25 - cal_val[n_cal] * 0.001 * 0.22) * 26214);
          DAC.writeUpdateCh(sch - 1, vOut);
          delay(1000);
          if (sch < 5)
          {
            int adc0 = adc1.readADC_SingleEnded(sch - 1);
            Serial.print(((adc0 * 0.001) - 1.25) / 0.22 * 1000, 2);
            Serial.print(", ");
          }
          else
          {
            int adc0 = adc2.readADC_SingleEnded(sch - 5);
            Serial.print(((adc0 * 0.001) - 1.25) / 0.22 * 1000, 2);
            Serial.print(", ");
          }
        }
        Serial.println('\n');
      }

      Serial.print('\n');
      break;

    case 'w':                // Set all channels to 0V

      DAC.writeUpdateCh(0, 1.256 * 26214.0);
      DAC.writeUpdateCh(1, 1.246 * 26214.0);
      DAC.writeUpdateCh(2, 1.242 * 26214.0);
      DAC.writeUpdateCh(3, 1.257 * 26214.0);

      DAC.writeUpdateCh(4, 1.261 * 26214.0);
      DAC.writeUpdateCh(5, 1.253 * 26214.0);
      DAC.writeUpdateCh(6, 1.235 * 26214.0);
      DAC.writeUpdateCh(7, 1.262 * 26214.0);

      break;


    case 'q':                // Querry on currents feedback 
      query();
      break;

    case 'y':                // Querry on volatge feedback
      queryVoltage();
      break;

    case 'r':                // Reset the arduino
      Serial.println("\n");
      resetFunc();
      break;

    case 'i':                // Update all the channels input current
          val0 = Serial.parseFloat();
          DAC.writeUpdateCh(0, val0); 
          val1=Serial.parseFloat();
          DAC.writeUpdateCh(1, val1);
          val2=Serial.parseFloat();
          DAC.writeUpdateCh(2, val2);
          val3=Serial.parseFloat();
          DAC.writeUpdateCh(3, val3);
          val4=Serial.parseFloat();
          DAC.writeUpdateCh(4, val4);
          val5=Serial.parseFloat();
          DAC.writeUpdateCh(5, val5);
          val6=Serial.parseFloat();
          DAC.writeUpdateCh(6, val6);
          val7=Serial.parseFloat();
          DAC.writeUpdateCh(7, val7);
          query();
       break;

 
    case 'a':              // Update channel 1 input current 
      element = 0;
      val = Serial.parseFloat();
      //Serial.print("Channel 1 : "); Serial.print(val); Serial.println(" mA");
      req_val = (val - p2[element]) / p1[element];
      setCh(element, req_val);
      break;

    case 'b':              // Update channel 2 input current 
      element = 1;
      val = Serial.parseFloat();
      //Serial.print("Channel 2 : "); Serial.print(val); Serial.println(" mA");
      req_val = (val - p2[element]) / p1[element];
      setCh(element, req_val);
      break;
    //req_val = (val - 2.72) / 0.65;

    case 'c':             // Update channel 3 input current 
      element = 2;
      val = Serial.parseFloat();
      //Serial.print("Channel 3 : "); Serial.print(val); Serial.println(" mA");
      req_val = (val - p2[element]) / p1[element];
      setCh(element, req_val);
      break;
    //req_val = (val + 32.72) / 0.6452;

    case 'd':            // Update channel 4 input current 
      element = 3;
      val = Serial.parseFloat();
      //Serial.print("Channel 4 : "); Serial.print(val); Serial.println(" mA");
      req_val = (val - p2[element]) / p1[element];
      setCh(element, req_val);
      break;

    case 'e':            // Update channel 5 input current 
      element = 4;
      val = Serial.parseFloat();
      //Serial.print("Channel 5 : "); Serial.print(val); Serial.println(" mA");
      req_val = (val - p2[element]) / p1[element];
      setCh(element, req_val);
      break;
    case 'f':            // Update channel 6 input current 
      element = 5;
      val = Serial.parseFloat();
      //Serial.print("Channel 6 : "); Serial.print(val); Serial.println(" mA");
      req_val = (val - p2[element]) / p1[element];
      setCh(element, req_val);
      break;
    case 'g':            // Update channel 7 input current 
      element = 6;
      val = Serial.parseFloat();
      //Serial.print("Channel 7 : "); Serial.print(val); Serial.println(" mA");
      req_val = (val - p2[element]) / p1[element];
      setCh(element, req_val);
      break;
    case 'h':            // Update channel 8 input current 
      element = 7;
      val = Serial.parseFloat();
      //Serial.print("Channel 8 : "); Serial.print(val); Serial.println(" mA");
      req_val = (val - p2[element]) / p1[element];
      setCh(element, req_val);
      break;

    case 'z':            // Update all the channels with one value
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

void queryVoltage() {
  Serial.println("Current state");
  for (int chAdr = 0; chAdr <= 3; chAdr++)
  {
    int adc0 = adc1.readADC_SingleEnded(chAdr);
    Serial.print("CH "); Serial.print(chAdr + 1); Serial.print(" set to : ");
    Serial.print(adc0, 4); Serial.println(" V");
    //Serial.print((adc0 * 0.001) , 6); Serial.println(" V");
  }
  for (int chAdr = 0; chAdr <= 3; chAdr++)
  {
    int adc0 = adc2.readADC_SingleEnded(chAdr);
    Serial.print("CH "); Serial.print(chAdr + 5); Serial.print(" set to : ");
    Serial.print(adc0, 4); Serial.println(" V");
    //Serial.print((adc0 * 0.001), 6); Serial.println(" V");
  }
}

void setCh(int element, float val) {
  //float val = Serial.parseFloat();
  float vOut;
  vOut = ((1.25 - val * 0.001 * 0.22) * 26214);
  //Serial.println(vOut);
  DAC.writeUpdateCh(element, vOut);
}

void feedback(int element){
    int adc0 = adc1.readADC_SingleEnded(element);
    Serial.print("CH "); Serial.print(element + 1); Serial.print(" set to : ");
    Serial.print(((adc0 * 0.001) - 1.25) / 0.22 * 1000, 6); Serial.println(" mA");
}
