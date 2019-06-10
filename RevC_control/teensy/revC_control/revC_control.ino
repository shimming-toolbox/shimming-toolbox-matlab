#include "hardware.h"
#include "util.h"

void setup() {
  initIO();
  Serial.begin(115200);
  spiInit();
  selectBoard(0);
  delay(100);
  Serial.println("ready to use");
  for (int b = 0; b < NUM_B; b++) {
    selectBoard(b);
    for (int c = 0; c < NUM_C; c++) {
      zeroPoint[b][c] = 0;
      gain[b][c] = -1.6;
    }
  }
  //Initialize DAC
  selectNone();
  delay(100);
  selectDAC();
  delay(100);
  selectNone();

}

void loop() {
  char incomingByte;
  if (Serial.available() > 0) {
    incomingByte = Serial.read();
    switch (incomingByte) {
      
      case 'c':
        calibratedaccompensation();
        break;

      case 's'://zero all channels
        userresetallshims( ) ;
        break;

      case 'q':// Print all channel currents in A
        usergetallchannelcurrents();
        break;
        
      case 'g'://used to set DAC values
        selectBoard(0);
        float crt_val = Serial.parseFloat();
        Serial.print("Write DAC: "); Serial.print(crt_val); Serial.print(" A");
        for (int i = 0; i < 8; i++) {
          LTC2656Write(WRITE_AND_UPDATE, channelMap[i], computeDacVal_I(crt_val + 0.1 * i, 0, i));
        }
        Serial.println();
        break;
        
      case 'h': // prints TRUE/FALSE \n
        usergetsystemheartbeat();
        break;
    }
  }
}
