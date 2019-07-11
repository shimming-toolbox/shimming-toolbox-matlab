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
  //Initialize DAC communication
  selectNone();
  delay(100);
  selectDAC();
  delay(100);
  selectNone();

      // reset DAC correction terms 
  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++)  
  {
    isChannelCalibrationSuccessful[iCh] = false ;
    dacGain[iCh]   = 1.0 ;
    dacOffset[iCh] = 0 ; 
  }

  userresetallshims( ) ;
  // system heartbeat prints TRUE to indicate system is responsive
  usergetsystemheartbeat() ;

}

void loop() {
  char incomingByte;
  if (Serial.available() > 0) {
    incomingByte = Serial.read();
    switch (incomingByte) {

      case 'h': // prints TRUE/FALSE \n
        usergetsystemheartbeat();
        break;

      case 's'://zero all channels
        userresetallshims( ) ;
        break;

      case 'a': // prints TRUE/FALSE \n set shim currents
        usersetandloadallshims();
        break;
        
      case 'q':// Print all channel currents in A
        usergetallchannelcurrents();
        break;

      case 'c':
        calibratedaccompensation();
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


    }
  }
}
