#include "hardware.h"
#include "util.h"

int8_t channels_used[NUM_B][NUM_C]  =
{
  {0, 1, 2, 3, 4, 5, 6, 7},
};

void setup() {
  initIO();
  Serial.begin(115200);
  spiInit();
  selectBoard(0);
  delay(100);
  while (!Serial) ;
  Serial.println("ready to use");
  for (int b = 0; b < NUM_B; b++) {
    selectBoard(b);
    for (int c = 0; c < NUM_C; c++) {
      zeroPoint[b][c] = 0;
      gain[b][c] = -1.6;
    }
  }

  for (int j = 0; j < NUM_B * NUM_C; j++) {
    channel_order[j] = -1;
    board_order[j] = -1;
  }
  int i = 0;
  for (int b = 0; b < NUM_B; b++) {
    for (int c = 0; c < NUM_C; c++) {
      if (channels_used[b][c] != -1) {
        channel_order[i] = channels_used[b][c];
        board_order[i] = b;
//        Serial.print("Channel no. :"); Serial.print(i); Serial.print(" --- channel order: "); Serial.print(channel_order[i]); Serial.print(" --- board order: ");Serial.print(board_order[i]); //Used for debug only
        i = i + 1;
//        Serial.println();
      } else {
        break;
      }
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
        userresetallshims() ;
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
