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

      case 'a': // prints TRUE/FALSE \n
        usersetandloadallshims();
        break;

      case 'c':  // prints TRUE/FALSE \n for each shim channel
        calibratedaccompensation();
        break;

      case 'e': // prints TRUE/FALSE \n
        usersetandloadshimbychannel();
        break;

      case 'f': // prints TRUE/FALSE \n
        usersetandloadshimbychannelasfloat();
        break;

      case 'h': // prints TRUE/FALSE \n
        usergetsystemheartbeat();
        break;

      case 'q':  // prints 5-digit precision channel current [units: A] \n for each shim channel
        usergetallchannelcurrents();
        break;

      case 'r': // prints TRUE/FALSE \n
        userresetallshims( ) ;
        break;

      case 'v':   // prints uint16_t channel voltage [units: mV] \n for each shim channel
        usergetallchannelvoltages();
        break;

      case 'u': // prints 5-digit precision DAC offset (float) \n  DAC gain \n for each shim channel
        usergetdaccompensationcoefficients() ;
        break;
      case 'w':
        //pushcurrentchannel();
        break;

        //      case 'z': // prints TRUE/FALSE \n
        //        resetFunc();
        //        break;

        //      case 'g'://used to set DAC values
        //        selectBoard(0);
        //        float crt_val = Serial.parseFloat();
        //        Serial.print("Write DAC: "); Serial.print(crt_val); Serial.print(" A");
        //        for (int i = 0; i < 8; i++) {
        //          LTC2656Write(WRITE_AND_UPDATE, channelMap[i], computeDacVal_I(crt_val + 0.1 * i, 0, i));
        //        }
        //        Serial.println();
        //        break;


    }
  }
}
