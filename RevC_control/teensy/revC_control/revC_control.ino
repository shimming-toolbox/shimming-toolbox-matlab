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
int add_no;
float crt_val;
void loop() {
  char incomingByte;
  if (Serial.available() > 0) {
    incomingByte = Serial.read();
    switch (incomingByte) {
      case 'c':
        for (int b = 0; b < NUM_B; b++) {
          for (int c = 0; c < NUM_C; c++) {
            calibrate_channel(b, c);
          }
          delay(500);
        }
        break;
      case 'a'://used to test selectNone
        selectNone();
        Serial.println("selectNone");
        break;
      case 's'://used to test selectDAC
        selectDAC();
        Serial.println("selectDAC");
        break;
      case 'd': //used to test selectADC
        selectADC();
        Serial.println("selectADC");
        break;
      case 'f': //used to test selectError
        selectError();
        Serial.println("selectError");
        break;
      case 'q':// Print all channel currents in A
        usergetallchannelcurrents();
        break;
      case 'g'://used to set DAC values
        selectBoard(0);
        crt_val = Serial.parseFloat();
        Serial.print("Write DAC: "); Serial.print(crt_val); Serial.print(" A");
        for (int i = 0; i < 8; i++) {
          LTC2656Write(WRITE_AND_UPDATE, channelMap[i], computeDacVal_I(crt_val + 0.1 * i, 0, i));
        }
        Serial.println();
        break;
      case 'z'://zero all channels
        for (int b = 0; b < NUM_B; b++) {
          for (int c = 0; c < NUM_C; c++) {
            selectBoard(b);
            LTC2656Write(WRITE_AND_UPDATE, channelMap[c], computeDacVal_I(0, b, c));
            delay(500);
          }
        }
        Serial.println("Zero all DAC");
        break;
      case 'o'://used to test mosi_sck pins from teensy board
        mosi_sck_hi();
        Serial.println("mosi_sck_hi");
        break;
      case 'p'://used to test mosi_sck pins from teensy board
        mosi_sck_lo();
        Serial.println("mosi_sck_lo");
        break;
      case 'n': //used to test real board address
        add_no = Serial.parseInt();
        Serial.print("Address: "); Serial.print(add_no);
        digitalWrite(boardSelect0, board_address[add_no][0]);
        digitalWrite(boardSelect1, board_address[add_no][1]);
        digitalWrite(boardSelect2, board_address[add_no][2]);
        Serial.println();
        break;
      case 'm'://used to map board address
        add_no = Serial.parseInt();
        Serial.print("Address: "); Serial.print(add_no);
        Serial.print(logic_address[add_no][0]);
        Serial.print(logic_address[add_no][1]);
        Serial.print(logic_address[add_no][2]);
        digitalWrite(boardSelect0, logic_address[add_no][0]);
        digitalWrite(boardSelect1, logic_address[add_no][1]);
        digitalWrite(boardSelect2, logic_address[add_no][2]);
        Serial.println();
        break;
      case 'h': // prints TRUE/FALSE \n
        usergetsystemheartbeat();
        break;
    }
  }
}
