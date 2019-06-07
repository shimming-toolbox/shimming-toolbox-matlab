#include "hardware.h"

void setup() {
  initIO();
  Serial.begin(115200);
  Serial.println("ready to use");
  spiInit();
}
int add_no, crt;
void loop() {

  char incomingByte;
  if (Serial.available() > 0) {
    incomingByte = Serial.read();
    //Serial.println(incomingByte);
    switch (incomingByte) {
      case 'a':
        selectNone();
        Serial.println("selectNone");
        break;
      case 's':
        selectDAC();
        Serial.println("selectDAC");
        break;
      case 'w':
        selectADC();
        Serial.println("selectADC");
        break;
      case 'd':
        selectBoards();
        Serial.println("selectBoard");
        break;
      case 'r':
        selectError();
        Serial.println("selectError");
        break;
      case 't':
        Serial.println("\n");
        Serial.print("ADC output:");
        print_all();
        break;
      case 'f':
        unselectBoards();
        Serial.println("unselectBoard");
        break;
      case 'g':
        //selectBoards();
        LTC2656Write(WRITE_AND_UPDATE,channelMap[7], 5000);
//        Serial.println("Write_DAC:");
        break;
      case 'z':
        LTC2656Write(WRITE_AND_UPDATE, channelMap[1], 0);
        break;
      case 'b':
        selectBoards();
        Serial.println("selectBoard");
        break;
      case 'h':
        mosi_sck_hi();
        Serial.println("mosi_sck_hi");
        break;
      case 'j':
        mosi_sck_lo();
        Serial.println("mosi_sck_lo");
        break;
      case 'e':
        add_no = Serial.parseInt();
        Serial.println();
        Serial.print("Address: "); Serial.print(add_no); Serial.print(" --- ");
        digitalWrite(boardSelect0, logic_address_real[add_no][0]);
        digitalWrite(boardSelect1, logic_address_real[add_no][1]);
        digitalWrite(boardSelect2, logic_address_real[add_no][2]);
        break;
      case 'q':
        add_no = Serial.parseInt();
        Serial.println();
        Serial.print("Address: "); Serial.print(add_no); Serial.print(" --- ");
        Serial.print(logic_address[add_no][0]);
        Serial.print(logic_address[add_no][1]);
        Serial.print(logic_address[add_no][2]);
        digitalWrite(boardSelect0, logic_address[add_no][0]);
        digitalWrite(boardSelect1, logic_address[add_no][1]);
        digitalWrite(boardSelect2, logic_address[add_no][2]);
        break;

    }
  }
}
