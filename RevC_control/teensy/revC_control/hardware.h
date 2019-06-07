#include <t3spi.h>

#define NUM_B 1//number of boards
#define NUM_C 8//number of channesl per board

typedef enum LTC26456_COMMAND {
  WRITE_TO_INPUT             = 0x00,
  UPDATE_DAC                 = 0x10,
  WRITE_TO_INPUT_UPDATE_ALL  = 0x20,
  WRITE_AND_UPDATE           = 0x30,
  POWER_DOWN                 = 0x40,
  POWER_DOWN_CHIP            = 0x50,
  SELECT_INTERNAL_REF        = 0x60,
  SELECT_EXTERNAL_REF        = 0x70,
  NOP                        = 0xF0
} LTC26456_COMMAND;


typedef enum LTC2656_ADDRESS {
  DAC_A   = 0x00,
  DAC_B   = 0x01,
  DAC_C   = 0x02,
  DAC_D   = 0x03,
  DAC_E   = 0x04,
  DAC_F   = 0x05,
  DAC_G   = 0x06,
  DAC_H   = 0x07,
  DAC_ALL = 0x0F
} LTC2656_ADDRESS;


T3SPI * SPI_MASTER;
T3SPI * SPI_SLAVE;
//communication buffers
volatile uint16_t data_tx[20] = {};
volatile uint8_t data_tx_8[20] = {};
volatile uint16_t data_rx[20] = {};
bool read_in_flight;
// volatile uint8_t data_rx[20] = {};
const bool bncOut = true;

LTC2656_ADDRESS channelMap[] = {DAC_E, DAC_F, DAC_G, DAC_H, DAC_A, DAC_B, DAC_C, DAC_D};
uint8_t channelMap_ADC[] = {0, 4, 1, 5, 2, 6, 3, 7};
const int boardMap[NUM_B] = {2};//{3,4,5,6};

float zeroPoint[NUM_B][NUM_C];
float gain[NUM_B][NUM_C];
bool calibrationStatus[NUM_B][NUM_C];
int8_t channel_order[NUM_B * NUM_C];
int8_t board_order[NUM_B * NUM_C];

const int loopsize = 1;
float output_currents[loopsize][8] =
{
  //{-0.00,-0.01,-0.17,0.03,-0.12,0.01,-0.27,-0.00,-0.17,0.36,-0.20,0.18,0.09,-0.09,-0.08,0.10,-0.58,0.26,0.11,-0.05,0.10,0.36,0.05,-0.00,-0.25,-0.18,-1.06,0.03,0.13,0.34,0.15,0.14,},
  {0.00, -0.02, -0.00, 0.33, -0.02, 0.03, 0, 0,},
};
/******************************************************/
/*********************** PROTOTYPES *******************/
/******************************************************/


void initIO(void);
void selectNone(void);
void selectDAC(void);
void selectADC(void);
void selectError(void);
void selectBoard(int board);
void spiInit(void);



void LTC2656Write(T3SPI * SPIx
                  , LTC26456_COMMAND action
                  , LTC2656_ADDRESS address
                  , uint16_t value);

uint16_t LTC1863ReadSlow(T3SPI * SPIx, uint8_t address);
/******************************************************/
/*********************** SELECTION UTIL* **************/
/******************************************************/

const int selectPin0 = 14;
const int selectPin1 = 15;
const int boardSelect2 = 16;
const int boardSelect1 = 17;
const int boardSelect0 = 18;
const int interruptPin = 6;
const int bncPin = 4;
const int CS_BB = 30;
// Alex Constants

const int sck_ms = 13;
const int mosi_pin = 11;

void initIO() {
  pinMode(selectPin0, OUTPUT);
  pinMode(selectPin1, OUTPUT);
  pinMode(boardSelect0, OUTPUT);
  pinMode(boardSelect1, OUTPUT);
  pinMode(boardSelect2, OUTPUT);
  pinMode(interruptPin, INPUT);
  pinMode(sck_ms, OUTPUT);
  pinMode(mosi_pin, OUTPUT);
  pinMode(CS_BB, OUTPUT);
  digitalWrite(CS_BB, 1);
  if (bncOut) {
    pinMode(bncPin, OUTPUT);
  } else {
    pinMode(bncPin, INPUT);
  }

}

void spiInit() {
  //SETUP digital COMS
  SPI_MASTER = new T3SPI(&KINETISK_SPI0);
  SPI_SLAVE = new T3SPI(&KINETISK_SPI1);

  SPI_MASTER->begin_MASTER(SCK, MOSI, MISO, CS0, CS_ActiveLOW);
  SPI_SLAVE->begin_SLAVE(SCK, MOSI, MISO, CS0);

  SPI_MASTER->setCTAR(CTAR_MODE0, 16, SPI_MODE0, LSB_FIRST, SPI_CLOCK_DIV16);
  //  SPI_MASTER->setCTAR(CTAR_MODE1, 8, SPI_MODE0, LSB_FIRST, SPI_CLOCK_DIV16);
  SPI_SLAVE->setCTAR_SLAVE(16, SPI_MODE1);
}

void selectNone() {
  digitalWrite(selectPin0, HIGH);
  digitalWrite(selectPin1, HIGH);
}

void selectDAC() {
  digitalWrite(selectPin0, LOW);
  digitalWrite(selectPin1, LOW);
}

void selectADC() {
  digitalWrite(selectPin0, HIGH);
  digitalWrite(selectPin1, LOW );
}

void selectError() {
  digitalWrite(selectPin0, HIGH);
  digitalWrite(selectPin1, LOW);
  //digitalWrite(selectPin0, LOW);
  //digitalWrite(selectPin1, HIGH);
}

void selectBoard(int board) {
  int b = boardMap[board];
  digitalWrite(boardSelect0, b & 0x01);
  digitalWrite(boardSelect1, b & 0x02);
  digitalWrite(boardSelect2, b & 0x04);
}

void LTC2656Write(LTC26456_COMMAND action, LTC2656_ADDRESS address, uint16_t value) {
  //NVIC_ENABLE_IRQ(IRQ_SPI1);
  selectNone();
  selectDAC();
  //  data_tx[0] = ((action | address) & 0xFF);
  //  Serial.println(action);
  //  Serial.println(address);
  //  Serial.println((action | address)& 0xFF);
  //  data_tx[1] = ((value >> 8) & 0xFF) << 8;
  //  data_tx[1] |= (value & 0xFF);
  //  SPI_MASTER->tx16(data_tx, 1, CTAR_MODE0, CS0);
  delay(2000);
  selectNone();
  delayMicroseconds(100);
}


uint16_t LTC1863ReadSlow(uint8_t address) {
  read_in_flight = true;
  NVIC_ENABLE_IRQ(IRQ_SPI1);
  selectNone();
  selectADC();
  digitalWrite(CS_BB, 0);
  data_tx[0] = ((0x80 | ((channelMap_ADC[address] << 4)) | 0x04) << 8) | (0x00);;
  SPI_MASTER->tx16(data_tx, 1, CTAR_MODE0, CS0);
  selectNone(); //toggle CS line to initialize the conversion
  delayMicroseconds(20); //wait the conversion time

  selectNone();
  selectADC();
  data_tx[0] = ((0x80 | ((address << 4)) | 0x04) << 8) | (0x00);;
  SPI_MASTER->tx16(data_tx, 1, CTAR_MODE0, CS0);
  while (SPI_SLAVE->packetCT == 0) {

  }
  digitalWrite(CS_BB, 1);
  read_in_flight = false;
  NVIC_DISABLE_IRQ(IRQ_SPI1);
  SPI_SLAVE->packetCT = 0;
  SPI_SLAVE->dataPointer = 0;
  selectNone(); //toggle CS line to initialize the conversion

  //  return (~data_rx[1]&0xFFFF)>>4; // 74HCT244
  return (data_rx[1] & 0xFFFF) >> 4; // 74HCT240
}

uint16_t LTC1863ReadSlow(uint8_t address, uint8_t n) {
  n = 1;
  uint32_t avg = 0;
  for (int i = 0; i < n; i++) {
    avg += LTC1863ReadSlow(address);
  }
  return uint16_t(avg / n);
}

//Interrupt Service Routine to handle incoming data
void spi1_isr(void) {
  cli();
  if (read_in_flight) {
    SPI_SLAVE->rx16(data_rx, 2);
  } else {
    uint16_t dump;
    SPI_SLAVE->rx16(&dump, 1);
  }
  sei();
}

uint16_t computeDacVal_I(float current) {
  Serial.println(uint16_t((65535.0 / (current / 1.66 + 2.5))));
  return uint16_t((65535.0 / (current / 1.66 + 2.5)));
}

void selectBoards() {
  digitalWrite(boardSelect0, HIGH);
  digitalWrite(boardSelect1, LOW);
  digitalWrite(boardSelect2, LOW);
}

void unselectBoards() {
  digitalWrite(boardSelect0, HIGH);
  digitalWrite(boardSelect1, LOW);
  digitalWrite(boardSelect2, HIGH);
}

void mosi_sck_hi() {
  digitalWrite(mosi_pin, HIGH);
  digitalWrite(sck_ms, HIGH);
}

void mosi_sck_lo() {
  digitalWrite(mosi_pin, LOW);
  digitalWrite(sck_ms, LOW);
}
boolean logic_address[8][3] = {{LOW, LOW, LOW},
  {LOW, LOW, HIGH},
  {LOW, HIGH, LOW},
  {LOW, HIGH, HIGH},
  {HIGH, LOW, LOW},
  {HIGH, LOW, HIGH},
  {HIGH, HIGH, LOW},
  {HIGH, HIGH, HIGH}
};
boolean logic_address_real[8][3] = {{LOW, LOW, HIGH},
  {HIGH, HIGH, HIGH},
  {LOW, HIGH, HIGH},
  {LOW, LOW, LOW},
  {HIGH, HIGH, LOW},
  {HIGH, LOW, LOW},
  {LOW, HIGH, LOW},
  {HIGH, LOW, HIGH}
};

float computeOutI(uint16_t dacVal) {
  return ((float(dacVal) * 4.096 / 4096.0) - 1.25) / 10 / 0.2;
}

void print_all() {
  Serial.println();
  for (int i = 0; i < 7; i++) {
    uint16_t data = LTC1863ReadSlow(i);
    //    Serial.print(i);
    //    Serial.print("\t");
    Serial.print(computeOutI(data), 4);
    Serial.print("\t");
  }
}
