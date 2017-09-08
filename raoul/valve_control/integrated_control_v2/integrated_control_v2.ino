
// include the library for LCD
#include <LiquidCrystal.h>

//initialize LCD
LiquidCrystal lcd(7, 8, 9, 10, 11, 12);

int solenoidPinInput = 6; //possitive pressure valve
int solenoidPinOutput = 5; //negative pressure valve
int buttonPin1 = 4; //button 1
int buttonPin2 = 3; //button 2
int pressurePin = A5; // pressure sensor
int statePin = 2; // green led

int decision = 0;
int buttonON = 1;
int buttonOFF = 0;

int button1State1;
int button1State2;
int button2State1;
int button2State2;


int skip_cmpt;
int sampling_per;
int testSerial;


int selectNumber = 1;
int selectNumberMax = 4;
int cycle = 0;
int cycleMax = 100;
float currentPressure;
int flagMode;
int currentCycle;
int setPressure;
int setLowPressure;
int setHighPressure;
int valP = 0;
float pinVoltage;
int cutime;

int SelVolume[4] = {0, 0, 1000, 2500};
int SelModes[4] = {0, 0, 1, 2};
int pressureLimit [6] = {0, 6, 30, 6, 30}; // changed initial pressure for dynamic to 6 good eliminate

void setup() {
  Serial.begin(9600);
  setupBoard();

  digitalWrite(statePin, HIGH);

  startMessage();
  
  digitalWrite(solenoidPinOutput, HIGH);
  delay(3000);
  digitalWrite(solenoidPinOutput, LOW);
}

void loop() {

  selectNumber = 1;
  selectMode();

  // 0L volume
  if (flagMode == 1)
  {
    setEmptyBalloon();
  }

  // static 1L or 2L volume
  if (flagMode == 2)
  {
    setStaticBalloon();
  }

  // dynamic 1L - 2L
  if (flagMode == 3)
  {
    setDynamicBalloon();
  }
}
