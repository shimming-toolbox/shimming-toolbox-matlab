void setupBoard() {
  //lcd setup
  lcd.begin(16, 2);

  //Pin setup
  pinMode(solenoidPinInput, OUTPUT);
  pinMode(solenoidPinOutput, OUTPUT);
  pinMode(statePin, OUTPUT);
  pinMode(buttonPin1, INPUT);
  pinMode(buttonPin2, INPUT);
  pinMode(pressurePin, INPUT);

}
