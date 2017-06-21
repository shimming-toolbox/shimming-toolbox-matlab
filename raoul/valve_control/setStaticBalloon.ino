void setStaticBalloon()
{
  currentPressure = checkPressure();
  setPressure = pressureLimit[selectNumber];
  lcd.clear();
  lcd.print("Mode : Static");
  lcd.setCursor(0, 1);
  lcd.print("Volume: ");
  lcd.setCursor(8, 1);
  lcd.print(SelVolume[selectNumber]);

  //startSerialComm();

  while (currentPressure < setPressure)
  {
    digitalWrite(solenoidPinInput, HIGH);
    currentPressure = checkPressure();
    delay(10);
  }
  digitalWrite(solenoidPinInput, LOW);
}

