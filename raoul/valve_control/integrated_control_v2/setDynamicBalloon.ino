void setDynamicBalloon()
{

  setLowPressure = pressureLimit[selectNumber];
  setHighPressure = pressureLimit[selectNumber + 1];

  lcd.clear();
  lcd.print("Mode : Dynamic");
  lcd.setCursor(0, 1);

  currentCycle = 1;
  while (currentCycle < cycleMax)
  {

    currentPressure = checkPressure();
    //Inhale
    while (currentPressure < setHighPressure)
    {
      digitalWrite(solenoidPinInput, HIGH);
      currentPressure = checkPressure();
      delay(10);
    }
    digitalWrite(solenoidPinInput, LOW);

    //delay(4000);

    //Exhale
    currentPressure = checkPressure();
    while (currentPressure > setLowPressure)
    {
      digitalWrite(solenoidPinOutput, HIGH);
      currentPressure = checkPressure();
      delay(10);
    }
    digitalWrite(solenoidPinOutput, LOW);
    //delay(4000);
    currentCycle = currentCycle + 1;
  }
}

