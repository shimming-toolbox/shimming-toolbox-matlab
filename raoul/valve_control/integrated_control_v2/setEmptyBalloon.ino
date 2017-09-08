void setEmptyBalloon()
{


  currentPressure = checkPressure();
  setPressure = pressureLimit[selectNumber];
  lcd.clear();
  lcd.print("Mode : Empty");


  //Empty the balloons
  while (currentPressure > setPressure)
  {

    digitalWrite(solenoidPinOutput, HIGH);
    currentPressure = checkPressure();
    delay(10);
  }
  digitalWrite(solenoidPinOutput, LOW);
}
