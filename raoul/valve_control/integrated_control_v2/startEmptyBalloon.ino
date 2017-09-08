void startEmptyBalloon()
{
	currentPressure = checkPressure();
  setPressure = pressureLimit[1];
  while (currentPressure > setPressure)
  {
    digitalWrite(solenoidPinOutput, HIGH);
    currentPressure = checkPressure();
    delay(10);
    lcd.clear();
    lcd.print("Current pressure:");
    lcd.setCursor(0,1);
    lcd.print(currentPressure);
    lcd.setCursor(5,1);
    lcd.print("mBar");
    
  }
 digitalWrite(solenoidPinOutput, LOW); 
}
