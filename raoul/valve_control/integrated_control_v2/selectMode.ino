void selectMode() {
  decision = 0; //no decision was taken
  lcd.clear();
  lcd.setCursor(0, 0);
  lcd.print("Push B1 to");
  lcd.setCursor(0, 1);
  lcd.print("select mode");


  while (decision == 0)
  {


    delay(50);
    button1State1 = digitalRead(buttonPin1);
    if (button1State1 == LOW)
    {
      button1State2 = HIGH;
    }
    if (button1State1 == LOW && button1State2 == HIGH)
    {
      button1State2 = LOW;
      selectNumber = selectNumber + 1;

      if (selectNumber > selectNumberMax)
      {
        selectNumber = 1;
      }
      lcd.clear();
      lcd.setCursor(0, 0);
      lcd.print("Set Mode:");
      lcd.setCursor(0, 1);
      //lcd.print(SelVolume[selectNumber]);
      //lcd.setCursor(6, 1);
      //lcd.print(" mL");
      lcd.print(SelModes[selectNumber]);
      delay(50);
    }
    delay(50);
    button2State1 = digitalRead(buttonPin2);
    if (button2State1 != LOW)
    {
      button2State2 = HIGH;
    }
    if (button2State1 == LOW && button2State2 == HIGH)
    {
      button2State2 = LOW;
      decision = 1;
      if (selectNumber == 1) // empty  
      {
        flagMode = 1;
        cycleMax = 1;
      }
      else if (selectNumber == 2) //volume 1L
      {
        cycleMax = 1;
        flagMode = 2;
      }
      else if (selectNumber == 3) //volume 2L
      {
        cycleMax = 1;
        flagMode = 2;
      }
      else if (selectNumber == 4) //volume 1-2L
      {
        cycleMax = 50;
        flagMode = 3;
      }

      lcd.setCursor(0, 1);
      lcd.print("ERROR !!!");
      lcd.clear();
      lcd.setCursor(0, 0);
      lcd.print("Push B2 ");
      lcd.setCursor(0, 1);
      lcd.print("to validate");
      delay(500);

      modeValidation();
    }
  }
}
