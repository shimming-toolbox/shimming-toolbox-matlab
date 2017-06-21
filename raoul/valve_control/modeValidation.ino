void modeValidation() {
  decision = 0; //no decision was taken

  while (decision == 0)
  {
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
    }
    delay(100);
  }
}
