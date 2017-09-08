float checkPressure() {
  valP = analogRead(pressurePin);    // read the input pin
  pinVoltage = valP * (5 / 1024.0);      //  Calculate the voltage on the A/D pin
  //  A reading of 1 for the A/D = 0.0048mV
  //  if we multiply the A/D reading by 0.00488 then
  //  we get the voltage on the pin.
  float currentPressure = pinVoltage * 100 / 5; //in mBar
  return currentPressure;
}
