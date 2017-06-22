float checkPressure() 
{

  pinBits    = analogRead(pressurePin);          // read the input pin
  pinVoltage = float ( pinBits ) * mVPerAdcBit ; // voltage on the ADC pin
  return pinVoltage * mBarPermV ;                // returns pressure in mBar

    /* return ssc.pressure() ;  */
}
