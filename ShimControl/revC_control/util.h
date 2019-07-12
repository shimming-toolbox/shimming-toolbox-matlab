

/******************************************************/
/******************* FUNCTIONS FOR REVC ***************/
/******************************************************/

uint16_t computeDacVal_V(float voltage, int b, int c) {
  return uint16_t((65535.0 * (voltage - zeroPoint[b][c]) / 5.0));
}

uint16_t computeDacVal_I(float current, int b, int c) {
  return uint16_t((65535.0 * (current / gain[b][c] + 2.5 - zeroPoint[b][c]) / 5.0));
}

float computeOutI(uint16_t dacVal) {
  return ((float(dacVal) * 4.096 / 4096.0) - 1.25) / 10 / 0.2;
}

void print_all() { // Print all channel currents in A
  for (int i = 0; i < 7; i++) {
    uint16_t data = LTC1863ReadSlow(i);
    Serial.print(computeOutI(data), 5);
  }
  Serial.println();
}

float measure_gain(uint8_t b, uint8_t c) {
  //jump to 2.0 first so output returns nutral;
  selectBoard(b);
  delay(1);
  LTC2656Write(WRITE_AND_UPDATE, channelMap[c], computeDacVal_V(2.0, 0, 0));
  delayMicroseconds(1000);
  uint16_t out_2v0 = LTC1863ReadSlow(c, 50);
  //    Serial.println(computeOutI(out_2v0),5);
  LTC2656Write(WRITE_AND_UPDATE, channelMap[c], computeDacVal_V(2.5, 0, 0));
  delayMicroseconds(1000);
  uint16_t out_2v5 = LTC1863ReadSlow(c, 50);
  //    Serial.println(computeOutI(out_2v5),5);

  return (computeOutI(out_2v5) - computeOutI(out_2v0)) / (0.5);
}

bool calibrate_channel(uint8_t b, uint8_t c) {

  //  Serial.println();
  //  Serial.print("Calibrating board: "); Serial.print(b); Serial.print(" CH: "); Serial.print(c+1);
  bool isCalibrationSuccessful = false ;
  zeroPoint[b][c] = 0;
  delay(1);
  gain[b][c] = measure_gain(b, c);
  if (abs(gain[b][c] + 1.62) > 0.5) {
    //    Serial.println(" failed (gain)");
    calibrationStatus[b][c] = false;
    isCalibrationSuccessful = false ;
  } else {
  }
  //  Serial.print("gain: ");
  //  Serial.println(gain[b][c]);
  for (int i = 0; i < 10; i++) {
    float output_offset_I = computeOutI(LTC1863ReadSlow(c));
    //    Serial.print("iteration: ");
    //    Serial.println(i);
    //    Serial.println(output_offset_I,5);
    if (abs(output_offset_I) <= 0.001) {
      calibrationStatus[b][c] = true;
      isCalibrationSuccessful = true ;
    }
    zeroPoint[b][c] = zeroPoint[b][c] + (output_offset_I / gain[b][c]);
    //    Serial.print("next: ");
    //    Serial.println(zeroPoint[b][c],5);
    LTC2656Write(WRITE_AND_UPDATE, channelMap[c], computeDacVal_I(0, b, c));
    delay(10);
  }
  //  Serial.println("failed (cal)");
  calibrationStatus[b][c] = false;
  zeroPoint[b][c] = 0;
  LTC2656Write(WRITE_AND_UPDATE, channelMap[c], computeDacVal_I(0, b, c));
  isCalibrationSuccessful = false;
  return isCalibrationSuccessful ;
}

/******************************************************/
/******************* FUNCTIONS FROM ACDC **************/
/******************************************************/
void resetallshims( )
{
  for (int b = 0; b < NUM_B; b++) {
    for (int c = 0; c < NUM_C; c++) {
      selectBoard(b);
      LTC2656Write(WRITE_AND_UPDATE, channelMap[c], computeDacVal_I(0, b, c));
    }
  }
}

void usergetsystemheartbeat()
{
  //simple query to ensure system is responsive
  Serial.println( true );
}

//void usergetallchannelcurrents()
//{
//  for (int b = 0; b < NUM_B; b++) {
//    selectBoard(b);
//    print_all();
//  }
//}

bool userresetallshims()
{
  // same as resetallshims() but prints+returns true when finished
  resetallshims();
  Serial.println(true); return true;
}

//bool calibratedaccompensation()
//{
//
//  bool isCalibrationSuccessful = true ;
//  int SHIM_NCHANNELS = NUM_B * NUM_C;
//  bool isChannelCalibrationSuccessful [SHIM_NCHANNELS] ;
//  int iCh = 0 ;
//  for (int b = 0; b < NUM_B; b++) {
//    for (int c = 0; c < NUM_C; c++) {
//      selectBoard(b);
//      isChannelCalibrationSuccessful[iCh] = calibrate_channel(b, c);
//      iCh = iCh + 1;
//    }
//  }
//  for ( iCh = 0; iCh < SHIM_NCHANNELS; iCh++)
//    if ( !isChannelCalibrationSuccessful[iCh])
//      isCalibrationSuccessful = false ;
//
//  Serial.println(isCalibrationSuccessful);
//  return isCalibrationSuccessful ;
//}

bool readfivedigitcurrent( uint16_t &current )
{
  // Reads 5 digits from serial comprising a current value scaled between
  // [0:65535]
  //
  // Returns TRUE if successful

  String inString = "";
  uint8_t nBytesRead = 0;
  int inByte;
  long D[5] ; // data buffer

  while ( nBytesRead < 5 )
  {
    if ( Serial.available() > 0 )
    {
      inByte = Serial.read();

      if ( !isDigit( (char)inByte ) )
      {
        Serial.println(false); return false;
      }
      inString = (char)inByte ;
      D[nBytesRead] = inString.toInt() ;
      nBytesRead = nBytesRead + 1 ;
    }
  }

  current = uint16_t(D[0] * 10000 + D[1] * 1000 + D[2] * 100 + D[3] * 10 + D[4]) ;

  if ( ( (D[0] * 10000 + D[1] * 1000 + D[2] * 100 + D[3] * 10 + D[4]) < 0 ) ||
       ( (D[0] * 10000 + D[1] * 1000 + D[2] * 100 + D[3] * 10 + D[4]) > 65535 ) )
    return false;
  else
    return true;
}

uint16_t ampstodac( uint8_t iChannel, float current )
{
  return uint16_t( DAC_BITSPERVOLT * ( current * dacGain[iChannel] * float(DAC_PREAMP_RESISTANCE) + float( DAC_VREF ) + dacOffset[iChannel] ) / 1000.0  ) ;
}

void setshimbufferbychannel( uint8_t iCh, float current )
{
  currentsBuffer[iCh] = current ;
  dacBuffer[iCh] = ampstodac( iCh, currentsBuffer[iCh] ) ;
  selectBoard(board_order[iCh]);
  LTC2656Write(WRITE_TO_INPUT, channelMap[channel_order[iCh]], dacBuffer[iCh]);
}

float uint16toamps( uint16_t uint16current )
{
  // Input: current scaled between [0, 65535],
  // Output: current as float in units of amperes
  return ( float(uint16current) - 32768.0 ) / DAC_BITSPERAMP ;
}

bool usersetallshims()
{
  // Reads sequentially from serial 5 digits X SHIM_NCHANNELS
  // where each consecutive 5 digits represents the channel's
  // shim current scaled to be between [0:65535]
  //
  // Updates shim buffer upon completion

  bool isCurrentReadSuccessful = false;
  uint16_t inputCurrents [ NUM_B * NUM_C ] ;

  for ( uint8_t iCh = 0; iCh < (NUM_B * NUM_C); iCh++ )
  {
    isCurrentReadSuccessful = readfivedigitcurrent( inputCurrents[iCh] ) ;

    if ( !isCurrentReadSuccessful )
    {
      Serial.println(false); return false;
    }
  }

  for ( uint8_t iCh = 0; iCh < (NUM_B * NUM_C); iCh++ )
  {
    setshimbufferbychannel( iCh, uint16toamps( inputCurrents[iCh] ) ) ;
  }

  return isCurrentReadSuccessful ;
}

void loadallshims()
{
  // Write shim buffer entries of all channels to DAC
  // update DAC register for remaining channel + update all simultaneously

  // update DAC registers except for last channel
  for ( uint8_t iCh = 0; iCh < ((NUM_B * NUM_C) - 1); iCh++ ) {
    selectBoard(board_order[iCh]);
    LTC2656Write(WRITE_TO_INPUT, channelMap[channel_order[iCh]], dacBuffer[iCh]);
  }
  // update DAC register for remaining channel + update all simultaneously
  selectBoard(board_order[(NUM_B * NUM_C) - 1]);
  LTC2656Write(WRITE_TO_INPUT_UPDATE_ALL, channelMap[channel_order[(NUM_B * NUM_C) - 1]], dacBuffer[(NUM_B * NUM_C) - 1]);
}


bool usersetandloadallshims()
{
  // Reads sequentially from serial 5 digits X SHIM_NCHANNELS
  // where each consecutive 5 digits represents the channel's
  // shim current scaled to be between [0:65535]
  //
  // Updates the shim buffer + outputs DAC values upon completion
  //
  // Returns TRUE if successful

  bool isReadSuccessful = false ;
  isReadSuccessful = usersetallshims() ;
  loadallshims() ;
  Serial.println(isReadSuccessful);
  return isReadSuccessful;

}

uint16_t querychannelvoltage( uint8_t iChannel )
{ // return single channel voltage [units: mV]
  selectBoard(board_order[iChannel]);
  return ADC_MILLIVOLTSPERBIT * LTC1863ReadSlow(channelMap_ADC[channel_order[iChannel]]);
}

float querychannelcurrent( uint8_t iChannel )
{ // return single channel current [units: A]
  return ( float(querychannelvoltage( iChannel )) - float(DAC_VREF) ) / float(DAC_PREAMP_RESISTANCE) ;
}

void rampallshims()
{
  // Ramp all shims up to their buffered current values in 100 increments over 1.0 s

  float nCurrentDivisions = 100 ;
  float pause = 1.0 / nCurrentDivisions ; // max at 1.0 s

  float currents0[ NUM_B * NUM_C ] ;
  float currents1[ NUM_B * NUM_C ] ;
  float currentIncrements[ NUM_B * NUM_C ] ;

  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++ )
  {
    currents0[ iCh ]  = querychannelcurrent( iCh ) ;
    currents1[ iCh ]  = currentsBuffer[ iCh ] ;
    currentIncrements[ iCh ] = ( currents1[iCh] - currents0[iCh] ) / nCurrentDivisions ;
  }

  for ( uint8_t iCurrent = 1; iCurrent <= nCurrentDivisions; iCurrent++ )
  {
    for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++ )
    {
      setshimbufferbychannel( iCh, currents0[ iCh ] + iCurrent * currentIncrements[ iCh ] ) ;
    }
    loadallshims() ;
    delayMicroseconds(pause) ;
  }
}

void rampdownallshims( void )
{
  // Ramps down all channels to 0 A in increments over 1 s

  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++ )
  {
    currentsBuffer[iCh] = 0.0 ;
    dacBuffer[iCh] = ampstodac( iCh, currentsBuffer[iCh] ) ;
  }
  rampallshims();
}

void loadshimbychannel( uint8_t iCh )
{
  // Write shim buffer of single channel to DAC
  selectBoard(board_order[iCh]);
  LTC2656Write(WRITE_AND_UPDATE, channelMap[channel_order[iCh]], dacBuffer[iCh]);
}

void setandloadshimbychannel( uint8_t iCh, float current )
{
    setshimbufferbychannel( iCh, current ) ;
    loadshimbychannel( iCh ) ;
}

void rampshimbychannel( uint8_t iCh ) 
{
// Ramp single channel up to its buffered current value in 100 increments over 1.0 s
    
    float nCurrentDivisions = 100 ;
    float pause = 1.0/nCurrentDivisions ; // max at 1.0 s 
    
    float current0 = querychannelcurrent( iCh ) ;
    float current1 = currentsBuffer[ iCh ] ;
    
    float currentIncrement = ( current1 - current0 )/nCurrentDivisions ;
   
    for( uint8_t iCurrent = 1; iCurrent <= nCurrentDivisions; iCurrent++ )
    { 
        setandloadshimbychannel( iCh, current0 + iCurrent*currentIncrement ) ;
        delayMicroseconds(pause);
    }
}


bool calibratedaccompensation()
{
  // Determine the DAC voltage offsets and gain corrections for each channel
  //
  // Prints a bool for each channel, TRUE if calibration was succesful
  //
  // Returns TRUE if all channels successful
  bool isCalibrationSuccessful = false ;
  bool isChannelCalibrationSuccessful [NUM_B * NUM_C] ;

  float offsetError0 [NUM_B * NUM_C] ; // uncorrected
  float offsetError1 [NUM_B * NUM_C] ; // corrected
  float gainError0 [NUM_B * NUM_C] ; // uncorrected
  float gainError1 [NUM_B * NUM_C] ; // corrected

  float currentRequested;
  float currentsRead [ NUM_B * NUM_C ];
  float currentCorrected ;

  // reset DAC correction terms
  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++)
  {
    isChannelCalibrationSuccessful[iCh] = false ;
    dacGain[iCh]   = 1.0 ;
    dacOffset[iCh] = 0 ;
  }

  // determine DAC offset correction
  currentRequested = 0.0 ;
  // attempt to set all channels to 0.0 A
  resetallshims();

  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++)
  {
    int16_t voltageRead = querychannelvoltage( iCh ) ;
    currentsRead[iCh]   = querychannelcurrent( iCh ) ;
    dacOffset[iCh]      = voltageRead - DAC_VREF ;
  }

  // reset channels, now with the dac offsets adjusted
  resetallshims() ;

  // error of original vs. adjusted currents:
  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++)
  {
    currentCorrected = querychannelcurrent( iCh )  ;

    offsetError0[iCh] = abs( currentRequested - currentsRead[iCh] ) ;
    offsetError1[iCh] = abs( currentRequested - currentCorrected ) ;
  }

  // Determine DAC gain compensation
  currentRequested = 1.0 ;

  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++)
    setshimbufferbychannel( iCh, currentRequested ) ;

  rampallshims();

  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++)
  {
    currentsRead[iCh] = querychannelcurrent( iCh ) ;
    dacGain[iCh]      = currentRequested / currentsRead[iCh] ;
  }
  // update shims with dac correction in place & buffer still set at 1.0 A
  rampallshims();

  // pause to ensure query is accurate (some latency exists)
  delay(1);

  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++)
  {
    currentCorrected = querychannelcurrent( iCh ) ;

    // error of original vs. adjusted currents:
    gainError0[iCh] = abs( currentRequested - currentsRead[iCh] );
    gainError1[iCh] = abs( currentRequested - currentCorrected ) ;
  }

  rampdownallshims() ;

  // check adjusted results have lower error:
  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++)
  {
    // add 1.0 mV (~ADC precision) as tolerance for offset error
    if ( ( offsetError1[iCh] <= offsetError0[iCh] + 1.0 ) & ( gainError1[iCh] <= gainError0[iCh] ) )
    {
      isChannelCalibrationSuccessful[iCh] = true ;
      Serial.println( isChannelCalibrationSuccessful[iCh] ) ;
    }
    else
    {
      isChannelCalibrationSuccessful[iCh] = false ;
      Serial.println( isChannelCalibrationSuccessful[iCh] ) ;
    }
  }

  uint8_t iCh ;

  for ( iCh = 0; iCh < NUM_B * NUM_C; iCh++)
    if ( !isChannelCalibrationSuccessful[iCh] )
      break;

  if ( iCh == NUM_B * NUM_C )
    isCalibrationSuccessful = true;

  return isCalibrationSuccessful ;

}

bool usersetandrampallshims()
{
  // Reads sequentially from serial 5 digits X SHIM_NCHANNELS
  // where each consecutive 5 digits represents the channel's
  // shim current scaled to be between [0:65535]
  //
  // Updates the shim buffer and ramps current up over 1.0 s
  //
  // Returns TRUE if successful
  bool isReadSuccessful = false ;
  isReadSuccessful = usersetallshims() ;
  rampallshims() ;
  Serial.println(isReadSuccessful);  return isReadSuccessful;
}

bool usersetandloadshimbychannel( void )
{
  // Read from serial the channel index [0:SHIM_NCHANNELS] and 5 digit uint16
  // current val.  and set the single channel's current buffer
  //
  // Returns TRUE if successful

  String inString = "";
  uint8_t nBytesRead = 0;
  int inByte;
  uint8_t iCh = 0;
  long D ; // data buffer
  uint16_t current ;
  bool isCurrentReadSuccessful ;

  // read channel index
  while ( nBytesRead < 1 )
  {
    if ( Serial.available() > 0 )
    {
      inByte = Serial.read();

      if ( !isDigit( (char)inByte ) )
      {
        Serial.println(false); return false;
      }

      inString = (char)inByte ;
      D = inString.toInt() ;
      nBytesRead = nBytesRead + 1 ;

    }
  }

  iCh = uint8_t(D) ;
  if ( ( iCh < 0 ) || iCh >= NUM_B * NUM_C )
  {
    Serial.println(false); return false;
  }

  isCurrentReadSuccessful = readfivedigitcurrent( current ) ;

  if ( !isCurrentReadSuccessful )
  {
    Serial.println(false); return false ;
  }
  else
  {
    setandloadshimbychannel( iCh, uint16toamps( current ) ) ;
    Serial.println(true); return true;
  }

}

void usergetallchannelcurrents() 
{
// Print all channel currents in A
    for (uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++)
        Serial.println( querychannelcurrent( iCh ), 5 ) ;
}



bool usersetandrampshimbychannelasfloat( void ) 
{
// This function takes significantly longer than usersetandloadshimbychannel()
// and effectively produces the same result (setting the DAC buffer for a single shim
// channel and ramping up to it over 1.0 s) However, it may be more convenient
// for debugging:
//
// Rather than reading the channel current as a 5-digit scaled unsigned int,
// the input should be the channel index [0:SHIM_NCHANNELS-1] followed by the
// requested current, input as a single float in amperes.
//
// Returns TRUE if successful
    uint8_t iCh ;
    float current ;

    iCh = uint8_t( Serial.parseInt() ) ;

    if ( ( iCh < 0 ) || iCh >= NUM_B * NUM_C )
    {
        Serial.println(false); return false;
    }
    
    current = Serial.parseFloat() ; // [units: A]
    if ( abs(current) <= AMP_MAXCURRENTPERCHANNEL )
    {    
        setshimbufferbychannel( iCh, current ) ;
        rampshimbychannel( iCh ) ;
        Serial.println(true); return true;
    }
    else 
    {    
        Serial.println(false); return false;
    }    
    
}

bool userrampdownallshims() 
{
// same as resetallshims() but prints+returns true when finished 
// same as rampdownallshim() but prints+returns true when finished 
    rampdownallshims();
    Serial.println(true); return true;
}

void usergetallchannelvoltages( void ) 
{
// Print all channel voltages in mV
    for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++ )
        Serial.println( querychannelvoltage( iCh ) ) ;  
}

void usergetdaccompensationcoefficients( void ) 
{
    for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++ ) 
    {
        Serial.println( dacOffset[ iCh ], 5 ); 
        Serial.println( dacGain[ iCh ], 5 ); 
    }
}

