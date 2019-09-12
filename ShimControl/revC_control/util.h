
/******************************************************/
/******************* FUNCTIONS FROM ACDC **************/
/******************************************************/
uint16_t ampstodac( uint8_t iChannel, float current )
{
  return uint16_t( DAC_BITSPERVOLT * ( current * dacGain[iChannel] * float(DAC_PREAMP_RESISTANCE) + float( DAC_VREF ) + dacOffset[iChannel] ) / 1000.0  ) ;
}

void resetallshims( )
{
  for (int b = 0; b < NUM_B; b++) {
    for (int c = 0; c < NUM_C; c++) {
      selectBoard(b);
      LTC2656Write(WRITE_AND_UPDATE, channelMap[c], ampstodac(c, 0));
    }
  }
}

void usergetsystemheartbeat()
{
  //simple query to ensure system is responsive
  Serial.println( true );
}

bool userresetallshims()
{
  // same as resetallshims() but prints+returns true when finished
  resetallshims();
  Serial.println(true); return true;
}

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

  return ADC_MILLIVOLTSPERBIT * LTC1863ReadSlow(channelMap_ADC[channel_order[iChannel]]) / float(R_L);
}

float querychannelcurrent( uint8_t iChannel )

{ 
  // return single channel current [units: A]
//  const float   R_L                    = 2.00;
//  const float   R_OS                   = 10.00;
  if (iChannel == 7) {
//    Serial.println(R_L);
    return ( float(querychannelvoltage( iChannel ))  - float(DAC_VREF)) / float(DAC_PREAMP_RESISTANCE) / float(R_L) ;
  }
  else
  {
//    Serial.println(R_OS);
    return ( float(querychannelvoltage( iChannel ))  - float(DAC_VREF)) / float(DAC_PREAMP_RESISTANCE) / float(R_OS);
  }
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

bool calibratedaccompensation()
{
  //   Determine the DAC voltage offsets and gain corrections for each channel
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
  delay(1000);

  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++)
  {
    int16_t voltageRead = querychannelvoltage( iCh ) ;
    currentsRead[iCh]   = querychannelcurrent( iCh ) ;
    dacOffset[iCh]      = voltageRead - DAC_VREF ;

  }

  // reset channels, now with the dac offsets adjusted
  resetallshims() ;
  delay(1000);


  // error of original vs. adjusted currents:
  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++)
  {
    currentCorrected = querychannelcurrent( iCh )  ;
    offsetError0[iCh] = abs( currentRequested - currentsRead[iCh] ) ;
    offsetError1[iCh] = abs( currentRequested - currentCorrected ) ;
  }

  // Determine DAC gain compensation
  currentRequested = 0.5 ;

  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++)
  { LTC2656Write(WRITE_AND_UPDATE, channelMap[channel_order[iCh]], ampstodac(iCh, currentRequested)) ;
  }
  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++)
  {
    currentsRead[iCh] = querychannelcurrent( iCh ) ;
    dacGain[iCh]      = currentRequested / currentsRead[iCh] ;
    LTC2656Write(WRITE_AND_UPDATE, channelMap[channel_order[iCh]], ampstodac(iCh, currentRequested)) ;
    // update shims with dac correction in place & buffer still set at 1.0 A
  }


  resetallshims() ;
  // pause to ensure query is accurate (some latency exists)
  delay(1000);

  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++)
  {

    currentCorrected = querychannelcurrent( iCh ) ;
    // error of original vs. adjusted currents:
    gainError0[iCh] = abs( currentRequested - currentsRead[iCh] );
    gainError1[iCh] = abs( currentRequested - currentCorrected ) ;
  }

  resetallshims();
  delay(1000);

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
//  Serial.println("-------------");
  for (uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++)
    Serial.println( querychannelcurrent( iCh ), 5 ) ;
//  Serial.println("-------------");
}

void usergetallchannelvoltages( void )
{
  // Print all channel voltages in mV
  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++ )
    Serial.println( querychannelvoltage( iCh ) ) ;
}

bool usersetandloadshimbychannelasfloat(void)
{ selectBoard(0);
  uint8_t iCh ;
  float current ;

  iCh = uint8_t( Serial.parseInt() );
  Serial.print("CH: "); Serial.println(iCh);
  iCh = iCh - 1;

  if ( ( iCh < 0 ) || iCh >= NUM_B * NUM_C )
  {
    Serial.println(false); return false;
  }

  current = Serial.parseFloat() ; // [units: A]
  Serial.print("Current: "); Serial.println(current);
  if ( abs(current) <= AMP_MAXCURRENTPERCHANNEL )
  {
    LTC2656Write(WRITE_AND_UPDATE, channelMap[channel_order[iCh]], ampstodac(iCh, current)) ;
    Serial.println(true); return true;
  }
  else
  {
    Serial.println(false); return false;
  }

}

void usergetdaccompensationcoefficients( void )
{
  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++ )
  {
    Serial.println( dacOffset[ iCh ], 5 );
    Serial.println( dacGain[ iCh ], 5 );
  }
}

//uint16_t ampstodac2( uint8_t iChannel, float current , float dacGain, float dacOffset)
//{
//  return uint16_t( DAC_BITSPERVOLT * ( current * dacGain * float(DAC_PREAMP_RESISTANCE) + float( DAC_VREF ) + dacOffset ) / 1000.0  ) ;
//}
//
//void calibratedaccompensation2(uint8_t iCh)
//{
//  // Determine the DAC voltage offsets and gain corrections for each channel
//  //
//  // Prints a bool for each channel, TRUE if calibration was succesful
//  //
//  // Returns TRUE if all channels successful
//  //bool isCalibrationSuccessful = false ;
//  //bool isChannelCalibrationSuccessful [NUM_B * NUM_C] ;
//  selectBoard(0);
//  float offsetError0 ; // uncorrected
//  float offsetError1 ; // corrected
//  float gainError0  ; // uncorrected
//  float gainError1  ; // corrected
//
//  float currentRequested;
//  float currentsRead ;
//  float currentCorrected ;
//
//  // reset DAC correction terms
//  //  {
//  //    isChannelCalibrationSuccessful[iCh] = false ;
//  //    dacGain[iCh]   = 1.0 ;
//  //    dacOffset[iCsh] = 0 ;
//  //  }
//  float dacGain1   = 1.0 ;
//  float dacOffset1 = 0 ;
//  // determine DAC offset correction
//  currentRequested = 0.0 ;
//  // attempt to set all channels to 0.0 A
//  resetallshims();
//
//  int16_t voltageRead = querychannelvoltage( iCh ) ;
//  currentsRead   = querychannelcurrent( iCh ) ;
//  dacOffset1      = voltageRead - DAC_VREF ;
//
//  Serial.print("DAC offset: "); Serial.println(dacOffset1);
//
//  // reset channels, now with the dac offsets adjusted
//  resetallshims() ;
//
//  // error of original vs. adjusted currents:
//  iCh = 7;
//
//  currentCorrected = querychannelcurrent( iCh )  ;
//
//  offsetError0 = abs( currentRequested - currentsRead ) ;
//  offsetError1 = abs( currentRequested - currentCorrected ) ;
//
//  Serial.print("offsetError0: "); Serial.println(offsetError0);
//  Serial.print("offsetError1: "); Serial.println(offsetError1);
//
//  //  // Determine DAC gain compensation
//  currentRequested = 1.0 ;
//
//  LTC2656Write(WRITE_AND_UPDATE, channelMap[channel_order[iCh]], ampstodac2(iCh, currentRequested, dacGain1, dacOffset1)) ;
//  //
//  currentsRead = querychannelcurrent( iCh ) ;
//  dacGain1       = currentRequested / currentsRead ;
//
//
//
//  Serial.print("dacGain1: "); Serial.println(dacGain1);
//  LTC2656Write(WRITE_AND_UPDATE, channelMap[channel_order[iCh]], ampstodac2(iCh, currentRequested, dacGain1, dacOffset1)) ;
//
//  currentCorrected = querychannelcurrent( iCh ) ;
//
//  // error of original vs. adjusted currents:
//  gainError0 = abs( currentRequested - currentsRead );
//  gainError1 = abs( currentRequested - currentCorrected ) ;
//
//  Serial.print("gainError0: "); Serial.println(gainError0);
//  Serial.print("gainError1: "); Serial.println(gainError1);
//  Serial.print("1 amp query:");
//  usergetallchannelcurrents();
//  delay(500);
//  rampdownallshims() ;
//  delay(500);
//  Serial.println();
//  LTC2656Write(WRITE_AND_UPDATE, channelMap[channel_order[iCh]], ampstodac2(iCh, -1, dacGain1, dacOffset1));
//  Serial.println("-1 amp query:");
//  usergetallchannelcurrents();
//  delay(500);
//  rampdownallshims() ;
//  delay(500);
//  Serial.println();
//  LTC2656Write(WRITE_AND_UPDATE, channelMap[channel_order[iCh]], ampstodac2(iCh, 2, dacGain1, dacOffset1));
//  Serial.println("2 amp query:");
//  usergetallchannelcurrents();
//  delay(500);
//  rampdownallshims() ;
//  delay(500);
//  Serial.println();
//  LTC2656Write(WRITE_AND_UPDATE, channelMap[channel_order[iCh]], ampstodac2(iCh, -2, dacGain1, dacOffset1));
//  Serial.println("-2 amp query:");
//  usergetallchannelcurrents();
//  delay(500);
//  rampdownallshims() ;
//  delay(500);
//  Serial.println();
//  LTC2656Write(WRITE_AND_UPDATE, channelMap[channel_order[iCh]], ampstodac2(iCh, 0, dacGain1, dacOffset1));
//  Serial.println("0 amp query:");
//  usergetallchannelcurrents();
//  delay(500);
//  rampdownallshims() ;
//  delay(500);
//  Serial.println();
//
//  rampdownallshims() ;
//
//
//  //  }
//  //  //  // update shims with dac correction in place & buffer still set at 1.0 A
//  //  //  rampallshims();
//  //  //
//  //  //  // pause to ensure query is accurate (some latency exists)
//  //  delay(1);
//  //  //
//  //  //  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++)
//  //  {
//  //    currentCorrected = querychannelcurrent( iCh ) ;
//  //
//  //    // error of original vs. adjusted currents:
//  //    gainError0[iCh] = abs( currentRequested - currentsRead[iCh] );
//  //    gainError1[iCh] = abs( currentRequested - currentCorrected ) ;
//  //  }
//  //  Serial.println(gainError0[iCh]);
//  //  Serial.println(gainError1[iCh]);
//  //  //
//  //  rampdownallshims() ;
//  //  //
//  //  //  // check adjusted results have lower error:
//  //  //  for ( uint8_t iCh = 0; iCh < NUM_B * NUM_C; iCh++)
//  //  //  {
//  //  //    // add 1.0 mV (~ADC precision) as tolerance for offset error
//  //  //    if ( ( offsetError1[iCh] <= offsetError0[iCh] + 1.0 ) & ( gainError1[iCh] <= gainError0[iCh] ) )
//  //  //    {
//  //  //      isChannelCalibrationSuccessful[iCh] = true ;
//  //  //      Serial.println( isChannelCalibrationSuccessful[iCh] ) ;
//  //  //    }
//  //  //    else
//  //  //    {
//  //  //      isChannelCalibrationSuccessful[iCh] = false ;
//  //  //      Serial.println( isChannelCalibrationSuccessful[iCh] ) ;
//  //  //    }
//  //  //  }
//  //  //
//  //  //  uint8_t iCh ;
//  //  //
//  //  //  for ( iCh = 0; iCh < NUM_B * NUM_C; iCh++)
//  //  //    if ( !isChannelCalibrationSuccessful[iCh] )
//  //  //      break;
//  //  //
//  //  //  if ( iCh == NUM_B * NUM_C )
//  //  //    isCalibrationSuccessful = true;
//  //  //
//  //  //  return isCalibrationSuccessful ;
//  //
//}
