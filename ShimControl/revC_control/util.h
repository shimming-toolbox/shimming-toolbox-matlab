

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

void usergetallchannelcurrents()
{
  for (int b = 0; b < NUM_B; b++) {
    selectBoard(b);
    print_all();
  }
}

bool userresetallshims()
{
  // same as resetallshims() but prints+returns true when finished
  resetallshims();
  Serial.println(true); return true;
}

bool calibratedaccompensation()
{

  bool isCalibrationSuccessful = true ;
  int SHIM_NCHANNELS = NUM_B * NUM_C;
  bool isChannelCalibrationSuccessful [SHIM_NCHANNELS] ;
  int iCh = 0 ;
  for (int b = 0; b < NUM_B; b++) {
    for (int c = 0; c < NUM_C; c++) {
      selectBoard(b);
      isChannelCalibrationSuccessful[iCh] = calibrate_channel(b, c);
      iCh = iCh + 1;
    }
  }
  for ( iCh = 0; iCh < SHIM_NCHANNELS; iCh++)
    if ( !isChannelCalibrationSuccessful[iCh])
      isCalibrationSuccessful = false ;

  Serial.println(isCalibrationSuccessful);
  return isCalibrationSuccessful ;
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
    
    while( nBytesRead < 5 )
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
    
    current = uint16_t(D[0]*10000 + D[1]*1000 + D[2]*100 + D[3]*10 + D[4]) ;
    
    if ( ( (D[0]*10000 + D[1]*1000 + D[2]*100 + D[3]*10 + D[4]) < 0 ) || 
            ( (D[0]*10000 + D[1]*1000 + D[2]*100 + D[3]*10 + D[4]) > 65535 ) )
        return false;
    else
        return true;
}

uint16_t ampstodac( uint8_t iChannel, float current ) 
{
    return uint16_t( DAC_BITSPERVOLT * ( current*dacGain[iChannel]*float(DAC_PREAMP_RESISTANCE) + float( DAC_VREF ) + dacOffset[iChannel] )/1000.0  ) ;
}

void setshimbufferbychannel( uint8_t iCh, float current )
{
    currentsBuffer[iCh] = current ;
    dacBuffer[iCh] = ampstodac( iCh, currentsBuffer[iCh] ) ;
    selectBoard(board_order[iCh]);
    LTC2656Write(WRITE_AND_UPDATE, channelMap[channel_order[iCh]],dacBuffer[iCh]);
}

float uint16toamps( uint16_t uint16current ) 
{
// Input: current scaled between [0, 65535], 
// Output: current as float in units of amperes 
    return ( float(uint16current) - 32768.0 )/DAC_BITSPERAMP ;
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

    for( uint8_t iCh = 0; iCh < (NUM_B * NUM_C); iCh++ )
    {
        isCurrentReadSuccessful = readfivedigitcurrent( inputCurrents[iCh] ) ;

        if ( !isCurrentReadSuccessful )
        {
            Serial.println(false); return false; 
        }
    }

    for( uint8_t iCh = 0; iCh < (NUM_B * NUM_C); iCh++ )
    {
        setshimbufferbychannel( iCh, uint16toamps( inputCurrents[iCh] ) ) ;
    }   

    return isCurrentReadSuccessful ; 
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
    //loadallshims() ;
    Serial.println(isReadSuccessful);
    return isReadSuccessful;
        
}

