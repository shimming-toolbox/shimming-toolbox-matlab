/*
  ShimControl.ino

// =========================================================================
// Updated::20181028::ryan.topfer@polymtl.ca
// =========================================================================
*/

#include "AD5668.h"
#include <SPI.h>
#include <Wire.h>
#include "Adafruit_ADS1015.h"
#include <stdlib.h>     /* atoi */
#include <math.h>

//Declaration of global variables 
Adafruit_ADS1015 adc1 (0x48);
Adafruit_ADS1015 adc2 (0x49); /* Use this for the 12-bit version */

#define mosiPin 11
#define sclkPin 13
#define ssPin 10
#define clrPin 8
#define ldacPin 7

/*
    Software SPI instance, "AD5668(mosiPin, sclkPin ,ssPin, clrPin, ldacPin);"
    mosiPin is connected to AD5668 Din pin (15) and sclkPin to AD5668 SCK
    pin (16). Remaining connections as explained above.
*/

// global variables 

// re: shim board
const uint8_t SHIM_NCHANNELS         = 8 ;

// power supply:
const float AMP_MAXCURRENTPERCHANNEL = 2.5 ; // [units: amps]
const float AMP_CURRENTRANGE         = 2*AMP_MAXCURRENTPERCHANNEL ; // [units: amps]

// 3 terms describing ADC hardware:
const uint8_t ADC_NCHANNELS          = 4 ;
const uint8_t ADC_RESOLUTION         = 12 ; // 12-bit
const float   ADC_RANGE_VOUT         = 2.048 ; // for GAIN_TWO setting [units: volts]
// 1 derived term for convenience:
const uint16_t ADC_MILLIVOLTSPERBIT  = round( 1000*2.0*ADC_RANGE_VOUT/( pow( 2.0, float(ADC_RESOLUTION) ) - 1.0 ) ) ; // = 1 [uints: mV/bit-count]

// 3 terms describing DAC hardware:
AD5668 DAC                           = AD5668( mosiPin, sclkPin , ssPin, clrPin, ldacPin ) ;
const uint8_t DAC_RESOLUTION         = 16 ; // 16-bit
const int16_t DAC_VREF               = 1250 ; // [units: mV]
const int16_t DAC_PREAMP_RESISTANCE  = 220 ; // [units: milli-Ohms]

// 3 derived terms for convenience:
const uint16_t DAC_RANGE_VOUT        = 2*DAC_VREF ; // [units: mV]
const float DAC_BITSPERAMP           = ( pow( 2.0, float(DAC_RESOLUTION) ) - 1.0 )/AMP_CURRENTRANGE ; // = 13107 [units: bit-counts/A] 
const uint16_t DAC_BITSPERVOLT       = (1000* pow( 2.0, float(DAC_RESOLUTION) ) - 1.0 )/DAC_RANGE_VOUT ; // = 26214 [units: bit-counts/V]

// 
float dacOffset [ SHIM_NCHANNELS ] ; // [units: mV]
float dacGain [ SHIM_NCHANNELS ] ; // [unitless]

float currentsBuffer [ SHIM_NCHANNELS ] ; // [units: A] 
uint16_t dacBuffer [ SHIM_NCHANNELS ] ; // same as currentsBuffer but converted to DAC counts

/* bool isPrintModeVerbose = true ; */

void setup() 
{
  Serial.begin(115200);   //Baudrate of the serial communication : Maximum
  delay(100);

  DAC.init();
  DAC.enableInternalRef(); 
  DAC.powerDAC_Normal(B11111111); // Power up all channels normal

  adc1.setGain(GAIN_TWO); //+/- 2.048V  1 bit = 1mV
  adc2.setGain(GAIN_TWO); //+/- 2.048V  1 bit = 1mV

  adc1.begin();
  adc2.begin();

  for( uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++ )
  {
    currentsBuffer[iCh] = 0.0 ;
    dacBuffer[iCh]      = 32768 ;
    dacOffset[iCh]      = 0 ; 
    dacGain[iCh]        = 1.0 ; 
  }

  resetallshims() ; 
  // system heartbeat prints TRUE to indicate system is responsive
  usergetsystemheartbeat() ; 
  // TODO : print FALSE if setup unexpectedly failed  
}

void(* resetFunc) (void) = 0; //declare reset function at address 0

void loop() 
{
  
  char incomingByte;

  if (Serial.available() > 0) 
    incomingByte = Serial.read();
  
  switch (incomingByte) 
  {
    
    case 'a':
      usersetandloadallshims(); 
      break;
    
    case 'c':                 
      calibratedaccompensation() ;
      break;
    
    case 'f': // read in channel index followed by channel current as float in units of amperes
      usersetandloadshimbychannelasfloat();
      break;
    
    case 'h':
      usergetsystemheartbeat(); 
      break;
    
    case 'i':
      usersetandloadshimbychannel(); 
      break;
    
    /* case 'm': */
    /* //Toggle between verbose print mode (true/false) */
    /*   userswitchprintmode(); */
    /*   break; */
    /*  */
    /* case 'n':                 */
    /* //Return current print mode (true when verbose) */
    /*   usergetprintmode(); */
    /*   break; */

    case 'q':                
      usergetallchannelcurrents();
      break;
    
    case 'r':                
      userresetallshims( ) ;
      break;
    
    case 'v':                
      usergetallchannelvoltages();
      break;

    case 'u':                 
        usergetcompensationcoefficients() ;
      break;
    
    case 'z':                // Reset the arduino
      Serial.println("\n");
      resetFunc();
      break;

  }
}


uint16_t ampstodac( uint8_t iChannel, float current ) 
{
    return uint16_t( DAC_BITSPERVOLT * ( current*dacGain[iChannel]*float(DAC_PREAMP_RESISTANCE) + float( DAC_VREF ) + dacOffset[iChannel] )/1000.0  ) ;
}

float uint16toamps( uint16_t uint16current ) 
{
// Input: current scaled between [0, 65535], 
// Output: current as float in units of amperes 
    return ( float(uint16current) - 32768.0 )/DAC_BITSPERAMP ;
}

bool calibratedaccompensation( void )
{
// Determine the DAC voltage offsets and gain corrections for each channel
// 
// Prints a bool for each channel, TRUE if calibration was succesful
//
// Returns TRUE if all channels successful
  bool isCalibrationSuccessful = false ;
  bool isChannelCalibrationSuccessful [SHIM_NCHANNELS] ;

  float offsetError0 [SHIM_NCHANNELS] ; // uncorrected
  float offsetError1 [SHIM_NCHANNELS] ; // corrected
  float gainError0 [SHIM_NCHANNELS] ; // uncorrected
  float gainError1 [SHIM_NCHANNELS] ; // corrected

  float currentRequested;
  float currentsRead [ SHIM_NCHANNELS ];   
  float currentCorrected ;

  // reset DAC correction terms 
  for ( uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++)  
  {
    isChannelCalibrationSuccessful[iCh] = false ;
    dacGain[iCh]   = 1.0 ;
    dacOffset[iCh] = 0 ; 
  }

  // determine DAC offset correction
  currentRequested = 0.0 ;
  // attempt to set all channels to 0.0 A
  resetallshims(); 
  
  for ( uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++)  
  {
    int16_t voltageRead = querychannelvoltage( iCh ) ;
    currentsRead[iCh]   = querychannelcurrent( iCh ) ;
    dacOffset[iCh]      = voltageRead - DAC_VREF ;  
  }

  // reset channels, now with the dac offsets adjusted
  resetallshims() ; 
  
  // error of original vs. adjusted currents: 
  for ( uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++)  
  {
    currentCorrected = querychannelcurrent( iCh )  ;

    offsetError0[iCh] = abs( currentRequested - currentsRead[iCh] ) ;
    offsetError1[iCh] = abs( currentRequested - currentCorrected ) ;
  }
  
  // Determine dac gain compensation 
  currentRequested = 1.0 ;
  
  for ( uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++)  
    setandloadshimbychannel( iCh, currentRequested ) ;
  
  for ( uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++)  
  {
    currentsRead[iCh] = querychannelcurrent( iCh ) ;
    dacGain[iCh]      = currentRequested/currentsRead[iCh] ; 

    // update shim, now with adjusted correction:
    setandloadshimbychannel( iCh, currentRequested ) ;

    // pause to ensure query is accurate (some latency exists)
    delay(1000);
    currentCorrected = querychannelcurrent( iCh ) ;

    // error of original vs. adjusted currents: 
    gainError0[iCh] = abs( currentRequested - currentsRead[iCh] );
    gainError1[iCh] = abs( currentRequested - currentCorrected ) ;
    // NOTE:
    // Given the necessary delay between setting/reading the channels
    // it would be faster to load all channels at once, however that would
    // result in a much larger demand on the power supply...  
    // Hence, calibrating channel-by-channel for now:
    setandloadshimbychannel( iCh, 0.0 ) ;
  }
  
  // check adjusted results have lower error: 
  for ( uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++)  
  {
    if( ( offsetError1[iCh] < offsetError0[iCh] ) & ( gainError1[iCh] < gainError0[iCh] ) )
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

    for( iCh = 0; iCh < SHIM_NCHANNELS; iCh++)
        if ( !isChannelCalibrationSuccessful[iCh] )
            break;

    if ( iCh == SHIM_NCHANNELS )
        isCalibrationSuccessful = true;
      
    return isCalibrationSuccessful ; 
        
}

void loadallshims( void ) 
{
// Write shim buffer entries of all channels to DAC
    for ( uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++ ) 
        DAC.writeUpdateCh( iCh, dacBuffer[iCh] );    
}

void loadshimbychannel( uint8_t iCh ) 
{
// Write shim buffer of single channel to DAC
    DAC.writeUpdateCh( iCh, dacBuffer[iCh] );    
}


float querychannelcurrent( uint8_t iChannel ) 
{// return single channel current [units: A]
    return ( float(querychannelvoltage( iChannel )) - float(DAC_VREF) )/float(DAC_PREAMP_RESISTANCE) ;
}

uint16_t querychannelvoltage( uint8_t iChannel ) 
{// return single channel voltage [units: mV]
    return ADC_MILLIVOLTSPERBIT * readchanneladc( iChannel ) ; 
}

uint16_t readchanneladc( uint8_t iChannel ) 
{
    if ( iChannel < ADC_NCHANNELS )                                  
        return adc1.readADC_SingleEnded( iChannel ) ;
    else if ( iChannel < 2*ADC_NCHANNELS )
        return adc2.readADC_SingleEnded( iChannel - ADC_NCHANNELS ) ;
    else
		Serial.println("Error: Invalid channel index. ") ;
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

void resetallshims( void ) 
{
    for ( uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++ ) 
    {
        currentsBuffer[iCh] = 0.0 ;
        dacBuffer[iCh] = ampstodac( iCh, currentsBuffer[iCh] ) ;
    }
    loadallshims();
}


void setandloadshimbychannel( uint8_t iCh, float current )
{
    setshimbufferbychannel( iCh, current ) ;
    loadshimbychannel( iCh ) ;
}

void setshimbufferbychannel( uint8_t iCh, float current )
{
    currentsBuffer[iCh] = current ;
    dacBuffer[iCh] = ampstodac( iCh, currentsBuffer[iCh] ) ;
}

// =========================================================================
// User callable functions
// =========================================================================

void usergetallchannelcurrents() 
{
// Print all channel currents in A
    for (uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++)
        Serial.println( querychannelcurrent( iCh ), 5 ) ;
}

void usergetallchannelvoltages( void ) 
{
// Print all channel voltages in mV
    for ( uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++ )
        Serial.println( querychannelvoltage( iCh ) ) ;  
}

/* void userswitchprintmode( void )  */
/* { */
/*      isPrintModeVerbose = !isPrintModeVerbose ;  */
/* } */
/*  */
/* void usergetprintmode( void )  */
/* { */
/*      Serial.println( isPrintModeVerbose ) ;  */
/* } */

void usergetsystemheartbeat( void ) 
{
//simple query to ensure system is responsive
    Serial.println( true ); 
}

void usergetcompensationcoefficients( void ) 
{
    for ( uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++ ) 
    {
        Serial.println( dacOffset[ iCh ] ); 
        Serial.println( dacGain[ iCh ] ); 
    }
}

bool userresetallshims( void ) 
{
// same as resetallshims() but prints+returns true when finished 
    resetallshims();
    Serial.println(true); return true;
}

bool usersetandloadallshims( void ) 
{
// Reads sequentially from serial 5 digits X SHIM_NCHANNELS 
// where each consecutive 5 digits represents the channel's
// shim current scaled to be between [0:65535]
// 
// Updates the shim buffer + DAC values upon completion
//
// Returns TRUE if successful

    bool isCurrentReadSuccessful = false;
    uint16_t inputCurrents [ SHIM_NCHANNELS ] ;

    for( uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++ )
    {
        isCurrentReadSuccessful = readfivedigitcurrent( inputCurrents[iCh] ) ;

        if ( !isCurrentReadSuccessful )
        {
            Serial.println(false); return false; 
        }
    }

    for( uint8_t iCh = 0; iCh < SHIM_NCHANNELS; iCh++ )
    {
        setshimbufferbychannel( iCh, uint16toamps( inputCurrents[iCh] ) ) ;
    }    
    
    loadallshims() ;
    Serial.println(true);  return true;
        
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
    while( nBytesRead < 1 )
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
    if ( ( iCh < 0 ) || iCh >= SHIM_NCHANNELS )
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
    
bool usersetandloadshimbychannelasfloat( void ) 
{
// This function takes significantly longer than usersetandloadshimbychannel() 
// and produces the same result (setting and loading a single shim channel).
// However, it may be more convenient for debugging:
//
// Rather than reading the channel current as a 5-digit scaled unsigned int, 
// the input should be the channel index [0:SHIM_NCHANNELS-1] followed by the 
// requested current, input as a single float in amperes.
//
// Returns TRUE if successful
    uint8_t iCh ;
    float current ;

    iCh = uint8_t( Serial.parseInt() ) ;

    if ( ( iCh < 0 ) || iCh >= SHIM_NCHANNELS )
    {
        Serial.println(false); return false;
    }
    
    current = Serial.parseFloat() ; // [units: A]
    if ( abs(current) <= AMP_MAXCURRENTPERCHANNEL )
    {    
        setandloadshimbychannel( iCh, current ) ;
        Serial.println(true); return true;
    }
    else 
    {    
        Serial.println(false); return false;
    }    
    
}
