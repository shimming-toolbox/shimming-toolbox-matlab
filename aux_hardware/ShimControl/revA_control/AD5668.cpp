/*
  AD5668.cpp, version 1.1, written by Robert Hart, October 15, 2015

  This library implements both hardware-SPI and software SPI constructors.

  This is an Arduino library for the Analog Devices AD5668 16-bit, 8 Channel
  Digital-to-Analog converter chip. It should also be usable with the Texas
  Instruments DAC8568, a pin and command compatible part. However, to maintain
  full compatibility with the AD5668, this library does not impliment the
  "Flexible Mode" internal reference command write sequences of the DAC8568
  (0b1001).

  Portions of this code were based upon prior work by crx and Stephan Bergemann,
  specifically in the structure of the writeDAC and init functions, as well as
  the use of a command definition list. This library impliments all of the
  command functions of the chip, and allows for tying the ~LDAC pin to ground.
  When tying ~LDAC to ground, DO NOT use the toggleLDAC function. Also, use of the
  toggleDAC function is probably only needed in conjuction with command function
  writeChannel (chip command 0). For a full explanation of the various functions,
  refer to the readme file and the examples.

  The AD5668 library is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the Free
  Software Foundation, version 3 of the License, which should be included in this
  this distribution. If not, see <http://www.gnu.org/licenses/>.
*/

#include "AD5668.h"
#include <SPI.h>

// Hardware SPI Constructor, must use hardware SPI bus pins
AD5668::AD5668(uint8_t ssPin, uint8_t clrPin, int8_t ldacPin) {
  _ssPin = ssPin;
  _clrPin = clrPin;
  _ldacPin = ldacPin;
  pinMode(_ssPin, OUTPUT);
  pinMode(_clrPin, OUTPUT);
  if (_ldacPin > 0) {
    pinMode(_ldacPin, OUTPUT);
  }
  _hwSPI = true;
}

// Software SPI Constructor, may use any digital pins
AD5668::AD5668(uint8_t mosiPin, uint8_t sclkPin, uint8_t ssPin, uint8_t clrPin, int8_t ldacPin) {
  _mosiPin = mosiPin;
  _sclkPin = sclkPin;
  _ssPin = ssPin;
  _clrPin = clrPin;
  _ldacPin = ldacPin;
  pinMode(_mosiPin, OUTPUT);
  pinMode(_sclkPin, OUTPUT);
  pinMode(_ssPin, OUTPUT);
  pinMode(_clrPin, OUTPUT);
  if (_ldacPin > 0) {
    pinMode(_ldacPin, OUTPUT);
  }
  _hwSPI = false;
}

void AD5668::init() {
  if (_hwSPI) {
    SPI.begin();
    SPI.setBitOrder(MSBFIRST);
  }
  else {
	digitalWrite(_mosiPin, LOW);
	digitalWrite(_sclkPin, LOW);
  }
  digitalWrite(_ssPin, LOW);
  if (_ldacPin > 0) {
    digitalWrite(_ldacPin, HIGH);
  }
  digitalWrite(_clrPin, LOW);
  delayMicroseconds(1);
  digitalWrite(_clrPin, HIGH);
  delayMicroseconds(1);
}

void AD5668::writeDAC(uint8_t command, uint8_t address, unsigned int data, uint8_t function) {
  byte b1 = command;
  byte b2 = address << 4 | data >> 12; //4 address bits and 4 MSBs of data
  byte b3 = data >> 4; // middle 8 bits of data
  byte b4 = 0xF0 & (data << 4) >> 8 | function;
  digitalWrite(_ssPin, LOW);
  delayMicroseconds(1);
  if (_hwSPI) {
    SPI.transfer(b1);
    SPI.transfer(b2);
    SPI.transfer(b3);
    SPI.transfer(b4);
  }
  else {
	shiftOut(_mosiPin, _sclkPin, MSBFIRST, b1);
    shiftOut(_mosiPin, _sclkPin, MSBFIRST, b2);
    shiftOut(_mosiPin, _sclkPin, MSBFIRST, b3);
    shiftOut(_mosiPin, _sclkPin, MSBFIRST, b4);
  }
  delayMicroseconds(1);
  digitalWrite(_ssPin, HIGH);
}

void AD5668::writeChannel(uint8_t channel, unsigned int value) {
  _channel = channel;
  _value = value;
  writeDAC(WRITE_INPUT_REGISTER, _channel, _value, 15);
}

void AD5668::updateChannel(uint8_t channel) {
  _channel = channel;
  writeDAC(UPDATE_OUTPUT_REGISTER, _channel, 0, 15);
}

void AD5668::writeChUpdateAll(uint8_t channel, unsigned int value) {
  _channel = channel;
  _value = value;
  writeDAC(WRITE_INPUT_REGISTER_UPDATE_ALL, _channel, _value, 15);
}

void AD5668::writeUpdateCh(uint8_t channel, unsigned int value) {
  _channel = channel;
  _value = value;
  writeDAC(WRITE_INPUT_REGISTER_UPDATE_N, _channel, _value, 15);
}

void AD5668::powerDAC_Normal(uint8_t channelsN) {
  _channelsN = channelsN;
  unsigned int modeChE_H = _channelsN >> 4;
  byte chA_D = 0x0F & _channelsN;
  writeDAC(POWER_DOWN_UP_DAC, 0, modeChE_H, chA_D);
}

void AD5668::powerDAC_Down1K(uint8_t channelsN) {
  _channelsN = channelsN;
  unsigned int modeChE_H = 0x10 | (_channelsN >> 4);
  byte chA_D = 0x0F & _channelsN;
  writeDAC(POWER_DOWN_UP_DAC, 0, modeChE_H, chA_D);
}

void AD5668::powerDAC_Down100K(uint8_t channelsN) {
  _channelsN = channelsN;
  unsigned int modeChE_H = 0x20 | (_channelsN >> 4);
  byte chA_D = 0x0F & _channelsN;
  writeDAC(POWER_DOWN_UP_DAC, 0, modeChE_H, chA_D);
}

void AD5668::powerDAC_DownTri(uint8_t channelsN) {
  _channelsN = channelsN;
  unsigned int modeChE_H = 0x30 | (_channelsN >> 4);
  byte chA_D = 0x0F & _channelsN;
  writeDAC(POWER_DOWN_UP_DAC, 0, modeChE_H, chA_D);
}

void AD5668::setClearCode(uint8_t clearCode) {
  _clearCode = clearCode;
  writeDAC(LOAD_CLEAR_CODE_REGISTER, 0, 0, _clearCode);
}

void AD5668::setSoftLDAC(uint8_t channelsN) {
  _channelsN = channelsN;
  unsigned int chE_H = _channelsN >> 4;
  byte chA_D = 0x0F & _channelsN;
  writeDAC(LOAD_LDAC_REGISTER, 0, chE_H, chA_D);
}

void AD5668::reset() {
  writeDAC(RESET_POWER_ON, 0, 0, 0);
}

void AD5668::enableInternalRef() {
  writeDAC(SETUP_INTERNAL_REF, 0, 0, 1);
}

void AD5668::disableInternalRef() {
  writeDAC(SETUP_INTERNAL_REF, 0, 0, 0);
}

void AD5668::toggleLDAC() {
  digitalWrite(_ldacPin, LOW);
  delayMicroseconds(1);
  digitalWrite(_ldacPin, HIGH);
  delayMicroseconds(1);
}
