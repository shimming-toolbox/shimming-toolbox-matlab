int startSerialComm() {
  //Initilize serial communication
  sampling_per = -1;
  while (1) {
    if (sampling_per <= 0)
    { sampling_per = Serial.parseInt();
      if (sampling_per > 0)
      { skip_cmpt = sampling_per / 10;
        Serial.println(skip_cmpt);
        break;
      }

    }
  }
}
