int ENA1=4; //define Enable pin
int ENA2=5; //define Enable pin

int PUL1=7; //define Pulse pin
int DIR1=6; //define Direction pin

int PUL2=9; //define Pulse pin
int DIR2=10; //define Direction pin

int TTL_S=12; //define Short TTL pin
int TTL_L=13; //define Long TTL pin

void ttl_short(){
    digitalWrite(TTL_S,HIGH);
    delayMicroseconds(1);
    digitalWrite(TTL_S,LOW);
 }
 
 
void ttl_long_on(){
    digitalWrite(TTL_L,HIGH);
 }
 
void ttl_long_off(){
    digitalWrite(TTL_L,LOW);
 }


void for_x(){          //Moves 1mm in the positive x direction
  for (int k=0; k<252; k++){
    digitalWrite(ENA1,HIGH);
    digitalWrite(DIR1,LOW);
    digitalWrite(PUL1,HIGH);
    delayMicroseconds(500);
    digitalWrite(PUL1,LOW);
    delayMicroseconds(500);
    digitalWrite(ENA1,LOW);
  }
  
}

void back_x(){          //Moves 1mm in the negative x direction
  for (int k=0; k<252; k++){
    digitalWrite(ENA1,HIGH);
    digitalWrite(DIR1,HIGH);
    digitalWrite(PUL1,HIGH);
    delayMicroseconds(500);
    digitalWrite(PUL1,LOW);
    delayMicroseconds(500);
    digitalWrite(ENA1,LOW);
  }
  
}

void for_y(){         //Moves 1mm in the positive y direction
  for (int k=0; k<252; k++){
    digitalWrite(ENA2,HIGH);
    digitalWrite(DIR2,LOW);
    digitalWrite(PUL2,HIGH);
    delayMicroseconds(500);
    digitalWrite(PUL2,LOW);
    delayMicroseconds(500);
    digitalWrite(ENA2,LOW);
  }
  
}

void back_y(){          //Moves 1mm in the negative y direction
  for (int k=0; k<252; k++){
    digitalWrite(ENA2,HIGH);
    digitalWrite(DIR2,HIGH);
    digitalWrite(PUL2,HIGH);
    delayMicroseconds(500);
    digitalWrite(PUL2,LOW);
    delayMicroseconds(500);
    digitalWrite(ENA2,LOW);
  }
  
}


void setup() {
  pinMode (PUL1, OUTPUT);
  pinMode (DIR1, OUTPUT);
  pinMode (PUL2, OUTPUT);
  pinMode (DIR2, OUTPUT);
  Serial.begin(9600);

  for (int i=0; i<10; i++)
  {
    for (int j=0; j<10; j++)
      
    {
      for_x();
      ttl_long_on();    //Start long TTL pulse (up)
      ttl_short();    //One 1us signal before starting acquisition
      
      delay(1000);    //Data acquisition time
      
      ttl_short();    //One 1us signal to finish acquisition
      ttl_long_off();   //Finish long TTL pulse (down)
      
    }
    
    for (int j=0; j<10; j++)
      
    {
      back_x();      
    }
    
    for_y();
    ttl_long_on();    //Start long TTL pulse (up)
    ttl_short();    //One 1us signal before starting acquisition
      
    delay(1000);    //Data acquisition time
      
    ttl_short();    //One 1us signal to finish acquisition
    ttl_long_off();   //Finish long TTL pulse (down)
    
    
  } 

    for (int j=0; j<10; j++)
      
    {
      back_y();      //Go back to starting point
    }
    

  

}

void loop() {

}
