#include <math.h>
#include "MPU6050.h"

#define VERSION "NEBUCHADNEZZAR_202200304"
#define BAUD_RATE 115200
#define LED 13

MPU6050lib mpu; 

int16_t accelCount[3];
int16_t gyroCount[3];
float q[4] = {1.0f, 0.0f, 0.0f, 0.0f};
float dt = 0;
float gyroBias[3] = {0, 0, 0}, accelBias[3] = {0, 0, 0}; // Bias corrections for gyro and accelerometer
float GyroMeasError = PI * (40.0f / 180.0f);     // gyroscope measurement error in rads/s (start at 60 deg/s), then reduce after ~10 s to 3
float beta = sqrt(3.0f / 4.0f) * GyroMeasError;  // compute beta
float GyroMeasDrift = PI * (2.0f / 180.0f);      // gyroscope measurement drift in rad/s/s (start at 0.0 deg/s/s)
float zeta = sqrt(3.0f / 4.0f) * GyroMeasDrift;  // compute zeta, the other free parameter in the Madgwick scheme usually set to a small or zero value
float SelfTest[6];
uint32_t lastUpdate = 0, firstUpdate = 0;         // used to calculate integration interval
uint32_t Now = 0;                                 // used to calculate integration interval
float aRes, gRes; // scale resolutions per LSB for the sensors
float pitch, yaw, roll;

// Setup communication via USB/UART
String inputString = "";
boolean stringComplete = false;
struct command{
  String mnemonic;
  uint16_t option = 0;
  float value = 0.0f;
};

void setup() {
  pinMode(LED, OUTPUT);
  digitalWrite(LED, LOW);
  // Setup serial interrupt
  Serial.begin(BAUD_RATE);
  while(!Serial);
  inputString.reserve(200);

  // Self Test MPU650
  uint8_t addr = mpu.readByte(MPU6050_ADDRESS, WHO_AM_I_MPU6050);
  mpu.MPU6050SelfTest(SelfTest);
  if(addr == 0x68)
  {
    if (SelfTest[0] < 1.0f && SelfTest[1] < 1.0f && SelfTest[2] < 1.0f && SelfTest[3] < 1.0f && SelfTest[4] < 1.0f && SelfTest[5] < 1.0f)
    {
      mpu.calibrateMPU6050(gyroBias, accelBias);
      aRes = mpu.getAres();
      gRes = mpu.getGres()*PI/180.0f;
    }
    else
    {
        Serial.print("[!] ERROR: Could not connect to MPU6050: 0x");
        Serial.println(addr, HEX);
        while (1) ; // Loop forever if communication doesn't happen
      }
  }
}

void loop() {

  // If data ready bit set, all data registers have new data
  if (mpu.readByte(MPU6050_ADDRESS, INT_STATUS) & 0x01) { // check if data ready interrupt
    mpu.readAccelData(accelCount);  // Read the x/y/z adc values
    mpu.readGyroData(gyroCount); // Read the x/y/z adc values
    Now = micros();
    dt = ((Now - lastUpdate) / 1000000.0f); // set integration time by time elapsed since last filter update
    lastUpdate = Now;
    MadgwickQuaternionUpdate(q, accelCount, gyroCount);
  }

  // Check if command is completely read into the command structure
  if (stringComplete)
  {
    command cmd = parseData();

    if (cmd.mnemonic == "VER")
    {
      Serial.println(VERSION);
    }
    
    if (cmd.mnemonic == "YPR")
    {
      sendYPR();
    }
  }

  // Include this rather than the serialEvent() if deployed to Leonardo
  while (Serial.available())
  {
    char inChar = (char)Serial.read();
    if (inChar != '\n'  && inChar != '\r') {
      inputString += inChar;
    }
    if (inChar == '\n') {
      stringComplete = true;
    }
  }

}

command parseData()
{
  String string_copy = inputString;
  
  inputString = "";
  stringComplete = false;

  int str_len = string_copy.length() + 1;
  char char_array[str_len];
  string_copy.toCharArray(char_array, str_len);
  command cmd; 
  char *strtokIndx;

  strtokIndx = strtok(char_array,",");
  cmd.mnemonic = String(strtokIndx);

  strtokIndx = strtok(NULL, ",");
  cmd.option = atoi(strtokIndx);

  strtokIndx = strtok(NULL, ",");
  cmd.value = atof(strtokIndx);

  return cmd;
}

void sendYPR()
{
  calculateYPR();
  Serial.print(yaw);
  Serial.print(", ");
  Serial.print(pitch);
  Serial.print(", ");
  Serial.println(roll);  
}

void calculateYPR()
{
  yaw   = atan2(2.0f * (q[1] * q[2] + q[0] * q[3]), q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3]);
  pitch = -asin(2.0f * (q[1] * q[3] - q[0] * q[2]));
  roll  = atan2(2.0f * (q[0] * q[1] + q[2] * q[3]), q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3]);
  pitch *= 180.0f / PI;
  yaw   *= 180.0f / PI;
  roll  *= 180.0f / PI;
}

void MadgwickQuaternionUpdate(float *qarr, int16_t aCount[3], int16_t gCount[3])
{
    float ax = aCount[0]*aRes, ay = aCount[1]*aRes, az = aCount[2]*aRes;
    float gx = gCount[0]*gRes, gy = gCount[1]*gRes, gz = gCount[2]*gRes;
    float q1 = qarr[0], q2 = qarr[1], q3 = qarr[2], q4 = qarr[3];         // short name local variable for readability
    float norm;                                               // vector norm
    float f1, f2, f3;                                         // objetive funcyion elements
    float J_11or24, J_12or23, J_13or22, J_14or21, J_32, J_33; // objective function Jacobian elements
    float qDot1, qDot2, qDot3, qDot4;
    float hatDot1, hatDot2, hatDot3, hatDot4;
    float gerrx, gerry, gerrz, gbiasx, gbiasy, gbiasz;        // gyro bias error

    // Auxiliary variables to avoid repeated arithmetic
    float _halfq1 = 0.5f * q1;
    float _halfq2 = 0.5f * q2;
    float _halfq3 = 0.5f * q3;
    float _halfq4 = 0.5f * q4;
    float _2q1 = 2.0f * q1;
    float _2q2 = 2.0f * q2;
    float _2q3 = 2.0f * q3;
    float _2q4 = 2.0f * q4;
    float _2q1q3 = 2.0f * q1 * q3;
    float _2q3q4 = 2.0f * q3 * q4;

    // Normalise accelerometer measurement
    norm = sqrt(ax * ax + ay * ay + az * az);
    if (norm == 0.0f) return; // handle NaN
    norm = 1.0f/norm;
    ax *= norm;
    ay *= norm;
    az *= norm;
    
    // Compute the objective function and Jacobian
    f1 = _2q2 * q4 - _2q1 * q3 - ax;
    f2 = _2q1 * q2 + _2q3 * q4 - ay;
    f3 = 1.0f - _2q2 * q2 - _2q3 * q3 - az;
    J_11or24 = _2q3;
    J_12or23 = _2q4;
    J_13or22 = _2q1;
    J_14or21 = _2q2;
    J_32 = 2.0f * J_14or21;
    J_33 = 2.0f * J_11or24;
  
    // Compute the gradient (matrix multiplication)
    hatDot1 = J_14or21 * f2 - J_11or24 * f1;
    hatDot2 = J_12or23 * f1 + J_13or22 * f2 - J_32 * f3;
    hatDot3 = J_12or23 * f2 - J_33 *f3 - J_13or22 * f1;
    hatDot4 = J_14or21 * f1 + J_11or24 * f2;
    
    // Normalize the gradient
    norm = sqrt(hatDot1 * hatDot1 + hatDot2 * hatDot2 + hatDot3 * hatDot3 + hatDot4 * hatDot4);
    hatDot1 /= norm;
    hatDot2 /= norm;
    hatDot3 /= norm;
    hatDot4 /= norm;
    
    // Compute estimated gyroscope biases
    gerrx = _2q1 * hatDot2 - _2q2 * hatDot1 - _2q3 * hatDot4 + _2q4 * hatDot3;
    gerry = _2q1 * hatDot3 + _2q2 * hatDot4 - _2q3 * hatDot1 - _2q4 * hatDot2;
    gerrz = _2q1 * hatDot4 - _2q2 * hatDot3 + _2q3 * hatDot2 - _2q4 * hatDot1;
    
    // Compute and remove gyroscope biases
    gbiasx += gerrx * dt * zeta;
    gbiasy += gerry * dt * zeta;
    gbiasz += gerrz * dt * zeta;
    gx -= gbiasx;
    gy -= gbiasy;
    gz -= gbiasz;
    
    // Compute the quaternion derivative
    qDot1 = -_halfq2 * gx - _halfq3 * gy - _halfq4 * gz;
    qDot2 =  _halfq1 * gx + _halfq3 * gz - _halfq4 * gy;
    qDot3 =  _halfq1 * gy - _halfq2 * gz + _halfq4 * gx;
    qDot4 =  _halfq1 * gz + _halfq2 * gy - _halfq3 * gx;

    // Compute then integrate estimated quaternion derivative
    q1 += (qDot1 -(beta * hatDot1)) * dt;
    q2 += (qDot2 -(beta * hatDot2)) * dt;
    q3 += (qDot3 -(beta * hatDot3)) * dt;
    q4 += (qDot4 -(beta * hatDot4)) * dt;

    // Normalize the quaternion
    norm = sqrt(q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4);    // normalise quaternion
    norm = 1.0f/norm;
    q[0] = q1 * norm;
    q[1] = q2 * norm;
    q[2] = q3 * norm;
    q[3] = q4 * norm;
}
