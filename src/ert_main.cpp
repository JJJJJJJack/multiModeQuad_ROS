//
// File: ert_main.cpp
//
// Code generated for Simulink model 'multiModeQuad_ROS'.
//
// Model version                  : 1.71
// Simulink Coder version         : 9.6 (R2021b) 14-May-2021
// C/C++ source code generated on : Thu Mar 20 00:46:27 2025
//
// Target selection: ert.tlc
// Embedded hardware selection: Generic->Unspecified (assume 32-bit Generic)
// Code generation objectives: Unspecified
// Validation result: Not run
//
#include <stdio.h>
#include <stdlib.h>
#include "multiModeQuad_ROS.h"
#include "multiModeQuad_ROS_private.h"
#include "rtwtypes.h"
#include "limits.h"
#include "rt_nonfinite.h"
#include "linuxinitialize.h"
#define UNUSED(x)                      x = x
#define NAMELEN                        16

// Function prototype declaration
void exitFcn(int sig);
void *terminateTask(void *arg);
void *baseRateTask(void *arg);
void *subrateTask(void *arg);
volatile boolean_T stopRequested = false;
volatile boolean_T runModel = true;
sem_t stopSem;
sem_t baserateTaskSem;
pthread_t schedulerThread;
pthread_t baseRateThread;
void *threadJoinStatus;
int terminatingmodel = 0;
void *baseRateTask(void *arg)
{
  runModel = (rtmGetErrorStatus(multiModeQuad_ROS_M) == (NULL)) &&
    !rtmGetStopRequested(multiModeQuad_ROS_M);
  while (runModel) {
    sem_wait(&baserateTaskSem);
    multiModeQuad_ROS_step();

    // Get model outputs here
    stopRequested = !((rtmGetErrorStatus(multiModeQuad_ROS_M) == (NULL)) &&
                      !rtmGetStopRequested(multiModeQuad_ROS_M));
    runModel = !stopRequested;
  }

  runModel = 0;
  terminateTask(arg);
  pthread_exit((void *)0);
  return NULL;
}

void exitFcn(int sig)
{
  UNUSED(sig);
  rtmSetErrorStatus(multiModeQuad_ROS_M, "stopping the model");
}

void *terminateTask(void *arg)
{
  UNUSED(arg);
  terminatingmodel = 1;

  {
    runModel = 0;
  }

  // Terminate model
  multiModeQuad_ROS_terminate();
  sem_post(&stopSem);
  return NULL;
}

int main(int argc, char **argv)
{
  UNUSED(argc);
  UNUSED(argv);
  void slros_node_init(int argc, char** argv);
  slros_node_init(argc, argv);
  rtmSetErrorStatus(multiModeQuad_ROS_M, 0);

  // Initialize model
  multiModeQuad_ROS_initialize();

  // Call RTOS Initialization function
  myRTOSInit(0.005, 0);

  // Wait for stop semaphore
  sem_wait(&stopSem);

#if (MW_NUMBER_TIMER_DRIVEN_TASKS > 0)

  {
    int i;
    for (i=0; i < MW_NUMBER_TIMER_DRIVEN_TASKS; i++) {
      CHECK_STATUS(sem_destroy(&timerTaskSem[i]), 0, "sem_destroy");
    }
  }

#endif

  return 0;
}

//
// File trailer for generated code.
//
// [EOF]
//
