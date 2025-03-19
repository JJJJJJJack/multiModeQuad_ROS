#ifndef _SLROS_INITIALIZE_H_
#define _SLROS_INITIALIZE_H_

#include "slros_busmsg_conversion.h"
#include "slros_generic.h"
#include "multiModeQuad_ROS_types.h"

extern ros::NodeHandle * SLROSNodePtr;
extern const std::string SLROSNodeName;

// For Block multiModeQuad_ROS/Flight mode
extern SimulinkSubscriber<std_msgs::Int16, SL_Bus_multiModeQuad_ROS_std_msgs_Int16> Sub_multiModeQuad_ROS_472;

// For Block multiModeQuad_ROS/Sub setpoint attitude
extern SimulinkSubscriber<geometry_msgs::PoseStamped, SL_Bus_multiModeQuad_ROS_geometry_msgs_PoseStamped> Sub_multiModeQuad_ROS_496;

// For Block multiModeQuad_ROS/Sub setpoint rate
extern SimulinkSubscriber<geometry_msgs::TwistStamped, SL_Bus_multiModeQuad_ROS_geometry_msgs_TwistStamped> Sub_multiModeQuad_ROS_497;

// For Block multiModeQuad_ROS/Sub setpoint thrust
extern SimulinkSubscriber<std_msgs::Float32, SL_Bus_multiModeQuad_ROS_std_msgs_Float32> Sub_multiModeQuad_ROS_500;

// For Block multiModeQuad_ROS/Sub setpoint velocity
extern SimulinkSubscriber<geometry_msgs::Twist, SL_Bus_multiModeQuad_ROS_geometry_msgs_Twist> Sub_multiModeQuad_ROS_426;

// For Block multiModeQuad_ROS/Publish
extern SimulinkPublisher<geometry_msgs::PoseStamped, SL_Bus_multiModeQuad_ROS_geometry_msgs_PoseStamped> Pub_multiModeQuad_ROS_438;

// For Block multiModeQuad_ROS/Publish1
extern SimulinkPublisher<sensor_msgs::Imu, SL_Bus_multiModeQuad_ROS_sensor_msgs_Imu> Pub_multiModeQuad_ROS_477;

// For Block multiModeQuad_ROS/Get Parameter
extern SimulinkParameterGetter<real64_T, double> ParamGet_multiModeQuad_ROS_459;

// For Block multiModeQuad_ROS/Get Parameter1
extern SimulinkParameterGetter<real64_T, double> ParamGet_multiModeQuad_ROS_460;

// For Block multiModeQuad_ROS/Get Parameter2
extern SimulinkParameterGetter<real64_T, double> ParamGet_multiModeQuad_ROS_461;

void slros_node_init(int argc, char** argv);

#endif
