#include "slros_initialize.h"

ros::NodeHandle * SLROSNodePtr;
const std::string SLROSNodeName = "multiModeQuad_ROS";

// For Block multiModeQuad_ROS/Flight mode
SimulinkSubscriber<std_msgs::Int16, SL_Bus_multiModeQuad_ROS_std_msgs_Int16> Sub_multiModeQuad_ROS_472;

// For Block multiModeQuad_ROS/Sub setpoint attitude
SimulinkSubscriber<geometry_msgs::PoseStamped, SL_Bus_multiModeQuad_ROS_geometry_msgs_PoseStamped> Sub_multiModeQuad_ROS_496;

// For Block multiModeQuad_ROS/Sub setpoint rate
SimulinkSubscriber<geometry_msgs::TwistStamped, SL_Bus_multiModeQuad_ROS_geometry_msgs_TwistStamped> Sub_multiModeQuad_ROS_497;

// For Block multiModeQuad_ROS/Sub setpoint thrust
SimulinkSubscriber<std_msgs::Float32, SL_Bus_multiModeQuad_ROS_std_msgs_Float32> Sub_multiModeQuad_ROS_500;

// For Block multiModeQuad_ROS/Sub setpoint velocity
SimulinkSubscriber<geometry_msgs::Twist, SL_Bus_multiModeQuad_ROS_geometry_msgs_Twist> Sub_multiModeQuad_ROS_426;

// For Block multiModeQuad_ROS/Publish
SimulinkPublisher<geometry_msgs::PoseStamped, SL_Bus_multiModeQuad_ROS_geometry_msgs_PoseStamped> Pub_multiModeQuad_ROS_438;

// For Block multiModeQuad_ROS/Publish1
SimulinkPublisher<sensor_msgs::Imu, SL_Bus_multiModeQuad_ROS_sensor_msgs_Imu> Pub_multiModeQuad_ROS_477;

// For Block multiModeQuad_ROS/Get Parameter
SimulinkParameterGetter<real64_T, double> ParamGet_multiModeQuad_ROS_459;

// For Block multiModeQuad_ROS/Get Parameter1
SimulinkParameterGetter<real64_T, double> ParamGet_multiModeQuad_ROS_460;

// For Block multiModeQuad_ROS/Get Parameter2
SimulinkParameterGetter<real64_T, double> ParamGet_multiModeQuad_ROS_461;

void slros_node_init(int argc, char** argv)
{
  ros::init(argc, argv, SLROSNodeName);
  SLROSNodePtr = new ros::NodeHandle();
}

