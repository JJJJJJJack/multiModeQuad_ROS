# multiModeQuad_ROS
A generic quadcopter model based on MATLAB/Simulink model that accept various control modes. The repo is a ros package that auto generated by Simulink codegen. To modify the model, change the simulink file instead of the source c++ code.

## Control command
Send `~/setpoint_velocity/cmd_vel_unstamped` (geometry_msgs/Twist) to the node to control the velocity of the vehicle.
Send `~/setpoint_attitude/attitude` (geometry_msgs/PoseStamped) to the node to control the attitude (quaternion) of the vehicle.
Send `~/setpoint_attitude/cmd_vel` (geometry_msgs/TwistStamped) to the node to control the angular rate of the vehicle.
Send `~/thrust` (std_msgs/Float32) to the node to control the angular rate of the vehicle.

