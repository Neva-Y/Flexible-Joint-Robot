%% DEFINE THE PLANT INITIAL CONDITIONS

% Initial link angle [radians]
theta_l_init = deg2rad(0);
% Initial link velocity [radians/second]
thetadot_l_init = deg2rad(0);
% Initial motor angle [radians]
theta_m_init = deg2rad(0);
% Initial motor velocity [radians/second]
thetadot_m_init = deg2rad(0);
% Initial motor current [Amps]
Ia_init = 0;

% Initial conditions stacked as a 5x1 vector
init = [theta_m_init; theta_l_init; thetadot_m_init; thetadot_l_init; Ia_init];

%% DEFINE THE SENSOR PARAMETERS

% FOR VISION:
% > The frames per second [fps]
% > Allowed values: 25, 50, 100
vision_fps = 25;

% FOR ABSOLUTE ENCODER
% > Number of counts per revolution
absolute_encoder_counts_per_revolution = 128;

% FOR POTENTIOMETER
% > Percentage linearity [percent]
potentiometer_linearity_percent = 20;

% FOR TACHOMETER
% > Number of counts per revolution
tachometer_counts_per_revolution = 2048;
% > Frequency of the velocity estimates
% > Allowed values: 10, 25, 50, 100
tachometer_output_frequency = 50;

% FOR GYROSCOPE
% > Variance of the noise that drives the drift process [rad/s^2]
gyro_sigma_drift = deg2rad(5);