%% Linear plant model
% Motor 1 parameters
Jm = 0.0021;
Ra = 2.6;
La = 0.18*10^-3;
Km =  0.5369;

% Tune for system
Ks = 1.55; % [1.55, 1.80]
b = 0.013; %  [0.013, 0.019]
x2e = deg2rad(45);

Jl = 0.001;
m = 0.1;
g = 9.81;
l = 0.1;


% Linearise at 45 degrees
% A = [0 0 1 0 0;
%     0 0 0 1 0;
%     -Ks/Jm Ks/Jm 0 0 Km/Jm;
%     Ks/Jl -Ks/Jl-m*g*l*cos(x2e)/Jl 0 0 0;
%     0 0 -Km/La 0 -Ra/La];
% 
% B = [0;0;0;0;1/La];
% 

%% LTI SYSTEM MATRICES
A = [...
	 0,       0,              1,  0, 0;...
	 0,       0,              0,  1, 0;...
	-Ks/Jm,   Ks/Jm,          0,  0, Km/Jm;...
	 Ks/Jl,  -Ks/Jl-m*g*l*cos(x2e)/Jl,  0,  0, 0;...
     0    ,   0,         -Km/La,  0, -Ra/La];
B = [0; 0; 0; 0; 1/La];

%% Discrete observer design
T = 1/50;
Gp = ss(A, B, eye(5), 0);
Dp = c2d(Gp, T, 'zoh');
C = [1 0 0 0 0;
     0 0 0 1 0];
W_disc = eye(5)*10000; % Model Noise
V_disc = diag([.0001, 0.1]);  % Sensor Error
[P, L_disct, O] = idare(Dp.a',C',W_disc,V_disc);
L_disc = L_disct';

%% STATE AND INPUT SIZES
nx = size(A,2);
nu = size(B,2);

%% SPECIFY THE WEIGHTING MATRICES
% > for penalising each state
Wy = 10;
% > for penalising the input
Wu = 0;
% > for penalising input changes
WDeltau = 1;
% > for penalising the terminal state
WyN = 20;

%% BUILD THE PLANT
% Specify the C and D
C_for_mpc = [ 0 , 1 , 0 , 0, 0 ];
D_for_mpc = 0;
% Build the continuous-time plant
LTI_plant_continuous_time = ss(A,B,C_for_mpc,D_for_mpc);
% Convert to a discrete-time plant
sample_time = 1/50;
LTI_plant_discrete_time = c2d(LTI_plant_continuous_time, sample_time, 'zoh' );

%% OUTPUT VARIABLE SIZE
ny = size(C_for_mpc,1);

%% PUT THE WEIGHTINGS INTO A STRUCT
W = struct();
W.ManipulatedVariables     = Wu;
W.ManipulatedVariablesRate = WDeltau;
W.OutputVariables          = Wy;

%% SPECIFY THE MANIPULATED VARIABLE CONSTRAINTS
MV = struct();
MV.Min         = -6; %[V]
MV.Max         =  6; %[V]
MV.MinECR      = 0;
MV.MaxECR      = 0;
MV.RateMin     = -Inf;
MV.RateMax     = Inf;
MV.RateMinECR  = 0;
MV.RateMaxECR  = 0;
MV.Target      = 'nominal';
MV.Name        = 'Volts';
MV.Units       = 'V';
MV.ScaleFactor = 1;

%% SPECIFY THE OUTPUT VARIABLE CONSTRAINTS
OV = struct();
OV(1).Min         = deg2rad(25);
OV(1).Max         = deg2rad(65);
OV(1).MinECR      = 0.1;
OV(1).MaxECR      = 0.1;
OV(1).Name        = 'theta_l';
OV(1).Units       = 'radians';
OV(1).ScaleFactor = 1;

%% CREATE THE MPC OBJECT
prediction_horizon = 20;
control_horizon = prediction_horizon;
mpc_obj = mpc(LTI_plant_discrete_time, sample_time, prediction_horizon, control_horizon, W, MV, OV);

%% SET THE TERMINAL WEIGHT
% > For the output variables
terminal_setting_for_OV = struct();
terminal_setting_for_OV.Weight = WyN;
% > For the manipulated variables
terminal_setting_for_MV = struct();
terminal_setting_for_MV.Weight = W.ManipulatedVariables;
% > Call the function to set
setterminal(mpc_obj, terminal_setting_for_OV, terminal_setting_for_MV);

%% SET THE STATE ESTIMATOR
setEstimator(mpc_obj,'custom');

%% SET THE OUTPUT DISTURBANCE TO ZERO
setoutdist(mpc_obj,'model',ss(zeros(ny,1)));

%% CREATE THE MPC INITIAL STATE VARIABLE
% > Specify the initial condition for each state
theta_l_init     = deg2rad(0);
theta_m_init     = deg2rad(0);
theta_l_dot_init = deg2rad(0);
theta_m_dot_init = deg2rad(0);
Ia_init = 0;
% > Create the variable for the MPC block
mpc_initial_state = mpcstate(mpc_obj);
mpc_initial_state.Plant = [theta_m_init, theta_l_init, theta_m_dot_init, theta_l_dot_init,  Ia_init];

%% DISPLAY THE MPC OBJECT
display(mpc_obj)
