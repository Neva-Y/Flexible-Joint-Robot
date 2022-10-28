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

A_aug = [zeros(1, 1) [0 1 0 0 0] ;
         zeros(5, 1) A ];

B_aug = [zeros(1,1); B];

% Discrete Augmented matrix 
T = 1/50;
Gp = ss(A, B, eye(5), 0);
Dp = c2d(Gp, T, 'zoh');
Phi = Dp.a; Gam = Dp.b;
[n, m] = size(Phi);
Phia = [1 [0 1 0 0 0]; zeros(n,1) Phi];
Gama = [0; Gam];

%% Discrete observer design
C = [1 0 0 0 0;
     0 0 0 1 0];
W_disc = eye(5)*10000; % Model Noise
V_disc = diag([.0001, 0.1]);  % Sensor Error
[P, L_disct, O] = idare(Dp.a',C',W_disc,V_disc);
L_disc = L_disct';

%% STATE AND INPUT SIZES
nx = size(A_aug,2);
nu = size(B_aug,2);

%% SPECIFY THE Q AND R WEIGHTING MATRICES
% > for penalising each state as a diagonal matrix
Q = diag([10000, 1000, 1000, 30000, 30000, 1]);
% > for penalising the input
R = 1;

%% SPECIFY THE SAMPLE TIME OF THE CONTROLLER
sample_time = 1/50;

%% COMPUTE THE RICATTI SOLUTION
[P_lqr, K_lqr] = idare(Phia, Gama, Q, R);

%% SPECIFY THE TERMINAL COST WEIGHTING MATRIX
P = P_lqr;

%% SPECIFY THE C FOR DEFINING THE y_k PART TO BE USED IN THE CONSTRAINTS
C_for_constraints = [ 0,  0 , 1 , 0 , 0, 0 ];
ny = size(C_for_constraints,1);

%% CONSTRUCT THE WEIGHTING MATRICES FOR MATLAB MPC FORMULATION
% > for the state terms in the cost function
Wy = blkdiag( zeros(ny,ny) , eye(nx,nx) , zeros(nx) );
% > for the terminal state term in the cost function
WyN = blkdiag( zeros(ny,ny) , zeros(nx,nx) , eye(nx) );
% > for input terms in the cost function
Wu = chol(R);
% > for input change terms in the cost function
WDeltau = 0;


%% BUILD THE PLANT
% Specify the C and D
C_for_mpc = [ C_for_constraints ; chol(Q) ; chol(P) ];
D_for_mpc = zeros(ny+nx+nx,nu);
% Build the continuous-time plant
LTI_plant_continuous_time = ss(A_aug, B_aug, C_for_mpc, D_for_mpc);
% Convert to a discrete-time plant
LTI_plant_discrete_time = c2d(LTI_plant_continuous_time, sample_time, 'zoh' );

%% PUT THE WEIGHTINGS INTO A STRUCT
W = struct();
W.ManipulatedVariables     = diag(Wu).';
W.ManipulatedVariablesRate = diag(WDeltau).';
W.OutputVariables          = diag(Wy).';

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
terminal_setting_for_OV.Weight = diag(WyN).';
% > For the manipulated variables
terminal_setting_for_MV = struct();
terminal_setting_for_MV.Weight = W.ManipulatedVariables;
% > Call the function to set
setterminal(mpc_obj, terminal_setting_for_OV, terminal_setting_for_MV);

%% SET THE STATE ESTIMATOR
setEstimator(mpc_obj,'custom');

%% SET THE OUTPUT DISTURBANCE TO ZERO
setoutdist(mpc_obj,'model',ss(zeros(ny+nx+nx,1)));

%% CREATE THE MPC INITIAL STATE VARIABLE
% > Specify the initial condition for each state
theta_l_init     = deg2rad(0);
theta_m_init     = deg2rad(0);
theta_l_dot_init = deg2rad(0);
theta_m_dot_init = deg2rad(0);
Integral_error_init = 0;
Ia_init = 0;

% > Create the variable for the MPC block
mpc_initial_state = mpcstate(mpc_obj);
mpc_initial_state.Plant = [Integral_error_init, theta_m_init, theta_l_init, theta_m_dot_init, theta_l_dot_init,  Ia_init];

%% DISPLAY THE MPC OBJECT
display(mpc_obj)

% %% COMPUTE FEED-FORWARD CONVERSION VIA PLANT INVERSION
% % Specify the link angle as the reference
% C_for_theta_l = [0, 1, 0, 0];
% % Perform the inversion
% N_ff = inv( [A, B; C_for_theta_l, 0] );
% % Set any inversion artifact to zero
% N_ff( abs(N_ff)<1e-8 ) = 0;
% % Extract the state and input parts
% Nx_ff = N_ff(1:4,5);
% Nu_ff = N_ff(5,5);
