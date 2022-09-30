% Linear plant model
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
A = [0 0 1 0 0;
    0 0 0 1 0;
    -Ks/Jm Ks/Jm 0 0 Km/Jm;
    Ks/Jl -Ks/Jl-m*g*l*cos(x2e)/Jl 0 0 0;
    0 0 -Km/La 0 -Ra/La];


B = [0;0;0;0;1/La];

D = 0;

tolerance = 1e-3;

% Controllability
controllability = rank(ctrb(A,B), tolerance)

% Observability
% Absolute Encoder
C_enc = zeros(5);
C_enc(1,1) = 1;
Obsv_enc = rank(obsv(A,C_enc), tolerance)

% Vision Sensor
C_vis = zeros(5);
C_vis(2,2) = 1;
Obsv_vis = rank(obsv(A,C_vis), tolerance)

% Poentiometer
C_pot = C_enc;
Obsv_pot = rank(obsv(A,C_enc), tolerance)


% Tachometer
C_tach = zeros(5);
C_tach(3,3) = 1;
Obsv_tach = rank(obsv(A,C_tach), tolerance)

% Gyroscope
C_gyro = zeros(5);
C_gyro(4,4) = 1;
Obsv_gyro = rank(obsv(A,C_gyro), tolerance)

% Vision and Encoder
rank(obsv(A,C_enc+C_vis), tolerance)

% Gyro and Encoder
rank(obsv(A,C_enc+C_gyro), tolerance)

% Tachometer and Vision
rank(obsv(A,C_tach+C_vis), tolerance)

% Vision and Gyro
rank(obsv(A,C_gyro+C_vis), tolerance)

% Gyro and Potentiometer
rank(obsv(A,C_pot+C_vis), tolerance)

% Controller Design (Need Integrator since IMP)
C = [0 1 0 0 0];

A_aug = [zeros(1, 1) [0 1 0 0 0] ;
         zeros(5, 1) A ];

B_aug = [zeros(1,1); B];


controller_poles = [-10+1.5j, -10-1.5j, -40, -45, -50, -55];
K = place(A_aug, B_aug, controller_poles);

Ki = K(1);
K = K(2:6);

% Observer Design

C = [1 0 0 0 0;
     0 0 1 0 0];

% Transpose for duality when doing observer design
observer_poles = [-800 -900 -1000 -1100 -1200]*3
L = place(A', C', observer_poles)'

est_poles = eig(A - L*C)