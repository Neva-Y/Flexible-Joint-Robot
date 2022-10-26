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
A_aug = [zeros(1, 1) [0 1 0 0 0] ;
         zeros(5, 1) A ];

B_aug = [zeros(1,1); B];


controller_poles = [-15+2.5j, -15-2.5j, -60, -70, -80, -90];
K = place(A_aug, B_aug, controller_poles);

Ki = K(1);
K = K(2:6);


% Discrete controller design
T = 1/50;
Gp = ss(A, B, eye(5), 0);
Dp = c2d(Gp, T, 'zoh');
Phi = Dp.a; Gam = Dp.b;
[n, n] = size(Phi);
Phia = [1 [0 1 0 0 0]; zeros(n,1) Phi];
Gama = [0; Gam];
dclz = exp(controller_poles*T);
K_disc = place(Phia, Gama, dclz);

Ki_disc = K_disc(1);
K_disc = K_disc(2:6);

% Nxu = [Phi - eye(n) Gam; [0 1 0 0 0] 0]\[zeros(n,1); 1]
% Nx = Nxu(1:n)
% Nu = Nxu(n+1)


% Observer Design

C = [1 0 0 0 0;
     0 1 0 0 0];

% Transpose for duality when doing observer design
observer_poles = [-9600 -10800 -12000 -13200 -14400];
L = place(A', C', observer_poles)';

est_poles = eig(A - L*C)

% Reduced order observer design (ENCODER & TACHOMETER)
% A_new = [A(1,1) A(1,3) A(1,2) A(1,4) A(1,5);
%          A(3,1) A(3,3) A(3,2) A(3,4) A(3,5);
%          A(2,1) A(2,3) A(2,2) A(2,4) A(2,5);
%          A(4,1) A(4,3) A(4,2) A(4,4) A(4,5);
%          A(5,1) A(5,3) A(5,2) A(5,4) A(5,5)]

%Reduced order observer design (ENCODER & VISION) Continuous
% A_new = A
% 
% A_measured = [A_new(1,:); A_new(2,:)]
% A_estimated = [A_new(3,:); A_new(4,:); A_new(5,:)]
% 
% A_11 = A_measured(:,1:2)
% A_12 = A_measured(:,3:5)
% A_21 = A_estimated(:,1:2)
% A_22 = A_estimated(:,3:5)
% 
% B_1 = [B(1); B(2)]
% B_2 = [B(3); B(4); B(5)]
% 
% J_r = place(A_22', A_12', )'


%Reduced order observer design (ENCODER & VISION) Discrete
A_new = Dp.a
A_measured = [A_new(1,:); A_new(2,:)]
A_estimated = [A_new(3,:); A_new(4,:); A_new(5,:)]

A_11 = A_measured(:,1:2)
A_12 = A_measured(:,3:5)
A_21 = A_estimated(:,1:2)
A_22 = A_estimated(:,3:5)
B_1 = [Dp.b(1); Dp.b(2)]
B_2 = [Dp.b(3); Dp.b(4); Dp.b(5)]
obsz = exp([-700 -900 -1200]*T);
J_r = place(A_22', A_12', obsz)'

% dclz = exp([-5.5-0.5i, -5.5+0.5i -90, -100, -110, -120]*T);
% K_disc = place(Phia, Gama, dclz);
% Ki_disc = K_disc(1);
% K_disc = K_disc(2:6);



% Sinusoidal disturbance model
w = 2*pi*1/4;
A_d = [0 w; -w 0]
C_d = [1 0]


% A_dist = [A B*C_d; zeros(2, 5) A_d];
% B_dist = [B; zeros(2, 1)];
% C_dist = [C zeros(2, 2)];
% D_dist = 0;

%K_sin = place(A, B, controller_poles(1:5));
%J_sin = place(A_dist', C_dist', [observer_poles -15400 -16400]/5)';
%G_tf = tf(ss(A -B*K_sin, B, [0 1 0 0 0], 0));



A_dist = [A zeros(5,2);
           C_d'*[0 1 0 0 0] A_d']
B_dist = [B; zeros(2,1)]
C_dist = [0 1 0 0 0 zeros(1,2)]
K_dist = place(A_dist, B_dist, [controller_poles -100]/2)
K_sin = K_dist(1:5)
Kd_sin = K_dist(6:7)


