
format long g 
%%%% State space model
%%
% Linearised Plant at 45 degrees
% Using motor 1
% MOTOR PARAMETERS
Jm = 0.0021;
Km = 0.5369;
Ra = 2.6;
La = 0.18*10^-3;

% % Using motor 2
% % MOTOR PARAMETERS
% Jm = 7.5e-6;
% Km = 27e-3;
% Ra = 2;
% La = 1.3e-3;

% LINK PARAMETERS
Jl = 0.001;
m = 0.1;
l = 0.1;
g = 9.8;
Ks = 1.55;

e = deg2rad(45);

% STATE SPACE DEFINITIONS
A = [0 0 1 0 0;
     0 0 0 1 0;
     -Ks/Jm Ks/Jm 0 0 Km/Jm;
     Ks/Jl -Ks/Jl-m*g*l*cos(e)/Jl 0 0 0;
     0 0 -Km/La 0 -Ra/La];

B = [0;0;0;0;1/La];

C = [0 1 0 0 0; 0 0 1 0 0];

sizeC = size(C);
sizeB = size(B);
D = zeros(sizeC(1),sizeB(2));


s = tf('s');
% obsv(A,C)
obsv_rank = rank(obsv(A,C),1e-15);

H = C*inv(s*eye(5) - A)*B;

desired_poles = [-8, -18, -20, -60, -70];

Ts = 1/50;

desired_z_poles = exp(desired_poles.*Ts);

K = place(A,B,desired_poles)

K_z = place(A,B,desired_z_poles);


A_t = A - B*K;

C_ff = eye(5);
Gcl_eq = C_ff*inv(- A_t)*B;

N = 1./Gcl_eq;


desired_obsv_poles = [-600, -700, -800, -900, -1000];


L = place(A',C',desired_obsv_poles)'


integrator_desired_poles = [-6, -8, -30, -40, -50, -60];

% A_aug = [A zeros(5, 1);
%          [0 1 0 0 0] zeros(1, 1)];
% 
% B_aug = [B; zeros(1,1)];
% 
% Kcont= place(A_aug, B_aug, integrator_desired_poles);
% K_cont = Kcont(1:5)
% Ki_cont = Kcont(6)
% 



A_aug = [zeros(1, 1) [0 1 0 0 0] ;
         zeros(5, 1) A ];

B_aug = [zeros(1,1); B];

Kcont= place(A_aug, B_aug, integrator_desired_poles);

K_cont = Kcont(2:6)
Ki_cont = Kcont(1)





