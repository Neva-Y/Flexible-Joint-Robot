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

% Integral; theta_m; theta_l; theta_m_dot; theta_l_dot; current
% Controller Design (Need Integrator since IMP)
A_aug = [zeros(1, 1) [0 1 0 0 0] ;
         zeros(5, 1) A ];

B_aug = [zeros(1,1); B];

% SRL Design
C_aug = [0 1 0 0 0 0];
%syms s
%G = C_aug*inv(s*eye(6)-A_aug)*B_aug
[b,a] = ss2tf(A_aug, B_aug, C_aug, 0);
G = tf(b,a);
[num,den] = tfdata(G);
syms s
G_sym = poly2sym(cell2mat(num),s)/poly2sym(cell2mat(den),s);
G_sym_neg = subs(G_sym, s, -s);
G_SRL = syms2tf(G_sym*G_sym_neg);
rlocus(G_SRL)
rho = 700;

%Q_SRL = rho* (C' * C);
%R_SRL = 1;



Q = diag([100000, 1000, 1000, 100, 100, 1]);
R = 1;
Kr = lqr(A_aug, B_aug, Q, R);
Ki = Kr(1);
K = Kr(2:6);


% Observer Design
C = [1 0 0 0 0;
     0 0 0 1 0];

% Transpose for duality when doing observer design
W = eye(5)*1000; % Model Noise
V = diag([.01, 10]);  % Sensor Error
[P, L_t, O] = icare(A',C',W,V);
L = L_t';


% Discrete controller design
T = 1/50;
Gp = ss(A, B, eye(5), 0);
Dp = c2d(Gp, T, 'zoh');
Phi = Dp.a; Gam = Dp.b;
[n, m] = size(Phi);
Phia = [1 [0 1 0 0 0]; zeros(n,1) Phi];
Gama = [0; Gam];

% Solve discrete ARE
Q_disc = diag([10000, 1000, 1000, 30000, 30000, 1]);
R_disc = 1;
[X, Kr_disc] = idare(Phia, Gama, Q_disc, R_disc);
Ki_disc = Kr_disc(1);
K_disc = Kr_disc(2:6);

% Discrete observer design
W_disc = eye(5)*10000; % Model Noise
V_disc = diag([.0001, 0.1]);  % Sensor Error
[P, L_disct, O] = idare(Dp.a',C',W_disc,V_disc);
L_disc = L_disct';

