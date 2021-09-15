%% Enviromental constants
g = 9.81;
V_a = 580/3.6; % [m/s]
V_g = 580/3.6; % [km/h] ???????????????
d = deg2rad(1.5);


%% System constants
delta_a_max = deg2rad(30);


%% Linear system
A = [   -0.322,     0.052,      0.028,      -1.12,      0.002;
        0,          0,          1,          -0.001,     0;
        -10.6,      0,          -2.87,      0.46,       -0.65;
        6.87,       0,          -0.04,      -0.32,      -0.02;
        0,          0,          0,          0,          -7.5];
    
B = [   0;
        0;
        0;
        0;
        7.5];
    
C = [   1,  0,  0,  0,  0;
        0,  1,  0,  0,  0;
        0,  0,  1,  0,  0;
        0,  0,  0,  1,  0];
    
D = [   0;
        0;
        0;
        0];
    

%% Control constants
a_phi1 = 2.879;
a_phi2 = -0.65;

% Successive loop closure
delta_a_max = delta_a_max;
e_phi_max = deg2rad(15);
zeta_phi = 0.707;
omega_n_phi = sqrt(abs(a_phi2)*delta_a_max/e_phi_max);

zeta_chi = 0.707;
omega_n_chi = omega_n_phi/10;

%   phi
k_p_phi = delta_a_max/e_phi_max * sign(a_phi2);
k_d_phi = (2*zeta_phi*omega_n_phi - a_phi1)/(a_phi2);
k_i_phi = 0;
%   chi
k_p_chi = 2*zeta_chi*omega_n_chi*V_g/g;
k_i_chi = omega_n_chi^2*V_g/g;