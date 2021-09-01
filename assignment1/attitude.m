% M-script for numerical integration of the attitude dynamics of a rigid 
% body represented by unit quaternions. The MSS m-files must be on your
% Matlab path in order to run the script.
%
% System:                      .
%                              q = T(q)w
%                              .
%                            I w - S(Iw)w = tau
% Control law:
%                            tau = constant
% 
% Definitions:             
%                            I = inertia matrix (3x3)
%                            S(w) = skew-symmetric matrix (3x3)
%                            T(q) = transformation matrix (4x3)
%                            tau = control input (3x1)
%                            w = angular velocity vector (3x1)
%                            q = unit quaternion vector (4x1)
%
% Author:                   2018-08-15 Thor I. Fossen and H�kon H. Helgesen

%% USER INPUTS
h = 0.1;                     % sample time (s)
t_end = 1000;
N  = t_end/h;                    % number of samples. Should be adjusted

% model parameters
m = 180;
r = 2;
I = m*r^2*eye(3);            % inertia matrix
I_inv = inv(I);

% constants
deg2rad = pi/180;   
rad2deg = 180/pi;

phi = -5*deg2rad;            % initial Euler angles
theta = 10*deg2rad;
psi = -20*deg2rad;

q = euler2q(phi,theta,psi);   % transform initial Euler angles to q

w = [0 0 0]';                 % initial angular rates

%table = zeros(N+1,14);        % memory allocation
table = zeros(N+1,20);        % memory allocation

print_results = true;


%% Problem 1.1 (Linearization)
eps_0 = [1;0;0;0];
omega_0 = [0;0;0];
x_0 = [eps_0; omega_0];

%[A, B] = linearize_satellite_dynamics(I, x_0);
A = [ 0 0 0 1 0 0;
      0 0 0 0 1 0;
      0 0 0 0 0 1;
      0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0];

B = [ 0 0 0;
      0 0 0;
      0 0 0;
      I_inv];


%% Problem 1.2 (Stability)
k_p = 2;
k_d = 40;
K_d = k_d*eye(3);
K = [k_p*eye(3), K_d];

%eig(A+B*K)


%% Problem 1.3 (Control to zero)
K_r = [ 1 0 0 0 0 0 ;
        0 1 0 0 0 0;
        0 0 1 0 0 0
        ];
q_d = euler2q(pi, 0, 0);
r = [q_d(2:end);0;0;0];


%% Problem 1.5 (Control with reference)
k_p = 20;
k_d = 400;
K_d = k_d*eye(3);


%% FOR-END LOOP
for i = 1:N+1
   t = (i-1)*h;                     % time

   [phi,theta,psi] = q2euler(q); % transform q to Euler angles
   [J,J1,J2] = quatern(q);       % kinematic transformation matrices
   
   % Calculation of attitude error
   eul_d = [0; 15*cos(0.1*t); 10*sin(0.05*t)]*deg2rad;
   q_d = euler2q(eul_d(1), eul_d(2), eul_d(3));
   q_d(2:end) = -q_d(2:end);
   q_err = quatprod(q_d,q);
   eps_err = q_err(2:end);
   
   % Calculation of angular velocity error
   % (Manually differentiated)
   eul_rate_d = [0; (-1)*(0.1)*15*sin(0.1*t); (0.05)*10*cos(0.05*t)]*deg2rad;
   [ex1, ex2, T_theta] = eulerang(eul_d(1), eul_d(2), eul_d(3));
   w_d = T_theta\eul_rate_d;
   w_err = w - w_d;
   
   % Control law
   %tau = [0.5 1 -1]';              
   %tau = -K*[q(2:end,1); w];       % Task 1.3
   %tau = K_r*r-K*[q(2:end,1); w];  % Task 1.3 with test reference
   %tau = -K_d*w - k_p*eps_err;      % Task 1.5
   tau = -K_d*w_err - k_p*eps_err;      % Task 1.6
   
   q_dot = J2*w;                        % quaternion kinematics
   w_dot = I_inv*(Smtrx(I*w)*w + tau);  % rigid-body kinetics
   
   table(i,:) = [t q' phi theta psi w' tau' eul_d' w_d'];  % store data in table
   
   q = q + h*q_dot;	             % Euler integration
   w = w + h*w_dot;
   
   q  = q/norm(q);               % unit quaternion normalization
end 

%% PLOT FIGURES
close all

t       = table(:,1);  
q       = table(:,2:5); 
phi     = rad2deg*table(:,6);
theta   = rad2deg*table(:,7);
psi     = rad2deg*table(:,8);
w       = rad2deg*table(:,9:11);  
tau     = table(:,12:14);
phi_d = rad2deg*table(:,15);
theta_d = rad2deg*table(:,16);
psi_d = rad2deg*table(:,17);
w_d = rad2deg*table(:,18:20);

problem = "problem_1,6_";


figure (1); clf;
hold on;
plot(t, phi, 'b');
plot(t, theta, 'r');
plot(t, psi, 'g');
hold off;
grid on;
legend('\phi', '\theta', '\psi');
title('Euler angles');
xlabel('time [s]'); 
ylabel('angle [deg]');
if print_results
    print(strcat(problem, "euler_angles"), "-dpng")
end


figure (2); clf;
hold on;
plot(t, w(:,1), 'b');
plot(t, w(:,2), 'r');
plot(t, w(:,3), 'g');
hold off;
grid on;
legend('x', 'y', 'z');
title('Angular velocities');
xlabel('time [s]'); 
ylabel('angular rate [deg/s]');
if print_results
    print(strcat(problem, "angular_rate"), "-dpng")
end


figure (3); clf;
hold on;
plot(t, tau(:,1), 'b');
plot(t, tau(:,2), 'r');
plot(t, tau(:,3), 'g');
hold off;
grid on;
legend('x', 'y', 'z');
title('Control input');
xlabel('time [s]'); 
ylabel('input [Nm]');
if print_results
    print(strcat(problem, "control_input"), "-dpng")
end


figure (4); clf;
hold on;
plot(t, phi-phi_d, 'b');
plot(t, theta-theta_d, 'r');
plot(t, psi-psi_d, 'g');
hold off;
grid on;
legend('\phi _{err}', '\theta _{err}', '\psi _{err}');
title('Heading tracking error');
xlabel('time [s]'); 
ylabel('angle [deg]');
if print_results
    print(strcat(problem, "heading_tracking_error"), "-dpng")
end


figure (5); clf;
hold on;
plot(t, w(:,1)-w_d(:,1), 'b');
plot(t, w(:,2)-w_d(:,2), 'r');
plot(t, w(:,3)-w_d(:,3), 'g');
hold off;
grid on;
legend('p _{err}', 'q_{err}', 'r _{err}');
title('Velocity tracking error');
xlabel('time [s]'); 
ylabel('angle [deg]');
if print_results
    print(strcat(problem, "velocity_tracking_error"), "-dpng")
end


figure (6); clf;
subplot(3,1,1)
hold on
title('Controller performance');
plot(t, w(:,1), 'b');
plot(t, w_d(:,1), 'k--');
xlabel('time [s]'); 
ylabel('p [deg]');
hold off

subplot(3,1,2)
hold on
plot(t, w(:,2), 'r');
plot(t, w_d(:,2), 'k--');
xlabel('time [s]'); 
ylabel('q [deg]');
hold off

subplot(3,1,3)
hold on
plot(t, w(:,3), 'g');
plot(t, w_d(:,3), 'k--');
xlabel('time [s]'); 
ylabel('r [deg]');
hold off
if print_results
    print(strcat(problem, "velocity_controller_performance"), "-dpng")
end
