% Project in TTK4190 Guidance, Navigation and Control of Vehicles 
%
% Author:           My name
% Study program:    My study program

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

task = "";

h  = 0.1;    % sampling time [s]
Ns = 80000;%10000;    % no. of samples

U_ref   = 9;            % desired surge speed (m/s)


% ship parameters
L_oa = 161;             % lenght of ship (m)
delta_max  = 40 * pi/180; % Rudder angle

% initial states
nu_b_0  = [0.1 0 0]';     % initial velocity
eta_n_0 = [0 0 deg2rad(-110)]';       % initial position & heading
delta = 0;  
n = 0;
x = [nu_b_0' eta_n_0' delta n]';

% Controllers
% Yaw PID-controller design parameters
T = 169.5493;
K = 0.0075;

wb = 0.06;
zeta = 1;

wn = wb/(sqrt(1-2*zeta^2 + sqrt(4*zeta^4 - 4*zeta^2 + 2)));
Kp = wn^2*T/K;
Ki = wn^3*T/(10*K);
Kd = (2*zeta*wn*T - 1)/K;
% Surge velocity PID
Kp_u = 1;
Ki_u = 0.3;

e_psi_int = 0;      % integration state
e_u_int = 0;

% reference models
% yaw
x_r_y = [0;0;0];

w_r_y = 0.03;
zeta_r_y = 1;
A_r_y =     [   0,      1,      0;
                    0,      0,      1;
                    -w_r_y^3,   -(2*zeta_r_y + 1)*w_r_y^2,  -(2*zeta_r_y + 1)*w_r_y];
B_r_y = [0; 0; w_r_y^3];

% guidance objective
waypoints = load("WP.mat").WP;
threshold_radius = 80;
lookahead_delta = 3*L_oa;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simdata = zeros(Ns+1,17);       % table of simulation data

for i=1:Ns+1
    
    t = (i-1) * h;              % time (s)
    
    % current state
    nu_b = x(1:3);
    eta_n = x(4:6);
    
    psi = eta_n(3);
    r = nu_b(3);
    
    % current disturbance
    Vc = 1;
    beta_Vc = deg2rad(45);
    
    uc = Vc*cos(beta_Vc);
    vc = Vc*sin(beta_Vc);
    nu_c_n = [ uc vc 0 ]';        % Velocity of ocean currents NED {n}
    
    nu_c_b = (Rzyx(0, 0, eta_n(3))')*nu_c_n;   % Velocity of ocean currents BODY {b}
    
    % wind disturbance
    if t >= 200
        Vw = 10;
        beta_Vw = deg2rad(135);
        rho_a = 1.247;
        cy = 0.95;
        cn = 0.15;
        A_Lw = 10*L_oa;
        
        q = 0.5*rho_a*Vw^2;
        gamma_w = eta_n(3) - beta_Vw - pi;
        
        CY_gamma = cy*sin(gamma_w);
        CN_gamma = cn*sin(2*gamma_w);
        
        Ywind = q*CY_gamma*A_Lw;
        Nwind = q*CN_gamma*A_Lw*L_oa;
        
        %disp(CY_gamma)
        %disp(Ywind)
    else
        Ywind = 0;
        Nwind = 0;
    end
    tau_wind = [0 Ywind Nwind]';
    
    % heading
    if (nu_b(2)^2 + nu_b(1)^2) ~= 0
        nu_n = Rzyx(0, 0, eta_n(3))*nu_b;   % nu in {n}
        chi = atan2(nu_n(2), nu_n(1));
        beta_c = chi - eta_n(3);
        
        nu_r = nu_b - nu_c_b;
        beta = atan2(nu_r(2), nu_r(1));
    else
        chi = NaN;
    end
    
    % references
    % guidance law (yaw reference)
    chi_d = guidance(eta_n(1:2), waypoints, threshold_radius, lookahead_delta);
    psi_ref = chi_d;
    
    % reference models
    % yaw
    x_r_y = x_r_y + (A_r_y*x_r_y + B_r_y*psi_ref)*h;
    
    psi_d = x_r_y(1);
    r_d = x_r_y(2);
    % foreward velocity
    u_d = U_ref;
        
    % control laws
    % yaw
    e_psi = ssa(psi - psi_d);
    e_psi_int = e_psi_int + h*e_psi;
    
    delta_c = -(Kp*e_psi + Kd*r + Ki*e_psi_int); % rudder angle command (rad)  
    % integrator anti-windup (yaw)
    if abs(delta_c) >= delta_max
       delta_c = delta_max*sign(delta_c);
       e_psi_int = e_psi_int - h*e_psi; % Stop integration.
    end
    
    % surge
    e_u = nu_b(1) - u_d;
    e_u_int = e_u_int + h*e_u;
    n_c = -(Kp_u*e_u + Ki_u*e_u_int);                   % propeller speed (rps) 
    %(OBS! Summeres med feedforeward i ship.m)
    
    % ship dynamics
    u = [delta_c n_c]';
    [xdot,u] = ship(x,u,nu_c_n,tau_wind,u_d);
    
    % store simulation data in a table (for testing)
    simdata(i,:) = [t x(1:3)' x(4:6)' x(7) x(8) u(1) u(2) u_d psi_d r_d beta_c beta x_r_y(1)];      
 
    % Euler integration
    x = euler2(xdot,x,h);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t       = simdata(:,1);                 % s
u       = simdata(:,2);                 % m/s
v       = simdata(:,3);                 % m/s
r       = (180/pi) * simdata(:,4);      % deg/s
x       = simdata(:,5);                 % m
y       = simdata(:,6);                 % m
psi     = (180/pi) * simdata(:,7);      % deg
delta   = (180/pi) * simdata(:,8);      % deg
n       = 60 * simdata(:,9);            % rpm
delta_c = (180/pi) * simdata(:,10);     % deg
n_c     = 60 * simdata(:,11);           % rpm
u_d     = simdata(:,12);                % m/s
psi_d   = (180/pi) * simdata(:,13);     % deg
r_d     = (180/pi) * simdata(:,14);     % deg/s

figure(1)
figure(gcf)
subplot(311)
plot(y,x,'linewidth',2); axis('equal')
legend("Actual", "Desired");
title('North-East positions (m)'); xlabel('(m)'); ylabel('(m)'); 
subplot(312)
plot(t,psi,t,psi_d,'linewidth',2);
legend("Actual", "Desired");
title('Actual and desired yaw angles (deg)'); xlabel('time (s)');
subplot(313)
plot(t,r,t,r_d,'linewidth',2);
legend("Actual", "Desired");
title('Actual and desired yaw rates (deg/s)'); xlabel('time (s)');

%saveas(gcf, "Figures\oppg1c_integrator-windup.png")

figure(2)
figure(gcf)
subplot(311)
plot(t,u,t,u_d,'linewidth',2);
legend("Actual", "Desired");
title('Actual and desired surge velocities (m/s)'); xlabel('time (s)');
subplot(312)
plot(t,n,t,n_c,'linewidth',2);
legend("Actual", "Desired");
title('Actual and commanded propeller speed (rpm)'); xlabel('time (s)');
subplot(313)
plot(t,delta,t,delta_c,'linewidth',2);
legend("Actual", "Desired");
title('Actual and commanded rudder angles (deg)'); xlabel('time (s)');

pathplotter(x,y);

if task ~= ""
    figure(3)
    figure(gcf)
end
if task == "1b"
   hold on
   plot(t,beta_c_data,t,beta_data,'linewidth',2);
   hold off
   legend("\beta _c","\beta")
%    title('Crab angle vs. sideslip (deg) (V_c = 0)'); xlabel('time (s)');
   title('Crab angle vs. sideslip (deg) (V_c = 1)'); xlabel('time (s)');
   grid on
end

if task == "1c"
   hold on
   %plot(t,u,'linewidth',2);
   plot(t,u,t,u_d,'linewidth',2);
   hold off
   title('Actual and desired surge velocities (m/s)'); xlabel('time (s)');
   grid on
end

if task == "2d"
    hold on
   plot(t,x_r_y_data,'linewidth',2);
   hold off
   title('Reference model'); xlabel('time (s)');
   grid on
end
