% Project in TTK4190 Guidance, Navigation and Control of Vehicles 
%
% Author:           My name
% Study program:    My study program

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task = "2a";

h  = 0.1;    % sampling time [s]
Ns = 10000;    % no. of samples

psi_ref = 10 * pi/180;  % desired yaw angle (rad)
U_ref   = 7;            % desired surge speed (m/s)

L_oa = 161;             % lenght of ship (m)

% initial states
nu_b_0  = [0.1 0 0]';     % initial velocity
eta_n_0 = [0 0 0]';       % initial position & heading
delta = 0;  
n = 0;
x = [nu_b_0' eta_n_0' delta n]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simdata = zeros(Ns+1,16);       % table of simulation data

for i=1:Ns+1
    
    t = (i-1) * h;              % time (s)
    
    % current state
    nu_b = x(1:3);
    eta_n = x(4:6);
    
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
    
    % reference models
    psi_d = psi_ref;
    r_d = 0;
    u_d = U_ref;
        
    % control law
    delta_c = 0.1;              % rudder angle command (rad)
    n_c = 10;                   % propeller speed (rps)
    
    % ship dynamics
    u = [delta_c n_c]';
    [xdot,u] = ship(x,u,nu_c_n,tau_wind);
    
    % store simulation data in a table (for testing)
    simdata(i,:) = [t x(1:3)' x(4:6)' x(7) x(8) u(1) u(2) u_d psi_d r_d beta_c beta];     
 
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
beta_c_data   = (180/pi) * simdata(:,15);     % deg
beta_data = (180/pi) * simdata(:,16);     % deg

figure(1)
figure(gcf)
subplot(311)
plot(y,x,'linewidth',2); axis('equal')
title('North-East positions (m)'); xlabel('(m)'); ylabel('(m)'); 
subplot(312)
plot(t,psi,t,psi_d,'linewidth',2);
title('Actual and desired yaw angles (deg)'); xlabel('time (s)');
subplot(313)
plot(t,r,t,r_d,'linewidth',2);
title('Actual and desired yaw rates (deg/s)'); xlabel('time (s)');

figure(2)
figure(gcf)
subplot(311)
plot(t,u,t,u_d,'linewidth',2);
title('Actual and desired surge velocities (m/s)'); xlabel('time (s)');
subplot(312)
plot(t,n,t,n_c,'linewidth',2);
title('Actual and commanded propeller speed (rpm)'); xlabel('time (s)');
subplot(313)
plot(t,delta,t,delta_c,'linewidth',2);
title('Actual and commanded rudder angles (deg)'); xlabel('time (s)');

figure(3)
figure(gcf)
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