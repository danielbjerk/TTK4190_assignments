function [A_lin, B_lin] = linearize_satellite_dynamics(I_g, x_0)
syms eta eps1 eps2 eps3 omega1 omega2 omega3 tau1 tau2 tau3
q = [eta; eps1; eps2; eps3];
omega = [omega1; omega2; omega3];
x = [eta; eps1; eps2; eps3; omega1; omega2; omega3];
tau = [tau1; tau2; tau3];

T_q = Tquat(q);
S = Smtrx(I_g * omega);

xdot = [T_q*omega; I_g\S*omega + I_g\tau];
    
A_lin = subs(diff(xdot,x), x, x_0);
B_lin = subs(diff(xdot,'tau'), x, x_0);
end