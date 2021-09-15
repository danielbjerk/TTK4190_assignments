%% Oppg. 2b)
sys = ss(A,B,C,D);
sys_tf = tf(sys);
delta_c_to_p = sys_tf(3);
delta_c_to_delta = tf([7.5], [1 7.5]);
delta_to_p = zpk(delta_c_to_p/delta_c_to_delta);


%% Oppg. 2c)
k_i_min = -0.2;
k_i_max = 2;
res = 0.2;
i_max = (k_i_max - k_i_min)/res;

i_zeros = zeros(3, i_max);
i_poles = zeros(3, i_max);

for i = 1:i_max
    k_i_phi = k_i_min + i*res;
    phi_c_to_phi = a_phi2*tf([k_p_phi k_i_phi],[1 (a_phi1+k_d_phi*a_phi2) k_p_phi*a_phi2 k_i_phi]);
    
    i_zeros(:,i) = zero(phi_c_to_phi);
    i_poles(:,i) = pole(phi_c_to_phi);
end


plot(reshape(i_poles, [1, numel(i_poles)]), 'o')
xlabel("Re")
ylabel("Im")
xlim([-2 0.5])
title("Root locus of k_{i,\phi} from -0.2 to 2")
grid on
%print("oppg2c_boundary", "-dpng")

%% Oppg. 2e)
chi_c = deg2rad(30);