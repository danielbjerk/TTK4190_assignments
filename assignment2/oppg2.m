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

step_5_degrees_data = load('dataset/task_2e/step_5_degrees.mat').ans;
step_10_degrees_data = load('dataset/task_2e/step_10_degrees.mat').ans;
step_15_degrees_data = load('dataset/task_2e/step_15_degrees.mat').ans;

tiledlayout(3, 1);

nexttile;
plot_task_2efg(5, step_5_degrees_data, 50, -30, 10);

nexttile;
plot_task_2efg(10, step_10_degrees_data, 50, -30, 15);

nexttile;
plot_task_2efg(15, step_15_degrees_data, 50, -30, 25);


%% Oppg. 2f)
step_5_degrees_data = load('dataset/task_2f/step_5_degrees.mat').ans;
step_10_degrees_data = load('dataset/task_2f/step_10_degrees.mat').ans;
step_15_degrees_data = load('dataset/task_2f/step_15_degrees.mat').ans;

tiledlayout(3, 1);

nexttile;
plot_task_2efg(5, step_5_degrees_data, 140, -30, 10);

nexttile;
plot_task_2efg(10, step_10_degrees_data, 140, -30, 15);

nexttile;
plot_task_2efg(15, step_15_degrees_data, 140, -30, 25);


%% Oppg. 2g)
with_anti_integrator_windup = load('dataset/task_2g/step_30_with_anti_integrator_windup.mat').ans;
without_anti_integrator_windup = load('dataset/task_2g/step_30_without_anti_integrator_windup.mat').ans;

tiledlayout(2, 1);

nexttile;
plot_task_2efg(30, without_anti_integrator_windup, 105, -30, 55);

nexttile;
plot_task_2efg(30, with_anti_integrator_windup, 88, -30, 30);

%%
function plot_task_2efg(reference_course, data, data_end, min_y, max_y)
    time = data(1, 1:data_end);
    course = data(3, 1:data_end) * 180 / pi;
    delta_aileron = data(2, 1:data_end) * 180 / pi;

    plot(time, course);
    hold on;
    plot(time, delta_aileron);
    set(gca,'ytick',linspace(min_y,max_y, (max_y - min_y) / 5 + 1))
    ylim([min_y - 5, max_y + 5]);
    grid on;
    xlabel("time");
    legend('$\chi$', '$\delta_a$','Interpreter','latex')
    title("Course angle setpoint " + reference_course + " degrees")
end


