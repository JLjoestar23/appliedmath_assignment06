%% Simulate Mass-String System

% simulation parameters
amplitude_Uf = 0.05; % m
omega_Uf = w_n(5,5); % rad/s

num_masses = 10;
total_mass = 1; % kg
Uf_func = @(t_in) amplitude_Uf*cos(omega_Uf*t_in);
dUfdt_func = @(t_in) -omega_Uf*amplitude_Uf*sin(omega_Uf*t_in);
tension_force = 10; % N
string_length = 2; % m
damping_coeff = 0.05; % %kg/s
dx = string_length / num_masses + 1; % m

% generate the struct
string_params = struct();
string_params.n = num_masses;
string_params.M = total_mass;
string_params.Uf_func = Uf_func;
string_params.dUfdt_func = dUfdt_func;
string_params.Tf = tension_force;
string_params.L = string_length;
string_params.c = damping_coeff;
string_params.dx = dx;

%% Run sim
V0 = zeros(1, num_masses*2); % initial conditions
tspan = [0 30]; % integration period
simulate_system(string_params, V0, tspan, false, 'discrete-wave');

%% Modal Analysis

[M, K] = construct_2nd_order_matrices(string_params);

[U, lambda] = eig(K, M);

w_n = sqrt(lambda);