function assignment_main()
    experiment01();
end

function experiment01()
    % Simulate Mass-String System
    
    % simulation parameters
    amplitude_Uf = 0.025; % m
    omega_Uf = pi;
    
    num_masses = 10;
    total_mass = 0.1; % kg
    
    tension_force = 10; % N
    string_length = 1; % m
    damping_coeff = 0.01; % %kg/s
    dx = string_length / num_masses; % m
    
    % generate the struct
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;
    
    
    % Modal Analysis
    
    [M, K] = construct_2nd_order_matrices(string_params);
    
    [U, lambda] = eig(K, M);
    
    omega_index = 3;
    
    w_n = sqrt(lambda(omega_index,omega_index));
    
    Uf_func = @(t_in) amplitude_Uf*cos(w_n*t_in);
    dUfdt_func = @(t_in) -w_n*amplitude_Uf*sin(w_n*t_in);
    
    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;
    
    % Run sim at natural frequency
    
    V0 = zeros(1, num_masses*2); % initial conditions
    tspan = [0 30*2*pi/w_n]; % integration period
    simulate_system(string_params, V0, tspan, false, 'discrete-wave', w_n);
end


function experiment02()
    % Traveling wave
    
    string_params.n = 200;
    V0 = zeros(1, num_masses*2); % initial conditions
    tspan = [0 30*2*pi/w_n]; % integration period
    simulate_system(string_params, V0, tspan, false, 'discrete-wave', w_n);
end