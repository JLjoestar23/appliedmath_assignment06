function assignment_main()
    experiment01();
    %experiment02();
    %experiment03();
    %experiment04();
end

function experiment01()
    % Simulate Mass-String System
    
    % simulation parameters
    amplitude_Uf = 0.025; % m
    omega_Uf = pi;
    
    num_masses = 20;
    total_mass = 0.1; % kg
    
    tension_force = 1; % N
    string_length = 1; % m
    damping_coeff = 0.001; % %kg/s
    %dx = string_length / num_masses; % m
    
    % generate the struct
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    %string_params.dx = dx;
    
    
    % Modal Analysis
    [M, K] = construct_2nd_order_matrices(string_params);
    
    [U, lambda] = eig(K, M);
    
    omega_index = 2;
    
    w_n = sqrt(lambda(omega_index,omega_index));
    
    Uf_func = @(t_in) amplitude_Uf*cos(w_n*t_in);
    dUfdt_func = @(t_in) -w_n*amplitude_Uf*sin(w_n*t_in);
    
    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;
    
    % Run sim at natural frequency
    
    V0 = zeros(1, num_masses*2); % initial conditions
    tspan = [0 40*2*pi/w_n]; % integration period
    simulate_system(string_params, V0, tspan, 2*w_n, false, false, 'discrete-wave');
end


function experiment02()
    % Traveling wave
    num_masses = 500;
    total_mass = 0.1; % kg
    
    tension_force = 10; % N
    string_length = 1; % m
    damping_coeff = 0; % %kg/s
    
    % generate the struct
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    
    w = 0.025;
    h = 0.25;
    Uf_func = @(t_in) b_spline_pulse(t_in, w, h);
    dUfdt_func = @(t_in) b_spline_pulse_derivative(t_in, w, h);
    
    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;
    
    rho = total_mass/string_length;
    c = sqrt(tension_force / rho);
    
    V0 = zeros(1, num_masses*2); % initial conditions
    tspan = [0 4*string_params.L/c]; % integration period
    simulate_system(string_params, V0, tspan, 40*c, false, false, 'colliding-wave');

end

function experiment03()
    % Compare analytical solution to numerical solution
    
    num_masses = 50;
    total_mass = 0.1; % kg
    tension_force = 10; % N
    string_length = 1; % m
    damping_coeff = 0; % %kg/s
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    
    % compute nth mode shape analytically
    x = linspace(0, string_length, 500);
    for n = 1:5
        mode_shape(n, :) = sin(pi*n*x/ string_length);
    end

    % compute nth mode shape numerically
    [M, K] = construct_2nd_order_matrices(string_params);
    
    [U, lambda] = eig(K, M);
    
    figure();
    for i=1:n
        % ensure numerical and analytical mode shapes have the same sign
        if U(end, i) * mode_shape(i, end) < 0
            U(:, i) = -U(:, i);
        end
        % 5 subplots for 5 mode shape comparisons
        subplot(5, 1, i);
        hold on;
        plot(x, mode_shape(i, :), '-', 'DisplayName', 'Analytical');
        plot(linspace(0, string_length, num_masses), U(:, i) / max(U(:, i)), '.', 'MarkerSize', 10, 'DisplayName', 'Numerical');
        legend();
        if i == 1
            title(sprintf("%dst Mode Shape Comparison", i));
        else
            title(sprintf("%dth Mode Shape Comparison", i));
        end
        hold off;
    end

    % frequencies from analytical solution
    rho = total_mass / string_length;
    c = sqrt(tension_force / rho);
    wn_idx = linspace(1, num_masses, 1000);
    wn_list = c*pi*wn_idx / string_length;

    figure();
    hold on;
    plot(linspace(1, num_masses, num_masses), sqrt(diag(lambda)), '.', 'MarkerSize', 20, 'DisplayName', 'Numerical');
    plot(wn_idx, wn_list, '-', 'LineWidth', 2, 'DisplayName', 'Analytical');
    xlabel('Mode Number');
    ylabel('Frequency (rad/s)');
    title('Frequency Comparison');
    legend();
    grid on;
    hold off;

end

function experiment04()
    
    figure();
    
    err = zeros(1,100);
    
    wn_idx = 1;   % resonant freq idx
    
    for n = 1:200
        num_masses = n;
        total_mass = 0.1; % kg
        tension_force = 10; % N
        string_length = 1; % m
        
        rho = total_mass / string_length;
        c = sqrt(tension_force / rho);
        wn_analytical = c*pi*wn_idx/string_length;
    
        % build matrices
        string_params = struct();
        string_params.n = num_masses;
        string_params.M = total_mass;
        string_params.Tf = tension_force;
        string_params.L = string_length;
    
        [M, K] = construct_2nd_order_matrices(string_params);
    
        % compute eigenvalues
        [~, lambda] = eig(K, M);
    
        % sort eigenvalues
        lambda_vals = sort(diag(lambda));
    
        % numerical 3rd mode frequency
        wn_numerical = sqrt(lambda_vals(wn_idx));
    
        err(n) = abs(wn_analytical - wn_numerical);
    end
    
    semilogy(1:n, err, '.', 'MarkerSize', 10);
    xlabel('Number of Masses');
    ylabel('Absolute Error');
    title('Numerical Error of Mode 1 Natural Frequency');
    grid on;

end