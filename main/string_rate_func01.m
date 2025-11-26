% INPUTS
% t: current time
% V: system state. V = [U;dUdt] where
% U and dUdt are n x 1 column vectors
% string_params: a struct containing the system parameters describing the string
% string_params.n: number of masses
% string_params.M: total mass attached to the string
% string_params.Uf_func: function describing motion of end point
% string_params.dUfdt_func: time derivative of Uf
% string_params.Tf: tension in string
% string_params.L: length of string
% string_params.c: damping coefficient
% string_params.dx: horizontal spacing between masses
function dVdt = string_rate_func01(t, V, string_params)
    n = string_params.n; %number of masses
    M = string_params.M; %total mass attached to the string
    Uf_func = string_params.Uf_func; %function describing motion of end point
    dUfdt_func = string_params.dUfdt_func; %time derivative of Uf
    Tf = string_params.Tf; %tension in string
    L = string_params.L; %length of string
    c = string_params.c; %damping coefficient
    dx = string_params.dx; %horizontal spacing between masses

    %unpack state variable
    U = V(1:n);
    dUdt = V((n+1):(2*n));
    Uf = Uf_func(t);
    dUfdt = dUfdt_func(t);

    % construct useful data structures for computation
    m = M/n; % mass of each individual ball
    
    % construct K (stiffness) matrix
    I_n = eye(n);
    K = 2 * I_n;
    K = K + circshift(I_n, [0, 1]) * -1;
    K = K + circshift(I_n, [0, -1]) * -1;
    K(1, end) = 0; % delete unwanted value in top right corner
    K(end, 1) = 0; % delete unwanted value in bottom right corner
    K = Tf/dx * K;
    
    % construct input matrix
    B = zeros(n, 1);
    B(end) = 1;
    B = Tf/dx * B;

    %compute acceleration
    d2Udt2 =  1/m * (-K*U + B*Uf);

    %assemble state derivative
    dVdt = [dUdt; d2Udt2];
end