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
    dx = L / (n-1); %horizontal spacing between masses
    %dx = string_params.dx; %horizontal spacing between masses

    %unpack state variable
    U = V(1:n);
    dUdt = V((n+1):(2*n));
    Uf = Uf_func(t);
    dUfdt = dUfdt_func(t);

    % compute acceleration
    d2Udt2 = zeros(n, 1);
    % compute acceleration for each mass
    for i=2:n-1
        d2Udt2(i) = n/M * (Tf/dx * (U(i-1) - 2*U(i) + U(i+1)) + c/dx * (dUdt(i-1) - 2*dUdt(i) + dUdt(i+1)));
    end
    % apply boundary conditions
    d2Udt2(1) = n/M * (Tf/dx * (-2*U(1) + U(2)) + c/dx * (-2*dUdt(1) + dUdt(2)));
    %d2Udt2(1) = -n/M * (Tf/dx * (U(n-1) - 2*U(n) + Uf) + c/dx * (dUdt(n-1) - 2*dUdt(n) + dUfdt));
    d2Udt2(end) = n/M * (Tf/dx * (U(n-1) - 2*U(n) + Uf) + c/dx * (dUdt(n-1) - 2*dUdt(n) + dUfdt));

    %assemble state derivative
    dVdt = [dUdt; d2Udt2];
end