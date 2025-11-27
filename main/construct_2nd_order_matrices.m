% build the mass and stiffness matrices that describe the 2nd order system.
% 
% INPUTS
% string_params: a struct containing the system parameters describing the string
% string_params.n: number of masses
% string_params.M: total mass attached to the string
% string_params.Uf_func: function describing motion of end point
% string_params.dUfdt_func: time derivative of Uf
% string_params.Tf: %tension in string
% string_params.L: %length of string
% string_params.c: %damping coefficient
% string_params.dx: %horizontal spacing between masses
% 
% OUTPUTS
% M_mat: the n x n mass (inertia) matrix
% K_mat: the n x n stiffness matrix
function [M_mat, K_mat] = construct_2nd_order_matrices(string_params)
    % construct diagonal M matrix
    M_mat = string_params.M / string_params.n * eye(string_params.n);
    
    % construct K matrix
    I_n = eye(string_params.n);
    K_mat = 2*I_n - circshift(I_n, [0, 1]) - circshift(I_n, [0, -1]);
    K_mat(1, end) = 0; % remove unwanted 1 from top right
    K_mat(end, 1) = 0; % remove unwanted 1 from bottom left
end