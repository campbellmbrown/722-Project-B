function [nat_freqs, mode_shapes, M, K] = FiniteElement(L, w, t, rho, E, M_t, elements)
% FiniteElement function
% Finds the finite element approximation for natural frequencies and mode
% shapes of the beam
% INPUTS ======================
% L - beam length (m)
% w - beam width (m)
% t - beam thickness (m)
% rho - beam density (kg/m^3)
% E - Young's modulus (pascals)
% M_t - tip point-mass (kg)
% elements - number of elements to split the beam into
% OUTPUTS =====================
% nat_freqs - natural frequencies (rad/s)
% mode_shapes - modes shapes
% M - mass matrix
% K - stiffness matrix
% AUTHOR ======================
% Campbell Brown - 738509729

A_cs = w*t;
I = w*(t^3)/12;
DOF_t = (1 + elements)*2;	% Total DOFs

% Find the element K and M matrices
[K_element, M_element] = BeamFEMatrices(rho, A_cs, E*I, L/elements);

% Pre-allocating
K = zeros(DOF_t);
M = zeros(DOF_t);

% Add the contribution of each element
for i = 1:elements
    ind_a = 1 + (i - 1)*2;  % Find the starting index
    ind_b = (i + 1)*2;      % Find the ending index
    ind = ind_a:ind_b;
    K(ind, ind) = K(ind, ind) + K_element;  % Add the element K
    M(ind, ind) = M(ind, ind) + M_element;  % Add the element M
end

% Displacement at first node is zero, so remove first row and column
M = M(2:DOF_t, 2:DOF_t);
K = K(2:DOF_t, 2:DOF_t);
DOF_t = DOF_t - 1;

% Add the point mass contribution to mass matrix
M(DOF_t-1, DOF_t-1) = M(DOF_t-1, DOF_t-1) + M_t;

% Solve the eigen-problem
[eigen_vectors, eigen_values] = eig(K, M);

% Get the natural frequencies (diagonal terms)
nat_freqs = sqrt(diag(eigen_values));
nat_freqs = nat_freqs(1:4);

% Finding the mode shapes
for i = 1:4
    % Getting every second value (relating to the displacement DOFs)
    ind = 2:2:(DOF_t-1);
    eigen_vector = eigen_vectors(ind,i);
    % Normalise the mode shape to have a max value of 1
    normalise_value = max(abs(eigen_vector));
    mode_shapes(1:length(ind),i) = eigen_vector/normalise_value;
end

end