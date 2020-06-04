function [nat_freqs, mode_shapes] = RayleighRitz(L, w, t, rho, E, M_t)
% RayleighRitz function
% Finds the Rayleigh-Ritz approximation for natural frequencies and mode
% shapes of the beam
% INPUTS ======================
% L - beam length (m)
% w - beam width (m)
% t - beam thickness (m)
% rho - beam density (kg/m^3)
% E - Young's modulus (pascals)
% M_t - tip point-mass (kg)
% OUTPUTS =====================
% nat_freqs - natural frequencies (rad/s)
% mode_shapes - modes shapes
% AUTHOR ======================
% Campbell Brown - 738509729

% Calculating cross-sectional area and second moment of area
A_cs = w*t;
I = w*(t^3)/12;

% Defining our trial functions
syms x
phi = [x, x^2, x^3, x^4];

% Pre-allocating
K = zeros(4);
M = zeros(4);

for i = 1:4
    for j = 1:4
        % Filling in the k matrix
        func_k = diff(phi(i), 2)*diff(phi(j), 2);
        K(i,j) = int(E*I*func_k, [0 L]);
        % Filling in the m matrix
        func_m = rho*A_cs*phi(i)*phi(j);
        M(i,j) = int(func_m, [0 L]) + M_t*subs(phi(i), x, L)*subs(phi(j), x, L);
    end
end

% Converting from sym to double
K = double(K);
M = double(M);

K2 = E*I*[...
    0,      0,      0,      0;...
    0,      4*L,    6*L^2,  8*L^3;...
    0,      6*L^2,  12*L^3, 18*L^4;...
    0,      8*L^3,  18*L^4, 144/5*L^5];
M2 = rho*A_cs*[...
    L^3/3,  L^4/4,  L^5/5,  L^6/6;...
    L^4/4,  L^5/5,  L^6/6,  L^7/7;...
    L^5/5,  L^6/6,  L^7/7,  L^8/8;...
    L^6/6,  L^7/7,  L^8/8,  L^9/9];
M2 = M2 + M_t*[...
    L^2,    L^3,    L^4,    L^5;...
    L^3,    L^4,    L^5,    L^6;...
    L^4,    L^5,    L^6,    L^7;...
    L^5,    L^6,    L^7,    L^8];

% We now have |k-lambda*M|=0. This is an eigen problem
[eigen_vectors, eigen_values] = eig(K, M);
% Find the natural frequencies
nat_freqs = sqrt(diag(eigen_values));

% Find the mode shapes
x_range = 0:0.01:L;
mode_shapes = zeros(4, length(x_range));
for i = 1:4
    for j = 1:4
        mode_shapes(i,:) = mode_shapes(i,:) + eigen_vectors(j,i)*double(subs(phi(j), x, x_range));
    end
    % Normalising to have max value of 1
    mode_shapes(i,:) = mode_shapes(i,:)/max(abs(mode_shapes(i,:)));
end

end