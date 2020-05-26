function [natFreqs, modeShapes] = RayleighRitz(L, w, t, rho, E, M_t)

% Calculating cross-sectional area and second moment of area
A_cs = w*t;
I = w*(t^3)/12;

% Defining our trial functions
syms x
phi = [x^2, x^3, x^4, x^5];
n = length(phi);

% Pre-allocating
k = zeros(n);
m = zeros(n);
natFreqs = zeros(1, n);
modeShapes = zeros(n);

for i = 1:n
    for j = 1:n
        % Filling in the k matrix
        func_k = diff(phi(i), 2)*diff(phi(j), 2);
        k(i,j) = int(E*I*func_k, [0 L]);
        
        % Filling in the m matrix
        func_m = rho*A_cs*phi(i)*phi(j);
        x = L;
        m(i,j) = int(func_m, [0 L]) + M_t*subs(phi(i))*subs(phi(j));
    end
end

% Converting from sym to double
k = double(k);
m = double(m);

% We now have |k-lambda*M|=0. This is an eigen problem
[eigenVectors, eigenValues] = eig(k, m);

for i = 1:n
    % Find the natural frequencies
    natFreqs(i) = sqrt(eigenValues(i, i));
    % Find the mode shapes
    modeShapes(:,i) = eigenVectors(:,i)/max(abs(eigenVectors(:,i)));
end

end