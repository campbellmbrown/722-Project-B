function [natFreqs, modeShapes] = RayleighRitz(L, w, t, rho, E, M_t)

A_cs = w*t;
I = w*(t^3)/12;

syms x
phi = [x^2, x^3, x^4, x^5, x^6];
n = length(phi);

for i = 1:n
    for j = 1:n
        func_k = diff(phi(i), 2)*diff(phi(j), 2);
        k(i,j) = int(E*I*func_k, [0 L]);
        
        func_m = rho*A_cs*phi(i)*phi(j);
        x = L;
        m(i,j) = int(func_m, [0 L]) + M_t*subs(phi(i))*subs(phi(j));
    end
end

k = double(k);
m = double(m);

% We now have |k-lambda*M|=0. This is an eigenvalue problem.
[eigenVectors, eigenValues] = eig(k, m);

for i = 1:n
    natFreqs(i) = sqrt(eigenValues(i, i));
    modeShapes(:,i) = eigenVectors(:,i)/max(abs(eigenVectors(:,i)));
end

end