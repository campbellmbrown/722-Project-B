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

k_2(1, 1) = 4*L;
k_2(1, 2) = 6*L^2;
k_2(1, 3) = 8*L^3;
k_2(1, 4) = 10*L^4;
k_2(2, 1) = 6*L^2;
k_2(2, 2) = 12*L^3;
k_2(2, 3) = 18*L^4;
k_2(2, 4) = 24*L^5;
k_2(3, 1) = 8*L^3;
k_2(3, 2) = 18*L^4;
k_2(3, 3) = 144/5*L^5;
k_2(3, 4) = 40*L^6;
k_2(4, 1) = 10*L^4;
k_2(4, 2) = 24*L^5;
k_2(4, 3) = 40*L^6;
k_2(4, 4) = 400/7*L^7;
k_2 = k_2*E*I;

m_2(1, 1) = rho*A_cs*L^5/5 + M_t*L^4;
m_2(1, 2) = rho*A_cs*L^6/6 + M_t*L^5;
m_2(1, 3) = rho*A_cs*L^7/7 + M_t*L^6;
m_2(1, 4) = rho*A_cs*L^8/8 + M_t*L^7;
m_2(2, 1) = rho*A_cs*L^6/6 + M_t*L^5;
m_2(2, 2) = rho*A_cs*L^7/7 + M_t*L^6;
m_2(2, 3) = rho*A_cs*L^8/8 + M_t*L^7;
m_2(2, 4) = rho*A_cs*L^9/9 + M_t*L^8;
m_2(3, 1) = rho*A_cs*L^7/7 + M_t*L^6;
m_2(3, 2) = rho*A_cs*L^8/8 + M_t*L^7;
m_2(3, 3) = rho*A_cs*L^9/9 + M_t*L^8;
m_2(3, 4) = rho*A_cs*L^10/10 + M_t*L^9;
m_2(4, 1) = rho*A_cs*L^8/8 + M_t*L^7;
m_2(4, 2) = rho*A_cs*L^9/9 + M_t*L^8;
m_2(4, 3) = rho*A_cs*L^10/10 + M_t*L^9;
m_2(4, 4) = rho*A_cs*L^11/11 + M_t*L^10;
end