L = 0.35;
w = 0.02;
t = 0.002;
rho = 7850;
E = 200e9;
m = 0.02;

[natFreqs, modeShapes] = RayleighRitz(L, w, t, rho, E, m);

%[natFreqs, modeShapes] = FiniteElement();