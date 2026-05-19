function [x, info] = gaussSeidelBlock2x2(c, s, b)
%GAUSSSEIDELBLOCK2X2 Blokowa metoda Gaussa-Seidela z blokami 2x2.
%
% Wykorzystujemy naturalne bloki odpowiadajace parom niewiadomych (u_i, v_i).
% Dla kazdego i rozwiazujemy dokladnie maly uklad:
%
%       [ c_i   s_i ] [u_i] = [f_i]
%       [-s_i   c_i ] [v_i]   [g_i]
%
% Poniewaz c_i^2 + s_i^2 = 1, macierz odwrotna bloku to:
%
%       [ c_i  -s_i ]
%       [ s_i   c_i ]
%
% Stad:
%       u_i = c_i f_i - s_i g_i,
%       v_i = s_i f_i + c_i g_i.
%
% Dla tej szczegolnej macierzy metoda blokowa rozwiazuje uklad
% praktycznie w jednym przejsciu.

timerStart = tic;

c = c(:);
s = s(:);
b = b(:);

p = length(c);
n = 2 * p;

if length(s) ~= p
    error('Wektory c i s musza miec taka sama dlugosc.');
end

if length(b) ~= n
    error('Wektor b musi miec dlugosc 2p.');
end

f = b(1:p);
g = b(p+1:end);

u = c .* f - s .* g;
v = s .* f + c .* g;

x = [u; v];

rTop = c .* u + s .* v - f;
rBottom = -s .* u + c .* v - g;
residual = norm([rTop; rBottom], 2);

bNorm = max(norm(b, 2), eps);
relativeResidual = residual / bNorm;

info.methodName = "Block Gauss-Seidel 2x2";
info.converged = true;
info.iterations = 1;
info.residualHistory = residual;
info.relativeResidualHistory = relativeResidual;
info.stepHistory = NaN;
info.finalResidual = residual;
info.finalRelativeResidual = relativeResidual;
info.elapsedTime = toc(timerStart);
info.stopReason = "Rozwiazano przez dokladne bloki 2x2";
info.rhoTheory = NaN;
end