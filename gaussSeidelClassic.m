function [x, info] = gaussSeidelClassic(c, s, b, x0, tol, maxIter)
%GAUSSSEIDELCLASSIC Klasyczna metoda Gaussa-Seidela dla A x = b.
%
% Rozwazamy macierz:
%       A = [ C   S
%            -S   C ],
% gdzie C = diag(c), S = diag(s).
%
% Dla x = [u; v] oraz b = [f; g] uklad ma postac:
%       C u + S v = f,
%      -S u + C v = g.
%
% Klasyczne iteracje Gaussa-Seidela w tym uporzadkowaniu:
%       u_new = (f - s .* v_old) ./ c,
%       v_new = (g + s .* u_new) ./ c.
%
% Wejscie:
%   c, s    - wektory diagonalne macierzy C i S
%   b       - prawa strona ukladu, wektor dlugosci 2p
%   x0      - przyblizenie poczatkowe, wektor dlugosci 2p
%   tol     - tolerancja zatrzymania
%   maxIter - maksymalna liczba iteracji
%
% Wyjscie:
%   x    - przyblizone rozwiazanie
%   info - struktura z informacjami o przebiegu metody

timerStart = tic;

c = c(:);
s = s(:);
b = b(:);
x0 = x0(:);

p = length(c);
n = 2 * p;

if length(s) ~= p
    error('Wektory c i s musza miec taka sama dlugosc.');
end

if length(b) ~= n
    error('Wektor b musi miec dlugosc 2p.');
end

if length(x0) ~= n
    error('Wektor x0 musi miec dlugosc 2p.');
end

if any(abs(c) <= eps)
    error('Nie mozna wykonac metody klasycznej: wystepuje c_i zbyt bliskie zeru.');
end

f = b(1:p);
g = b(p+1:end);

u = x0(1:p);
v = x0(p+1:end);

bNorm = max(norm(b, 2), eps);
rhoTheory = max(abs(s ./ c).^2);

residualHistory = zeros(maxIter + 1, 1);
relativeResidualHistory = zeros(maxIter + 1, 1);
stepHistory = zeros(maxIter + 1, 1);

rTop = c .* u + s .* v - f;
rBottom = -s .* u + c .* v - g;
residual = norm([rTop; rBottom], 2);
relativeResidual = residual / bNorm;

residualHistory(1) = residual;
relativeResidualHistory(1) = relativeResidual;
stepHistory(1) = NaN;

converged = false;
stopReason = "Osiagnieto maksymalna liczbe iteracji";
iterations = 0;

initialResidual = max(residual, eps);

if relativeResidual < tol
    converged = true;
    stopReason = "Warunek residualu spelniony dla x0";
else
    for k = 1:maxIter
        uOld = u;
        vOld = v;

        % Klasyczny Gauss-Seidel:
        % najpierw aktualizujemy u z uzyciem starego v,
        % nastepnie aktualizujemy v z uzyciem nowego u.
        u = (f - s .* vOld) ./ c;
        v = (g + s .* u) ./ c;

        step = norm([u - uOld; v - vOld], 2) / ...
            max(1, norm([u; v], 2));

        rTop = c .* u + s .* v - f;
        rBottom = -s .* u + c .* v - g;
        residual = norm([rTop; rBottom], 2);
        relativeResidual = residual / bNorm;

        residualHistory(k + 1) = residual;
        relativeResidualHistory(k + 1) = relativeResidual;
        stepHistory(k + 1) = step;

        iterations = k;

        if isnan(residual) || isinf(residual)
            converged = false;
            stopReason = "Residual stal sie NaN albo Inf";
            break;
        end

        if residual > max(1e12, 1e12 * initialResidual)
            converged = false;
            stopReason = "Wykryto rozbieznosc residualu";
            break;
        end

        if relativeResidual < tol
            converged = true;
            stopReason = "Spelniono warunek na residual wzgledny";
            break;
        end

        if step < tol && relativeResidual < sqrt(tol)
            converged = true;
            stopReason = "Spelniono warunek malego kroku i malego residualu";
            break;
        end

        if step < tol && relativeResidual >= sqrt(tol)
            converged = false;
            stopReason = "Stagnacja: krok maly, ale residual nadal za duzy";
            break;
        end
    end
end

residualHistory = residualHistory(1:iterations + 1);
relativeResidualHistory = relativeResidualHistory(1:iterations + 1);
stepHistory = stepHistory(1:iterations + 1);

x = [u; v];

info.methodName = "Classical Gauss-Seidel";
info.converged = converged;
info.iterations = iterations;
info.residualHistory = residualHistory;
info.relativeResidualHistory = relativeResidualHistory;
info.stepHistory = stepHistory;
info.finalResidual = residualHistory(end);
info.finalRelativeResidual = relativeResidualHistory(end);
info.elapsedTime = toc(timerStart);
info.stopReason = stopReason;
info.rhoTheory = rhoTheory;
end