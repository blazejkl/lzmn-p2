function problem = generateProblem(testID)
%GENERATEPROBLEM Generuje dane testowe dla ukladu A x = b.
%
% Macierz ma postac:
%       A = [ C   S
%            -S   C ],
% gdzie C = diag(c), S = diag(s), c_i^2 + s_i^2 = 1.
%
% Wejscie:
%   testID - numer testu od 1 do 6
%
% Wyjscie:
%   problem - struktura zawierajaca A, b, x_exact, c, s oraz opis testu

rng(testID);

switch testID
    case 1
        problem.testName = "Maly latwy przypadek";
        problem.comment = "Maly wymiar, male s/c, szybka zbieznosc klasycznego GS.";
        p = 3;
        theta = [0.05; 0.10; 0.20];

    case 2
        problem.testName = "Sredni losowy zbie zny";
        problem.comment = "Sredni wymiar, losowe theta ponizej pi/4, klasyczny GS powinien zbiegac.";
        p = 50;
        theta = 0.05 + (0.60 - 0.05) * rand(p, 1);

    case 3
        problem.testName = "Duzy test wydajnosciowy";
        problem.comment = "Duzy wymiar, test czasu i skalowania metod.";
        p = 1000;
        theta = 0.05 + (0.50 - 0.05) * rand(p, 1);

    case 4
        problem.testName = "Przypadek blisko granicy zbieznosci";
        problem.comment = "Theta blisko pi/4 od dolu, klasyczny GS powinien zbiegac wolno.";
        p = 50;
        theta = linspace(0.72, 0.78, p)';

    case 5
        problem.testName = "Przypadek rozbiezny klasycznego GS";
        problem.comment = "Theta powyzej pi/4, max |s/c|^2 > 1, klasyczny GS powinien zawiesc.";
        p = 20;
        theta = linspace(0.86, 1.00, p)';

    case 6
        problem.testName = "Prawie osobliwe C";
        problem.comment = "c_i bardzo male, det(C) formalnie niezerowy, ale klasyczny GS jest niestabilny.";
        p = 20;
        theta = (pi/2 - 1e-8) * ones(p, 1);

    otherwise
        error('Nieznany testID. Dozwolone wartosci: 1, 2, ..., 6.');
end

c = cos(theta);
s = sin(theta);

C = spdiags(c, 0, p, p);
S = spdiags(s, 0, p, p);

A = [C S; -S C];

n = 2 * p;

x_exact = randn(n, 1);
b = A * x_exact;

f = b(1:p);
g = b(p+1:end);

rhoTheory = max(abs(s ./ c).^2);

problem.p = p;
problem.n = n;
problem.c = c;
problem.s = s;
problem.C = C;
problem.S = S;
problem.A = A;
problem.x_exact = x_exact;
problem.b = b;
problem.f = f;
problem.g = g;
problem.rhoTheory = rhoTheory;
end