clear;
clc;
close all;

% Projekt: klasyczna i blokowa metoda Gaussa-Seidela
% dla ukladu A x = b, gdzie A = [C S; -S C].

tol = 1e-10;
maxIter = 10000;
testIDs = 1:12;

metricsList = struct([]);
histories = struct([]);

fprintf('START TESTOW\n');
fprintf('===========\n\n');

for idx = 1:numel(testIDs)

    testID = testIDs(idx);
    problem = generateProblem(testID);

    fprintf('Test %d: %s\n', testID, problem.testName);
    fprintf('Opis: %s\n', problem.comment);
    fprintf('p = %d, n = %d, rho = %.4e\n', ...
        problem.p, problem.n, problem.rhoTheory);

    x0 = zeros(problem.n, 1);

    % 1. Klasyczna metoda Gaussa-Seidela
    [xClassic, infoClassic] = gaussSeidelClassic( ...
        problem.c, problem.s, problem.b, x0, tol, maxIter);

    % 2. Blokowa metoda Gaussa-Seidela 2x2
    [xBlock, infoBlock] = gaussSeidelBlock2x2( ...
        problem.c, problem.s, problem.b);

    % 3. Rozwiazanie referencyjne MATLAB-a
    xMatlab = problem.A \ problem.b;
    matlabError = norm(xMatlab - problem.x_exact, 2) / ...
        max(1, norm(problem.x_exact, 2));

    fprintf('Blad MATLAB A\\b wzgledem x_exact: %.4e\n', matlabError);

    % 4. Metryki
    metricsClassic = computeMetrics(problem, xClassic, infoClassic);
    metricsBlock = computeMetrics(problem, xBlock, infoBlock);

    if isempty(metricsList)
        metricsList = metricsClassic;
        metricsList(2) = metricsBlock;
    else
        metricsList(end + 1) = metricsClassic; %#ok<SAGROW>
        metricsList(end + 1) = metricsBlock; %#ok<SAGROW>
    end

    % 5. Historie residuali do wykresow
    histories(idx).testName = string(problem.testName); %#ok<SAGROW>
    histories(idx).rhoTheory = problem.rhoTheory;
    histories(idx).n = problem.n;
    histories(idx).classicResidual = infoClassic.residualHistory;
    histories(idx).blockResidual = infoBlock.residualHistory;

    fprintf('Klasyczny GS: converged = %d, iter = %d, residual = %.4e\n', ...
        infoClassic.converged, infoClassic.iterations, infoClassic.finalResidual);

    fprintf('Blokowy GS:   converged = %d, iter = %d, residual = %.4e\n', ...
        infoBlock.converged, infoBlock.iterations, infoBlock.finalResidual);

    fprintf('\n');

end

% 6. Tabela wynikow
resultsTable = struct2table(metricsList);

disp('TABELA WYNIKOW:');
disp(resultsTable);

% 7. Zapis tabeli do CSV
writetable(resultsTable, 'results_table.csv');

% 8. Wykresy
plotResults(resultsTable, histories);

fprintf('\nZAKONCZONO.\n');
fprintf('Zapisano tabele: results_table.csv\n');
fprintf('Zapisano wykresy w folderze: figures\n');