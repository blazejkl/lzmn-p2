function plotResults(resultsTable, histories)
%PLOTRESULTS Tworzy czytelne wykresy porownujace metody.
%
% Wykresy:
%   1. Residual wzgledny klasycznego GS tylko dla przypadkow zbie znych.
%   2. Residual wzgledny klasycznego GS dla przypadkow niezbie znych.
%   3. Liczba iteracji wzgledem rho = max |s_i/c_i|^2.
%   4. Czas wykonania metod.
%   5. Wzgledny blad rozwiazania.

outDir = fullfile(pwd, 'figures');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

methodNames = string(resultsTable.MethodName);
idxClassic = methodNames == "Classical Gauss-Seidel";
idxBlock = methodNames == "Block Gauss-Seidel 2x2";

%% 1. Residual wzgledny - przypadki zbie zne

figure;
hold on;

legendLabels = strings(0);

for i = 1:numel(histories)
    testName = string(histories(i).testName);

    idxRow = idxClassic & string(resultsTable.TestName) == testName;

    if any(idxRow) && resultsTable.Converged(idxRow)
        rh = histories(i).classicResidual;

        if ~isempty(rh) && all(isfinite(rh))
            rhRel = rh ./ max(rh(1), eps);
            semilogy(0:numel(rhRel)-1, max(rhRel, realmin), ...
                'LineWidth', 1.4);
            legendLabels(end+1) = testName; %#ok<AGROW>
        end
    end
end

xlabel('Iteracja');
ylabel('Residual wzgledny wzgledem startu');
title('Klasyczny Gauss-Seidel - przypadki zbie zne');
grid on;
set(gca, 'YScale', 'log');

if ~isempty(legendLabels)
    legend(legendLabels, 'Interpreter', 'none', 'Location', 'best');
end

saveFigureSafe(gcf, fullfile(outDir, 'residual_convergent.png'));

%% 2. Residual wzgledny - przypadki niezbie zne / problematyczne

figure;
hold on;

legendLabels = strings(0);

for i = 1:numel(histories)
    testName = string(histories(i).testName);

    idxRow = idxClassic & string(resultsTable.TestName) == testName;

    if any(idxRow) && ~resultsTable.Converged(idxRow)
        rh = histories(i).classicResidual;

        if ~isempty(rh)
            rh = rh(isfinite(rh));

            if ~isempty(rh)
                rhRel = rh ./ max(rh(1), eps);
                semilogy(0:numel(rhRel)-1, max(rhRel, realmin), ...
                    'LineWidth', 1.4);
                legendLabels(end+1) = testName; %#ok<AGROW>
            end
        end
    end
end

xlabel('Iteracja');
ylabel('Residual wzgledny wzgledem startu');
title('Klasyczny Gauss-Seidel - przypadki niezbie zne / niestabilne');
grid on;
set(gca, 'YScale', 'log');

if ~isempty(legendLabels)
    legend(legendLabels, 'Interpreter', 'none', 'Location', 'best');
end

saveFigureSafe(gcf, fullfile(outDir, 'residual_nonconvergent.png'));

%% 3. Liczba iteracji wzgledem rho

rho = resultsTable.RhoTheory(idxClassic);
iters = resultsTable.Iterations(idxClassic);
conv = resultsTable.Converged(idxClassic);

figure;
hold on;

plot(rho(conv), iters(conv), 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(rho(~conv), iters(~conv), 'x', 'MarkerSize', 10, 'LineWidth', 2);

xline(1, '--', 'rho = 1', 'LineWidth', 1.2);

xlabel('\rho = max_i |s_i/c_i|^2');
ylabel('Liczba iteracji');
title('Wplyw \rho na liczbe iteracji klasycznego GS');
legend('Zbiezne', 'Niezbie zne / zatrzymane', 'Granica \rho=1', ...
    'Location', 'best');
grid on;

saveFigureSafe(gcf, fullfile(outDir, 'iterations_vs_rho.png'));

%% 4. Czas wykonania metod

testNames = unique(string(resultsTable.TestName), 'stable');
numTests = numel(testNames);

timeClassic = NaN(numTests, 1);
timeBlock = NaN(numTests, 1);

for i = 1:numTests
    idxTest = string(resultsTable.TestName) == testNames(i);

    idxC = idxTest & idxClassic;
    idxB = idxTest & idxBlock;

    if any(idxC)
        timeClassic(i) = resultsTable.ElapsedTime(idxC);
    end

    if any(idxB)
        timeBlock(i) = resultsTable.ElapsedTime(idxB);
    end
end

figure;
bar([timeClassic timeBlock]);

xlabel('Numer testu');
ylabel('Czas [s]');
title('Porownanie czasu wykonania metod');
legend('Klasyczny GS', 'Blokowy GS 2x2', 'Location', 'best');
grid on;

xticks(1:numTests);
xticklabels("T" + string(1:numTests));

saveFigureSafe(gcf, fullfile(outDir, 'time_vs_test.png'));

%% 5. Wzgledny blad rozwiazania

errClassic = NaN(numTests, 1);
errBlock = NaN(numTests, 1);

for i = 1:numTests
    idxTest = string(resultsTable.TestName) == testNames(i);

    idxC = idxTest & idxClassic;
    idxB = idxTest & idxBlock;

    if any(idxC)
        errClassic(i) = resultsTable.RelativeSolutionError(idxC);
    end

    if any(idxB)
        errBlock(i) = resultsTable.RelativeSolutionError(idxB);
    end
end

figure;
hold on;

semilogy(1:numTests, max(errClassic, realmin), 'o-', ...
    'LineWidth', 1.4, 'MarkerSize', 7);
semilogy(1:numTests, max(errBlock, realmin), 's-', ...
    'LineWidth', 1.4, 'MarkerSize', 7);

xlabel('Numer testu');
ylabel('Wzgledny blad rozwiazania');
title('Porownanie wzglednego bledu rozwiazania');
legend('Klasyczny GS', 'Blokowy GS 2x2', 'Location', 'best');
grid on;
set(gca, 'YScale', 'log');

xticks(1:numTests);
xticklabels("T" + string(1:numTests));

saveFigureSafe(gcf, fullfile(outDir, 'solution_error_comparison.png'));

end


function saveFigureSafe(figHandle, filePath)
%SAVEFIGURESAFE Bezpieczny zapis wykresu do pliku PNG.

try
    exportgraphics(figHandle, filePath, 'Resolution', 200);
catch
    saveas(figHandle, filePath);
end

end