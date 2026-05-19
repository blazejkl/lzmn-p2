function metrics = computeMetrics(problem, x, info)
%COMPUTEMETRICS Oblicza metryki bledu i zbieznosci dla metody.
%
% Wejscie:
%   problem - struktura wygenerowana przez generateProblem
%   x       - rozwiazanie zwrocone przez metode
%   info    - struktura zwrocona przez metode
%
% Wyjscie:
%   metrics - struktura gotowa do zamiany na tabele

x = x(:);

solutionError = norm(x - problem.x_exact, 2);
relativeSolutionError = solutionError / max(1, norm(problem.x_exact, 2));

metrics.TestName = string(problem.testName);
metrics.MethodName = string(info.methodName);
metrics.p = problem.p;
metrics.n = problem.n;
metrics.RhoTheory = problem.rhoTheory;
metrics.Converged = logical(info.converged);
metrics.Iterations = info.iterations;
metrics.FinalResidual = info.finalResidual;
metrics.FinalRelativeResidual = info.finalRelativeResidual;
metrics.SolutionError = solutionError;
metrics.RelativeSolutionError = relativeSolutionError;
metrics.ElapsedTime = info.elapsedTime;
metrics.StopReason = string(info.stopReason);
metrics.Comment = string(problem.comment);
end