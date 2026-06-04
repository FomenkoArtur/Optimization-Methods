clear; clc; close all;
format long

fprintf('\n===============================================================\n');
fprintf('  ВИЗУАЛЬНАЯ ДИАГНОСТИКА: w < 0 в методе Полака-Рибьера\n');
fprintf('===============================================================\n');

%% ── Тест 1: Функция Розенброка ────────────────────────────────────────
a = 1; b = 100;
syms z [2 1]
f_rosen = (a - z(1))^2 + b * (z(2) - z(1)^2)^2;

fprintf('\n🔹 ФУНКЦИЯ РОЗЕНБРОКА:\n');
fprintf('   f(x) = (%.0f-x₁)² + %.0f·(x₂-x₁²)²\n', a, b);
fprintf('   Минимум: (%.1f, %.1f), f* = 0\n\n', a, a^2);

x0 = [-100.5, -100.5]; E = 1e-4;
% Вызов с флагом 'final' для отрисовки только последних случаев
[x_min, f_min, ~, n, w_log, viz_rosen] = PolakRibiereViz(x0, E, f_rosen, 'final');

fprintf('   Результат: x = [%.6f, %.6f], f = %.6e, итераций: %d\n', ...
        x_min(1), x_min(2), f_min, n);
fprintf('   Случаев w < 0: %d\n', length(w_log));

if ~isempty(w_log)
    fprintf('\n   Последние 3 случая w < 0:\n');
    fprintf('   %-6s | %-12s | %-12s | %-12s | %-10s\n', ...
            'Iter', 'w_raw', '||grad||', '||s_prev||', 'cos(grad,s)');
    fprintf('   ----------------------------------------------------------\n');
    for i = max(1, length(w_log)-2):length(w_log)
        grad_x = w_log(i, 2:3);
        s_prev = w_log(i, 6:7);
        cos_angle = (grad_x * s_prev') / (norm(grad_x)*norm(s_prev) + eps);
        fprintf('   %-6d | %-12.3e | %-12.3e | %-12.3e | %-10.3f\n', ...
                w_log(i,1), w_log(i,8), norm(grad_x), norm(s_prev), cos_angle);
    end
end

%% ── Тест 2: Функция Бута ─────────────────────────────────────────────
fprintf('\n🔹 ФУНКЦИЯ БУТА:\n');
f_booth = (z(1) + 2*z(2) - 7)^2 + (2*z(1) + z(2) - 5)^2;
fprintf('   f(x) = (x₁+2x₂-7)² + (2x₁+x₂-5)²\n');
fprintf('   Минимум: (1, 3), f* = 0\n\n');

x0 = [0, 0]; E = 1e-4;
[x_min_b, f_min_b, ~, n_b, w_log_b, viz_booth] = PolakRibiereViz(x0, E, f_booth, 'final');

fprintf('   Результат: x = [%.6f, %.6f], f = %.6e, итераций: %d\n', ...
        x_min_b(1), x_min_b(2), f_min_b, n_b);
fprintf('   Случаев w < 0: %d\n', length(w_log_b));

if isempty(w_log_b)
    fprintf('Ни одного случая w < 0 — метод работает стабильно.\n');
end
