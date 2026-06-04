clear; clc; format long

a_roz = 1; b_roz = 100;
a_ack = 20; b_ack = 0.2; c_ack = 2*pi;

% УБРАНО: syms z [2 1] и все f_sym_* (символьные выражения)

f_vec_roz  = @(x) (a_roz-x(1))^2 + b_roz*(x(2)-x(1)^2)^2;
f_vec_ack  = @(x) -a_ack*exp(-b_ack*sqrt((x(1)^2+x(2)^2)/2)) ...
                  - exp((cos(c_ack*x(1))+cos(c_ack*x(2)))/2) + a_ack + exp(1);
f_vec_himm = @(x) (x(1)^2+x(2)-11)^2 + (x(1)+x(2)^2-7)^2;
f_vec_boot = @(x) (x(1)+2*x(2)-7)^2 + (2*x(1)+x(2)-5)^2;
f_vec_rast = @(x) 20 + (x(1)^2-10*cos(2*pi*x(1))) + (x(2)^2-10*cos(2*pi*x(2)));

% Сетки для отрисовки поверхности
grids = {
    [-2.5,  2.5,  -1.5,  3.5,  300];   % Rosenbrock
    [-5.0,  5.0,  -5.0,  5.0,  250];   % Ackley
    [-5.0,  5.0,  -5.0,  5.0,  250];   % Himmelblau
    [-5.0,  5.0,  -5.0,  5.0,  250];   % Booth
    [-5.5,  5.5,  -5.5,  5.5,  250];   % Rastrigin
};

% Начальные точки — далеко от минимума, но внутри сетки (чтобы видеть траектории)
tests = {                               % УБРАНО: столбец f_sym из каждой строки
    f_vec_roz,  [-100.0,  100.0], 1e-4, ...
        'Rozenbrock',  'f=(1-x1)^2+100*(x2-x1^2)^2',                [1.0, 1.0], 0;
    f_vec_ack,  [-1.5,  3.0], 1e-4, ...
        'Ackley',      'f=-20*exp(...)-exp(...)+20+e',                [0.0, 0.0], 0;
    f_vec_himm, [0.0, -4.0], 1e-4, ...
        'Himmelblau',  'f=(x1^2+x2-11)^2+(x1+x2^2-7)^2',            [3.0, 2.0], 0;
    f_vec_boot, [-4.5, -4.0], 1e-4, ...
        'Booth',       'f=(x1+2x2-7)^2+(2x1+x2-5)^2',               [1.0, 3.0], 0;
    f_vec_rast, [-4.5,  4.5], 1e-4, ...
        'Rastrigin',   'f=20+x1^2-10cos(2pi*x1)+x2^2-10cos(2pi*x2)',[0.0, 0.0], 0;
};

method_names = {'SteepestDescent','ConjugateGradient','PolakRibiere'};
colors  = {[0.85,0.10,0.10]; [0.05,0.65,0.10]; [0.10,0.30,0.95]};
markers = {'o','s','^'};

fprintf('\n');
fprintf('==========================================================================\n');
fprintf('        СРАВНЕНИЕ МЕТОДОВ МИНИМИЗАЦИИ ФНП (методы первого порядка)\n');
fprintf('==========================================================================\n');

for i = 1:size(tests,1)
    f_vec  = tests{i,1};           % УБРАНО: f_sym = tests{i,2}; сдвинуты индексы столбцов
    x0     = tests{i,2};
    E      = tests{i,3};
    fname  = tests{i,4};
    fdesc  = tests{i,5};
    x_opt  = tests{i,6};
    f_opt  = tests{i,7};
    g      = grids{i};

    [x_sd, f_sd, t_sd, n_sd, hist_sd] = SteepestDescentMethod (x0, E, f_vec); % УБРАНО: f_sym
    [x_cg, f_cg, t_cg, n_cg, hist_cg] = ConjugateGradientMethod(x0, E, f_vec); % УБРАНО: f_sym
    [x_pr, f_pr, t_pr, n_pr, hist_pr] = PolakRibiereMethod     (x0, E, f_vec); % УБРАНО: f_sym

    fprintf('\n');
    fprintf('--------------------------------------------------------------------------\n');
    fprintf('  ФУНКЦИЯ : %s\n', fname);
    fprintf('  ФОРМУЛА : %s\n', fdesc);
    fprintf('  Начальная точка : x0=[%.2f, %.2f] | Точность: %.0e\n', x0(1), x0(2), E);
    fprintf('--------------------------------------------------------------------------\n');
    fprintf('  %-22s | x1=%-12.6f x2=%-12.6f | F=%-14.6e | T=%-11.7f | N=%-6d\n', ...
        'SteepestDescent',   x_sd(1), x_sd(2), f_sd, t_sd, n_sd);
    fprintf('  %-22s | x1=%-12.6f x2=%-12.6f | F=%-14.6e | T=%-11.7f | N=%-6d\n', ...
        'ConjugateGradient', x_cg(1), x_cg(2), f_cg, t_cg, n_cg);
    fprintf('  %-22s | x1=%-12.6f x2=%-12.6f | F=%-14.6e | T=%-11.7f | N=%-6d\n', ...
        'PolakRibiere',      x_pr(1), x_pr(2), f_pr, t_pr, n_pr);
    fprintf('--------------------------------------------------------------------------\n');

    if strcmp(fname, 'Himmelblau')
        fprintf('  Аналитические минимумы (4 точки):\n');
        fprintf('    x1* = [ 3.0000,  2.0000]\n');
        fprintf('    x2* = [-2.8051,  3.1313]\n');
        fprintf('    x3* = [-3.7793, -3.2832]\n');
        fprintf('    x4* = [ 3.5844, -1.8481]\n');
        fprintf('    f*  = 0\n');
    else
        fprintf('  Аналитический минимум: x*=[%.4f, %.4f], f*=%.4f\n', x_opt(1), x_opt(2), f_opt);
    end

    % ── Построение фигуры ──────────────────────────────────────────────────
    figure('Name', sprintf('Сходимость: %s', fname), ...
           'NumberTitle', 'off', 'Color', 'white');

    histories = {hist_sd, hist_cg, hist_pr};

    % Поверхность (только внутри сетки)
    x1lin = linspace(g(1), g(2), g(5));
    x2lin = linspace(g(3), g(4), g(5));
    [X1g, X2g] = meshgrid(x1lin, x2lin);
    Zg = zeros(size(X1g));
    for r = 1:size(X1g,1)
        for c = 1:size(X1g,2)
            Zg(r,c) = f_vec([X1g(r,c), X2g(r,c)]);
        end
    end

    surf(X1g, X2g, log1p(Zg), 'EdgeColor', 'none', 'FaceAlpha', 0.80);
    colormap(gca, jet);
    hold on;

    % Траектории — без клиппинга координат, чтобы видеть реальный путь
    for k = 1:numel(histories)
        H = histories{k};
        if isempty(H) || size(H,1) < 2; continue; end

        % z-координата по реальным значениям функции (не клипируем)
        Zh = arrayfun(@(r) log1p(f_vec(H(r,:))), 1:size(H,1));

        plot3(H(:,1), H(:,2), Zh, ...
            '-', ...
            'Color',     colors{k}, ...
            'LineWidth', 1.8, ...
            'DisplayName', method_names{k});

        % Маркеры только в начале и в конце траектории
        plot3(H(1,1),   H(1,2),   Zh(1),   markers{k}, ...
            'Color', colors{k}, 'MarkerFaceColor', colors{k}, ...
            'MarkerSize', 7, 'HandleVisibility', 'off');
        plot3(H(end,1), H(end,2), Zh(end), markers{k}, ...
            'Color', colors{k}, 'MarkerFaceColor', 'white', ...
            'MarkerSize', 7, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end

    % Начальная точка
    z_x0 = log1p(f_vec(x0));
    plot3(x0(1), x0(2), z_x0, ...
        'ks', 'MarkerSize', 11, 'LineWidth', 2, ...
        'MarkerFaceColor', [1 0.9 0], 'DisplayName', 'x^0');

    % Аналитический минимум — z вычисляется через f_opt
    z_opt = log1p(f_opt);
    plot3(x_opt(1), x_opt(2), z_opt, ...
        'k*', 'MarkerSize', 13, 'LineWidth', 2, 'DisplayName', 'x^*');

    xlabel('x_1',      'FontSize', 12);
    ylabel('x_2',      'FontSize', 12);
    zlabel('log(1+f)', 'FontSize', 12);
    title(sprintf('Траектории сходимости | %s | x0=[%.1f,%.1f] | E=%.0e', ...
                  fname, x0(1), x0(2), E), 'FontSize', 13);
    legend('Location', 'northeast', 'FontSize', 9);
    grid on;
    view(-40, 30);
    rotate3d on;

    % Оси — расширяем до фактических границ траекторий,
    % но не меньше границ сетки поверхности
    all_x1 = g(1:2);
    all_x2 = g(3:4);
    for k = 1:numel(histories)
        H = histories{k};
        if isempty(H); continue; end
        all_x1 = [all_x1, H(:,1)'];  %#ok<AGROW>
        all_x2 = [all_x2, H(:,2)'];  %#ok<AGROW>
    end
    all_x1 = [all_x1, x0(1)];
    all_x2 = [all_x2, x0(2)];

    pad = 0.15;
    xlim([min(all_x1)-pad, max(all_x1)+pad]);
    ylim([min(all_x2)-pad, max(all_x2)+pad]);
end

fprintf('\n');
fprintf('==========================================================================\n');
fprintf('  ВСЕ ТЕСТЫ ЗАВЕРШЕНЫ\n');
fprintf('==========================================================================\n');