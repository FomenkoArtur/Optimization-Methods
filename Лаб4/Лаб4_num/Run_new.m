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

grids = {
    [-2.5,  2.5,  -1.5,  3.5,  300];   % Rosenbrock
    [-5.0,  5.0,  -5.0,  5.0,  250];   % Ackley
    [-5.0,  5.0,  -5.0,  5.0,  250];   % Himmelblau
    [-5.0,  5.0,  -5.0,  5.0,  250];   % Booth
    [-5.5,  5.5,  -5.5,  5.5,  250];   % Rastrigin
};

tests = {                               % УБРАНО: столбец f_sym из каждой строки
    f_vec_roz,  [-2.0,  3.0], 1e-4, ...
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

% [ИЗМЕНЕНО] Добавлен 'LBFGS' в конец списка методов
method_names = {'Newton', 'NewtonRaphson', 'Secant', 'Broyden', 'DFP', 'BFGS', 'LBFGS'};

% [ИЗМЕНЕНО] Добавлен цвет для LBFGS, исправлена запятая у BFGS
colors = {
    [0.85, 0.10, 0.10];  % Newton        — красный
    [0.05, 0.65, 0.10];  % NewtonRaphson — зелёный
    [0.10, 0.30, 0.95];  % Secant        — синий
    [0.90, 0.50, 0.10];  % Broyden       — оранжевый
    [0.60, 0.10, 0.80];  % DFP           — фиолетовый
    [0.10, 0.70, 0.70];  % BFGS          — циан
    [0.80, 0.70, 0.05]   % LBFGS         — жёлто-золотой
};

% [ИЗМЕНЕНО] Добавлен маркер '<' для LBFGS
markers = {'o', 's', '^', 'd', 'v', '>', '<'};

fprintf('\n');
fprintf('==========================================================================\n');
fprintf('     СРАВНЕНИЕ МЕТОДОВ МИНИМИЗАЦИИ ФНП (методы второго порядка)\n');
fprintf('==========================================================================\n');

for i = 1:size(tests, 1)
    f_vec  = tests{i, 1};               % УБРАНО: f_sym = tests{i,2}; сдвинуты индексы столбцов
    x0     = tests{i, 2};
    E      = tests{i, 3};
    fname  = tests{i, 4};
    fdesc  = tests{i, 5};
    x_opt  = tests{i, 6};
    f_opt  = tests{i, 7};
    g      = grids{i};

    [x_newton,  f_newton,  t_newton,  n_newton,  hist_newton]  = NewtonMethod(x0, E, f_vec);      % УБРАНО: f_sym
    [x_nr,      f_nr,      t_nr,      n_nr,      hist_nr]      = NewtonRaphsonMethod(x0, E, f_vec); % УБРАНО: f_sym
    [x_secant,  f_secant,  t_secant,  n_secant,  hist_secant]  = SecantMethod(x0, E, f_vec);      % УБРАНО: f_sym
    [x_broyden, f_broyden, t_broyden, n_broyden, hist_broyden] = BroydenMethod(x0, E, f_vec);     % УБРАНО: f_sym
    [x_dfp,     f_dfp,     t_dfp,     n_dfp,     hist_dfp]     = DFPMethod(x0, E, f_vec);         % УБРАНО: f_sym
    [x_bfgs,    f_bfgs,    t_bfgs,    n_bfgs,    hist_bfgs]    = BFGSMethod(x0, E, f_vec);        % УБРАНО: f_sym
    % [ИЗМЕНЕНО] Добавлен вызов LBFGSMethod
    [x_lbfgs,   f_lbfgs,   t_lbfgs,   n_lbfgs,   hist_lbfgs]  = LBFGSMethod(x0, E, f_vec);

    fprintf('\n');
    fprintf('--------------------------------------------------------------------------\n');
    fprintf('  ФУНКЦИЯ : %s\n', fname);
    fprintf('  ФОРМУЛА : %s\n', fdesc);
    fprintf('  Начальная точка : x0=[%.2f, %.2f] | Точность: %.0e\n', x0(1), x0(2), E);
    fprintf('--------------------------------------------------------------------------\n');
    fprintf('  %-22s | x1=%-12.6f x2=%-12.6f | F=%-14.6e | T=%-11.7f | N=%-6d\n', ...
        'Newton',            x_newton(1),  x_newton(2),  f_newton,  t_newton,  n_newton);
    fprintf('  %-22s | x1=%-12.6f x2=%-12.6f | F=%-14.6e | T=%-11.7f | N=%-6d\n', ...
        'NewtonRaphson',     x_nr(1),      x_nr(2),      f_nr,      t_nr,      n_nr);
    fprintf('  %-22s | x1=%-12.6f x2=%-12.6f | F=%-14.6e | T=%-11.7f | N=%-6d\n', ...
        'Secant',            x_secant(1),  x_secant(2),  f_secant,  t_secant,  n_secant);
    fprintf('  %-22s | x1=%-12.6f x2=%-12.6f | F=%-14.6e | T=%-11.7f | N=%-6d\n', ...
        'Broyden',           x_broyden(1), x_broyden(2), f_broyden, t_broyden, n_broyden);
    fprintf('  %-22s | x1=%-12.6f x2=%-12.6f | F=%-14.6e | T=%-11.7f | N=%-6d\n', ...
        'DFP',               x_dfp(1),     x_dfp(2),     f_dfp,     t_dfp,     n_dfp);
    fprintf('  %-22s | x1=%-12.6f x2=%-12.6f | F=%-14.6e | T=%-11.7f | N=%-6d\n', ...
        'BFGS',              x_bfgs(1),    x_bfgs(2),    f_bfgs,    t_bfgs,    n_bfgs);
    % [ИЗМЕНЕНО] Добавлена строка вывода для LBFGS
    fprintf('  %-22s | x1=%-12.6f x2=%-12.6f | F=%-14.6e | T=%-11.7f | N=%-6d\n', ...
        'LBFGS',             x_lbfgs(1),   x_lbfgs(2),   f_lbfgs,   t_lbfgs,   n_lbfgs);
    fprintf('--------------------------------------------------------------------------\n');

    % MODIFIED SECTION: Print specific info for Himmelblau as per lecture
    if strcmp(fname, 'Himmelblau')
        fprintf('  Аналитические минимумы (4 точки):\n');
        fprintf('    x1* = [3, 2]\n');
        fprintf('    x2* = [-2.805, 3.131]\n');
        fprintf('    x3* = [-3.779, -3.283]\n');
        fprintf('    x4* = [3.584, -1.848]\n');
        fprintf('    f* = 0\n');
    else
        fprintf('  Аналитический минимум: x*=[%.4f, %.4f], f*=%.4f\n', x_opt(1), x_opt(2), f_opt);
    end

    figure('Name', sprintf('Сходимость: %s', fname), ...
           'NumberTitle', 'off', 'Color', 'white', 'Position', [100, 100, 1200, 800]);

    % [ИЗМЕНЕНО] Коллекция историй расширена до 7 методов
    histories = {hist_newton, hist_nr, hist_secant, hist_broyden, hist_dfp, hist_bfgs, hist_lbfgs};

    % [ИЗМЕНЕНО] Обрезка траекторий расширена до 7 методов
    hist_clipped = cell(7, 1);
    for k = 1:7
        H = histories{k};
        if isempty(H) || size(H, 1) < 1
            hist_clipped{k} = [];
            continue;
        end
        H(:, 1) = max(g(1), min(g(2), H(:, 1)));
        H(:, 2) = max(g(3), min(g(4), H(:, 2)));
        hist_clipped{k} = H;
    end

    % Построение сетки функции для поверхности
    x1lin = linspace(g(1), g(2), g(5));
    x2lin = linspace(g(3), g(4), g(5));
    [X1g, X2g] = meshgrid(x1lin, x2lin);
    Zg = zeros(size(X1g));
    for r = 1:size(X1g, 1)
        for c = 1:size(X1g, 2)
            Zg(r, c) = f_vec([X1g(r, c), X2g(r, c)]);
        end
    end

    % Поверхность с логарифмическим масштабом (log1p как в оригинале)
    surf(X1g, X2g, log1p(Zg), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    colormap(gca, jet);
    hold on;

    % Отрисовка траекторий всех 7 методов (numel подхватывает автоматически)
    for k = 1:numel(histories)
        H = hist_clipped{k};
        if isempty(H) || size(H, 1) < 2; continue; end

        % Вычисление log1p(f) для точек траектории
        Zh = arrayfun(@(r) log1p(f_vec(H(r, :))), 1:size(H, 1));

        plot3(H(:, 1), H(:, 2), Zh, ...
            '-o', ...
            'Color',           colors{k}, ...
            'LineWidth',       1.6, ...
            'MarkerSize',      3, ...
            'MarkerFaceColor', colors{k}, ...
            'Marker',          markers{k}, ...
            'DisplayName',     method_names{k});
    end

    % Стартовая точка (чёрный квадрат)
    z_x0 = log1p(f_vec(x0));
    plot3(x0(1), x0(2), z_x0, ...
        'ks', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'x^0');

    % Аналитический минимум (чёрная звезда)
    plot3(x_opt(1), x_opt(2), 0, ...
        'k*', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'x^*');

    % Подписи осей и заголовок
    xlabel('x_1',      'FontSize', 12);
    ylabel('x_2',      'FontSize', 12);
    zlabel('log(1+f)', 'FontSize', 12);
    title(sprintf('Траектории сходимости | %s | x0=[%.1f,%.1f] | E=%.0e', ...
                  fname, x0(1), x0(2), E), ...
          'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 9, 'Box', 'on');
    grid on;
    view(-40, 30);        % Тот же угол обзора, что в оригинале
    rotate3d on;

    % Фиксированные границы осей (как в оригинале)
    xlim([g(1), g(2)]);
    ylim([g(3), g(4)]);

    zlim([0, max(log1p(Zg(:))) * 1.05]);
end

fprintf('\n');
fprintf('==========================================================================\n');
% [ИЗМЕНЕНО] Обновлён счётчик: 7 методов
fprintf('  ВСЕ ТЕСТЫ ЗАВЕРШЕНЫ | Протестировано 7 методов на 5 функциях\n');
fprintf('==========================================================================\n');