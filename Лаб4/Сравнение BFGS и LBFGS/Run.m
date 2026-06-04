clear; clc; format long
dim = 50;
E   = 1e-4;

% Построение символьного выражения Розенброка для dim=50
syms z [dim 1]
f_sym = sym(0);
for i = 1:2:dim-1
    f_sym = f_sym + 100*(z(i+1) - z(i)^2)^2 + (1 - z(i))^2;
end

x0_1 = -ones(dim, 1);
x0_2 =  2 * ones(dim, 1);
x0_3 = randn(dim, 1);

starts      = {x0_1, x0_2, x0_3};
start_names = {
    'x0 = -1', ...
    'x0 =  2', ...
    sprintf('x0 = rand [%.3f, %.3f, %.3f, ...]', x0_3(1), x0_3(2), x0_3(3))
};

fprintf('\n');
fprintf('==========================================================================\n');
fprintf('   ТЕСТ НА ФУНКЦИИ РОЗЕНБРОКА | dim=%d | E=%.0e\n', dim, E);
fprintf('==========================================================================\n');
fprintf('  Случайная начальная точка x0_3 (все координаты):\n  [');
fprintf(' %.4f', x0_3);
fprintf(' ]\n');

for s = 1:3
    x0 = starts{s};
    [x_bfgs_num,  f_bfgs_num,  t_bfgs_num,  n_bfgs_num,  ~] = BFGSMethod_num(x0,  E, @Rosenbrock);
    [x_bfgs_sym,  f_bfgs_sym,  t_bfgs_sym,  n_bfgs_sym,  ~] = BFGSMethod_sym(x0,  E, f_sym);
    [x_lbfgs_num, f_lbfgs_num, t_lbfgs_num, n_lbfgs_num, ~] = LBFGSMethod_num(x0, E, @Rosenbrock);
    [x_lbfgs_sym, f_lbfgs_sym, t_lbfgs_sym, n_lbfgs_sym, ~] = LBFGSMethod_sym(x0, E, f_sym);

    err_bfgs_num  = norm(x_bfgs_num  - ones(dim, 1));
    err_bfgs_sym  = norm(x_bfgs_sym  - ones(dim, 1));
    err_lbfgs_num = norm(x_lbfgs_num - ones(dim, 1));
    err_lbfgs_sym = norm(x_lbfgs_sym - ones(dim, 1));

    fprintf('\n');
    fprintf('--------------------------------------------------------------------------\n');
    fprintf('  НАЧАЛЬНАЯ ТОЧКА : %s\n', start_names{s});
    fprintf('--------------------------------------------------------------------------\n');
    fprintf('  %-16s | F=%-14.6e | ||x-x*||=%-12.6e | T=%-11.7f | N=%-6d\n', ...
        'BFGS_num',  f_bfgs_num,  err_bfgs_num,  t_bfgs_num,  n_bfgs_num);
    fprintf('  %-16s | F=%-14.6e | ||x-x*||=%-12.6e | T=%-11.7f | N=%-6d\n', ...
        'BFGS_sym',  f_bfgs_sym,  err_bfgs_sym,  t_bfgs_sym,  n_bfgs_sym);
    fprintf('  %-16s | F=%-14.6e | ||x-x*||=%-12.6e | T=%-11.7f | N=%-6d\n', ...
        'LBFGS_num', f_lbfgs_num, err_lbfgs_num, t_lbfgs_num, n_lbfgs_num);
    fprintf('  %-16s | F=%-14.6e | ||x-x*||=%-12.6e | T=%-11.7f | N=%-6d\n', ...
        'LBFGS_sym', f_lbfgs_sym, err_lbfgs_sym, t_lbfgs_sym, n_lbfgs_sym);
    fprintf('--------------------------------------------------------------------------\n');
    fprintf('  Аналитический минимум: x*=[1,...,1], f*=0\n');
end

fprintf('\n');
fprintf('==========================================================================\n');
fprintf('  ВСЕ ТЕСТЫ ЗАВЕРШЕНЫ | 4 метода x 3 начальные точки | dim=%d\n', dim);
fprintf('==========================================================================\n');