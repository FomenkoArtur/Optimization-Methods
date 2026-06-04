function [x_min, f_min, t, n, history] = LBFGSMethod_sym(x0, E, f)

    tic;
    n = 0;
    x = x0(:);
    dim = length(x0);
    kmax = 1e4 * dim;
    h = 0.025;

    vars = symvar(f);
    grad = gradient(f, vars);
    f_num   = matlabFunction(f,    'Vars', {vars});
    grad_num = matlabFunction(grad, 'Vars', {vars});

    history = zeros(kmax + 1, dim);
    history(1, :) = x';
    hist_idx = 1;

    m = 10;

    S = zeros(dim, m);
    Y = zeros(dim, m);

    rho = zeros(1, m);

    mem_count = 0;


    mem_ptr = 0;

    g = grad_num(x');
    g = g(:);

    while n <= kmax
        n = n + 1;

        p = TwoLoopRecursion(g, S, Y, rho, mem_count, mem_ptr, m);

        func  = @(alpha) f_num( (x + alpha * p)');
        alpha = GoldenSection(0, h, E, func);

        x_old = x;
        g_old = g;

        x = x + alpha * p;

        hist_idx = hist_idx + 1;
        history(hist_idx, :) = x';

        g = grad_num(x');
        g = g(:);

        dx = x - x_old;
        dg = g - g_old;

        denom = dx' * dg; % должно быть > 0 (условие кривизны)

        if abs(denom) > 1e-10 && norm(dx) > 1e-12

            mem_ptr = mod(mem_ptr, m) + 1;
            S(:, mem_ptr) = dx;
            Y(:, mem_ptr) = dg;
            rho(mem_ptr)  = 1.0 / denom;

            if mem_count < m
                mem_count = mem_count + 1;
            end
        end

        grad_norm = norm(g);
        if grad_norm < E
            break;
        end
    end

    x_min = x;
    f_min = f_num(x_min');
    t = toc;
    history = history(1:hist_idx, :);
end


function p = TwoLoopRecursion(g, S, Y, rho, mem_count, mem_ptr, m)

    if mem_count == 0
        p = -g;
        return;
    end

    q      = g;
    alphas = zeros(1, mem_count);

    for j = 1 : mem_count
        idx = mod(mem_ptr - j + m, m) + 1;
        alphas(j) = rho(idx) * (S(:, idx)' * q);
        q = q - alphas(j) * Y(:, idx);
    end

    last_idx = mem_ptr;
    gamma = (S(:, last_idx)' * Y(:, last_idx)) / ...
            (Y(:, last_idx)' * Y(:, last_idx));
    r = gamma * q;


    for j = mem_count : -1 : 1
        idx  = mod(mem_ptr - j + m, m) + 1;
        beta = rho(idx) * (Y(:, idx)' * r);
        r    = r + (alphas(j) - beta) * S(:, idx);
    end

    p = -r;
end