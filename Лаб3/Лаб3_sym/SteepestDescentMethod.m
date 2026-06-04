function [x_min, f_min, t, n, history] = SteepestDescentMethod(x0, E, f)
    tic;
    x   = x0(:);
    h   = 0.025;
    dim = length(x0);
    kmax = 1e4*dim;
    history = zeros(kmax + 1, dim);
    history(1, :) = x';
    hist_idx = 1;

    vars = symvar(f);
    grad_sym = gradient(f, vars);
    f_num = matlabFunction(f, 'Vars', {vars});
    grad_num = matlabFunction(grad_sym, 'Vars', {vars});
    
    grad_x0 = grad_num(x');
    y = -grad_x0;
    if norm(y) > 1e5
        y = y / norm(y);
    end
    func = @(alpha) f_num((x + alpha * y)');
    alpha = GoldenSection(0, h, E, func);
    x_prev = x;
    x = x + alpha * y;
    hist_idx = hist_idx + 1;
    history(hist_idx, :) = x';
    n = 1;
    if norm(x - x_prev) < E             % ИЗМЕНЕНО: было norm(grad_num(x')) < E
        x_min = x;
        f_min = f_num(x_min');
        t  = toc;
        history = history(1:hist_idx, :);
        return;
    end
    while n <= kmax
        grad_x = grad_num(x');
        y      = -grad_x;
        if norm(y) > 1e5
            y = y / norm(y);
        end
        func = @(alpha) f_num((x + alpha * y)');
        alpha = GoldenSection(0, h, E, func);
        x_prev = x;
        x = x + alpha * y;
        hist_idx = hist_idx + 1;
        history(hist_idx, :) = x';
        n = n + 1;
        if norm(x - x_prev) < E         % ИЗМЕНЕНО: было norm(grad_num(x')) < E
            break;
        end
    end
    x_min = x;
    f_min = f_num(x_min');
    t     = toc;
    history = history(1:hist_idx, :);
end