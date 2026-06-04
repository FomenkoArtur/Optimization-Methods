function [x_min, f_min, t, n, history] = SteepestDescentMethod(x0, E, f)
    tic;
    x   = x0(:);
    h   = 0.025;
    dim = length(x0);
    kmax = 1e4*dim;
    history = zeros(kmax + 1, dim);
    history(1, :) = x';
    hist_idx = 1;

    grad_x0 = gradient(f, x, E/10);
    y = -grad_x0;
    if norm(y) > 1e5
        y = y / norm(y);
    end
    func = @(alpha) f(x + alpha * y);
    alpha = GoldenSection(0, h, E, func);
    x_prev = x;
    x = x + alpha * y;
    hist_idx = hist_idx + 1;
    history(hist_idx, :) = x';
    n = 1;
if norm(x - x_prev) < E             % ИЗМЕНЕНО: было norm(grad_num(x')) < E
        x_min = x;
        f_min = f(x_min);
        t  = toc;
        history = history(1:hist_idx, :);
return;
end
while n <= kmax
        grad_x = gradient(f, x, E/10);
        y      = -grad_x;
        if norm(y) > 1e5
            y = y / norm(y);
        end
        func = @(alpha) f(x + alpha * y);
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
    f_min = f(x_min);
    t     = toc;
    history = history(1:hist_idx, :);
end