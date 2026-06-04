function [x_min, f_min, t, n, history] = ConjugateGradientMethod(x0, E, f)
        tic;
        x   = x0(:);
        h   = 0.025;
        dim = length(x0);
        kmax = 1e4*dim;
        history = zeros(kmax + 1, dim);
        history(1, :) = x';
        hist_idx = 1;

        grad_x0 = gradient(f, x, E/10);
        s = -grad_x0;
        if norm(s) > 1e5
            s = s / norm(s);
        end
        func = @(alpha) f(x + alpha * s);
        alpha = GoldenSection(0, h, E, func);
        x_prev = x;
        x = x + alpha * s;
        hist_idx = hist_idx + 1;
        history(hist_idx, :) = x';
        grad_norm_prev = norm(grad_x0);
        s_prev  = s;
        n = 1;
        cycle_count = 0;
    if norm(x - x_prev) < E                  % ИЗМЕНЕНО: было norm(grad_num(x')) < E
            x_min = x;
            f_min = f(x_min);
            t  = toc;
            history = history(1:hist_idx, :);
    return;
    end
    while n <= kmax
    if cycle_count >= dim
                grad_x = gradient(f, x, E/10);
                s = -grad_x;
                if norm(s) > 1e5
                    s = s / norm(s);
                end
                func = @(alpha) f(x + alpha * s);
                alpha = GoldenSection(0, h, E, func);
                x_prev = x;
                x = x + alpha * s;
                hist_idx = hist_idx + 1;
                history(hist_idx, :) = x';
                grad_norm_prev = norm(grad_x);
                s_prev = s;
                n = n + 1;
                cycle_count = 0;
    if norm(x - x_prev) < E          % ИЗМЕНЕНО: было norm(grad_num(x')) < E
        break;
    end
        continue;
    end
            grad_x = gradient(f, x, E/10);
            grad_norm_curr = norm(grad_x);
            w = (grad_norm_curr^2) / (grad_norm_prev^2);
            s = -grad_x + w * s_prev;
            if norm(s) > 1e5
                s = s / norm(s);
            end
            func  = @(alpha) f(x + alpha * s);
            alpha = GoldenSection(0, h, E, func);
            x_prev = x;
            x = x + alpha * s;
            hist_idx = hist_idx + 1;
            history(hist_idx, :) = x';
            grad_norm_prev = grad_norm_curr;
            s_prev = s;
            n = n + 1;
            cycle_count = cycle_count + 1;

    if norm(x - x_prev) < E             % ИЗМЕНЕНО: было norm(grad_num(x')) < E
        break;
    end

    end
        x_min = x;
        f_min = f(x_min);
        t     = toc;
        history = history(1:hist_idx, :);
end