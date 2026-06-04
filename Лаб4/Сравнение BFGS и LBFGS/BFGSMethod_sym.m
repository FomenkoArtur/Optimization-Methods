function [x_min, f_min, t, n, history] = BFGSMethod_sym(x0, E, f)
    tic;
    n = 0;
    x = x0(:);
    dim = length(x0);
    kmax = 1e4 * dim;
    h = 0.025;

    vars = symvar(f); 
    grad = gradient(f, vars);
    
    f_num = matlabFunction(f, 'Vars', {vars});
    grad_num = matlabFunction(grad, 'Vars', {vars});

    history = zeros(kmax + 1, dim);
    history(1, :) = x';
    hist_idx = 1;

    A = eye(dim);
    g = grad_num(x'); 
    g = g(:);
    
    while n <= kmax
        n = n + 1;
        
        p = -A * g;
        
        func = @(alpha) f_num( (x + alpha * p)' );
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
                
        denom = dx' * dg;
        
        if abs(denom) > 1e-10
            scalar_term = 1 + (dg' * A * dg) / denom;
            term1 = scalar_term * (dx * dx') / denom;

            term2 = (dx * (dg' * A) + (A * dg) * dx') / denom;
            
            delta_A = term1 - term2;
            
            A = A + delta_A;
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