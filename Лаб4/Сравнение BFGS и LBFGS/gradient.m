function g = gradient(func, x, h)
    x = x(:);
    n = length(x);
    g = zeros(n, 1);

    for i=1:n
        dx = zeros(n, 1);
        dx(i) = h;

        g(i) = (func(x + dx) - func(x - dx)) / (2 * h);
    end
end