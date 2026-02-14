function fixed_gain_filters
    t_max = 100;
    w = 2 * pi / t_max;

    t = [0 : 1 : t_max];
    N = length(t);

    x = 3 * sin(w * t);
    v = (rand(1, N) - 0.5);
    z = x + v;

    alpha = 0.2;
    beta = alpha * alpha / (2 - alpha);

    x_hat = zeros(1, N);
    x_hat(1) = z(1);
    for k = 2:N
        x_p = x_hat(k-1);
        x_hat(k) = x_p + alpha * (z(k) - x_p);
    end
    figure('name', ['Values, alpha = ', num2str(alpha)]);
    plot(t, x, '.', t, z, 'o', t, x_hat, '*')
    grid on

    x_hat = zeros(1, N);
    v_hat = zeros(1, N);
    x_hat(1) = z(1);
    v_hat(1) = 0;
    for k = 2:N
        T = t(k) - t(k-1);
        x_p = x_hat(k-1) + v_hat(k-1) * T;
        v_p = v_hat(k-1);
        delta = z(k) - x_p;
        x_hat(k) = x_p + alpha * delta;
        v_hat(k) = v_p + beta * delta / T;
    end
    figure('name', ['Values, alpha = ', num2str(alpha), ', beta = ', num2str(beta)]);
    plot(t, x, '.', t, z, 'o', t, x_hat, '*')
    grid on
end