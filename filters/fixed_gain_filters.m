function fixed_gain_filters
    t_max = 100;
    w = 2 * pi / t_max;
    
    t = [0 : 1 : t_max];
    x = 3 * sin(w * t);
    v = (rand(1, length(t)) - 0.5);
    z = x + v;
    
    alpha = 0.2;
    beta = alpha * alpha / (2 - alpha);
    
    x_hat = zeros(1, length(t));
    
    x_hat(1) = z(1);
    for k = 2:length(t)
        x_p = x_hat(k-1);
        x_hat(k) = x_p + alpha * (z(k) - x_p);
    endfor
    
    figure("name", ["Values, alpha = " num2str(alpha)]);
    plot(t, x, ".", t, z, "o", t, x_hat, "*");
    set(gca(), "xgrid", "on", "ygrid", "on");
    
    x_hat = zeros(1, length(t));
    v_hat = zeros(1, length(t));
    
    x_hat(1) = z(1);
    v_hat(1) = 0;
    for k = 2:length(t)
        T = t(k) - t(k-1);
        x_p = x_hat(k-1) + v_hat(k-1) * T;
        v_p = v_hat(k-1);
        delta = z(k) - x_p;
        x_hat(k) = x_p + alpha * delta;
        v_hat(k) = v_p + beta * delta / T;
    endfor
    
    figure("name", ["Values, alpha = " num2str(alpha) ", beta = " num2str(beta)]);
    plot(t, x, ".", t, z, "o", t, x_hat, "*");
    set(gca(), "xgrid", "on", "ygrid", "on");
endfunction