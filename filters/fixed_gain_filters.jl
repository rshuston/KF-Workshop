using Plots

t_max = 100
w = 2 * π / t_max

t = 0:1:t_max
N = length(t)

#x = 3 .* sin.(w .* t)
x = @. 3 * sin(w * t)
v = rand(N) .- 0.5
z = x .+ v

alpha = 0.2
beta = alpha^2 / (2 - alpha)

x_hat = zeros(N)
x_hat[1] = z[1]
for k in 2:N
    x_p = x_hat[k-1]
    x_hat[k] = x_p + alpha * (z[k] - x_p)
end
p_alpha = plot(t, [x z x_hat],
    seriestype = :scatter, markersize=2,
    title = "Values, alpha = $alpha")

x_hat = zeros(N)
v_hat = zeros(N)
x_hat[1] = z[1]
v_hat[1] = 0
for k in 2:N
    T = t[k] - t[k-1]
    x_p = x_hat[k-1] + v_hat[k-1] * T
    v_p = v_hat[k-1]
    delta = z[k] - x_p
    x_hat[k] = x_p + alpha * delta
    v_hat[k] = v_p + beta * delta / T
end
p_alpha_beta = plot(t, [x z x_hat],
    seriestype=:scatter, markersize=2,
    title="Values, alpha = $alpha, beta = $beta")

plot(p_alpha, p_alpha_beta,
    layout = (2,1),
    label = ["x" "z" "x_hat"])
