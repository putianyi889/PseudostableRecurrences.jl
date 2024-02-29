# PseudostableRecurrences.jl

## Introduction
Consider the recurrence
```math
x_1=x_2=x_3=\sqrt{2}, x_{n+3}-\frac{10}{3}x_{n+2}+3x_{n+1}-\frac{2}{3}x_n=0.
```
Its eigenvalues are ``2``, ``\frac{1}{3}`` and ``1``, so the solution is constant - ``x_n=\sqrt{2}``. However, the forward recurrence is numerically unstable since an eigenvalue is greater than 1. The following experiment confirms that.

```@example
using Plots
N = 100
x = zeros(N)
x[1] = sqrt(2)
x[2] = sqrt(2)
x[3] = sqrt(2)
for k in 1:N-3
    x[k+3] = 10/3*x[k+2] - 3*x[k+1] + 2/3*x[k]
end
plot(4:N, abs.(x.-sqrt(2))[4:N], label="numerical error", yaxis=:log, background_color=:transparent, foreground_color=:gray)
```

The backward recurrence, on the other hand, will be dominated by the eigenvalue ``\frac{1}{3}`` instead. While there may be other algorithms that are stable in this case, they often requires the asymptotic behavior of the sequence and vary across different problems.

The idea of pseudo-stablisation is to use a high precision to achieve a given error tolerance despite the instability. The precision choice is smart such that no extra care is taken by the end user, resulting in a "stable" behavior. Needless to say, pseudo-stablisation only works on problems where the initial values and the recurrence coefficients can be computed to arbitrary precision.

The current iteration of the strategy focuses on linear problems, where the end error scales linearly with the initial error which is considered as machine epsilon. To get the scaling factor, a test recurrence with random initial values is run. The random initial values ensure that the initial error is significant and the leading eigenvalue immediately takes over. The estimated scaling factor is then the ratio between the end norm of the test run and the initial norm.

The test recurrence may sound like a performance trade-off, but since it is run under a low precision, the cost is negligible.

For how `PseudostableRecurrence.jl` solves this recurrence and a few further examples, navigate to [Examples](@ref).