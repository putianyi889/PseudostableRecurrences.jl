# Examples

## Linear-recursive Sequence

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

## Pseudo-stablisation methods
The idea of pseudo-stablisation is to use a high precision to achieve a given error tolerance despite the instability. The precision choice is smart such that no extra care is taken by the end user, resulting in a "stable" behavior. Needless to say, pseudo-stablisation only works on problems where the initial values and the recurrence coefficients can be computed to arbitrary precision.

The current iteration of the strategy focuses on linear problems, where the end error scales linearly with the initial error which is considered as machine epsilon. To get the scaling factor, a test recurrence with random initial values is run. The random initial values ensure that the initial error is significant and the leading eigenvalue immediately takes over. The estimated scaling factor is then the ratio between the end norm of the test run and the initial norm.

The test recurrence may sound like a performance trade-off, but since it is run under a low precision, the cost is negligible.

### Linear-recursive Sequence
First, we shall recognise that the linear recursive sequence is basically a 1D stencil recurrence. The recurrence needs to be converted to an assignment first, that is, ``x_{n}=\frac{10}{3}x_{n-1}-3x_{n-2}+\frac{2}{3}x_{n-3}``. Such a recurrence can be defined by the stencil with the coefficients
```@example 1
stencil = (CartesianIndex(-3), CartesianIndex(-2), CartesianIndex(-1))
coefs = (n -> 2//3, n -> -3, n -> 10//3)
```

!!! note "Inhomogeneous case"
    To define a linear inhomogeneous recurrence, the coefficient associated with the zero cartesian index is used. For example, the recurrence
    ```math
    x_{n}=\frac{10}{3}x_{n-1}-3x_{n-2}+\frac{2}{3}x_{n-3}+\frac{1}{n}
    ```
    should be defined by
    ```julia
    stencil = (CartesianIndex(-3), CartesianIndex(-2), CartesianIndex(-1), CartesianIndex(0))
    coefs = (n -> 2//3, n -> -3, n -> 10//3, n -> 1//n)
    ```

We then need to define the initial values
```@example 1
f_init(T) = [sqrt(T(2)), sqrt(T(2)), sqrt(T(2)), T(0)]
```
where the place for the next step should be reserved. This is a function that can generate values based on type `T`.

Now we only need to further provide the size of the recurrence and we are ready to go.
```@example 1
using PseudostableRecurrences
P = StencilRecurrencePlan{Real}(stencil, coefs, f_init, (100,)) # 'Real' specifies the domain of the entries, as opposed to 'Complex', etc.
stable_recurrence(P) # defaults to stable_recurrence(P, Float64)
```