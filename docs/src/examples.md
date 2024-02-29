# Examples

## Linear-recursive Sequence
Consider the recurrence
```math
x_1=x_2=x_3=\sqrt{2}, x_{n+3}-\frac{10}{3}x_{n+2}+3x_{n+1}-\frac{2}{3}x_n=0.
```

First, we shall recognise that the linear recursive sequence is basically a 1D stencil recurrence. The recurrence needs to be converted to an assignment first, that is, ``x_{n}=\frac{10}{3}x_{n-1}-3x_{n-2}+\frac{2}{3}x_{n-3}``. Such a recurrence can be defined by the stencil with the coefficients
```@setup 1
using PseudostableRecurrences
```
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
P = StencilRecurrencePlan{Real}(stencil, coefs, f_init, (100,)) # 'Real' specifies the domain of the entries, as opposed to 'Complex', etc.
stable_recurrence(P) # defaults to stable_recurrence(P, Float64)
```

## Special integrals
Consider the integral
```math
I(m,n) = \int_0^\pi\left(\frac{\cos x}{1+\sin x}\right)^me^{inx}\mathrm{d}x
```
which has the recurrence
```math
I(m,n)=\begin{cases}
    2i/n,& m=0, n\text{ odd}\\
    \displaystyle(-1)^{m/2}\left(\pi-4\sum_{k=1}^{m/2}\frac{(-1)^k}{2k-1}\right),& m\text{ even}, n=0\\
    0,& \!\!\!\begin{array}{l}
        m=0\text{ or }n=0,\\
        \text{excluding above cases}
    \end{array}\\
    \displaystyle\frac{m+n-1}{ni}I(m,n-1)+\frac{m}{n}I(m-1,n-1)+\frac{2}{n}\begin{cases}
        i,\\
        -1,
    \end{cases}&\!\!\!\begin{array}{l}
        m+n\text{ odd}\\
        m+n\text{ even}.
    \end{array}
\end{cases}
```
Our goal is to obtain the numerical integrals for a range of ``m`` and ``n``. We may use a matrix to store the results, with `m+1` and `n+1` being the row and column indices.

The problem is again a stencil recurrence, where
```@setup 2
using PseudostableRecurrences, BenchmarkTools, QuadGK
```
```@example 2
stencil = (CartesianIndex(-1,-1), CartesianIndex(0,-1), CartesianIndex(0,0))
```
and
```@example 2
coef1(m,n) = (m-1)//(n-1)
coef2(m,n) = (m+n-3)//(n-1)//im
coef3(m,n) = 2//(n-1)*ifelse(isodd(m+n), im, -1)
coef = (coef1, coef2, coef3)
```
Note that Julia indices starts at 1, hence the `m-1`, `n-1`, etc.. Recall that the `CartesianIndex(0,0)` and `coef3` corresponds to the inhomogeneous term
```math
\frac{2}{n}\begin{cases} i,\\ -1, \end{cases}&\!\!\!\begin{array}{l} m+n\text{ odd}\\ m+n\text{ even}. \end{array}
```

We need to initialise the first column as well.
```@example 2
function f_init(T, m)
    A = zeros(Complex{T},m,2)
    A[1,1] = π
    for mm in 3:2:m
        A[mm,1] = -A[mm-2,1] + 4//(mm-2)
    end
    A
end
```
Here, `m` specifies the size of the first dimension. For an ``n``-dimensional stencil recurrence, each slice has ``n-1`` dimensions, so the initialization function has to take ``n-1`` size arguments.

!!! note "Initialising the first row"
    `PseudostableRecurrences.jl` doesn't support initialising the first row in the backend. However, if a stencil exceeds the boundary, the external entries will be ignored. This allows one to modify the coefficients to achieve the initialization of the first row. In this example, the recurrence happens to initialise the first row correctly, so there is no need to modify it.

Then we are ready to run the recurrence.
```@example 2
P = StencilRecurrencePlan{Complex}(stencil, coef, f_init, (101,101), (1,2))
ret = hcat(stable_recurrence(P)...)
```
By design, `stable_recurrence` returns a vector of slices, so `hcat` is used to re-assemble the result. The last argument of `StencilRecurrencePlan` specifies the index where the recurrence starts at. By default, it will start from where the stencil doesn't go out of bound. In this example, however, we do want to apply the stencil to the first row in order to initialise it.

Now we compare the results against that from [`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl) to verify the correctness:
```@example 2
maximum(abs, ret[92:end, 92:end] - [quadgk(x->(cos(x)/(1+sin(x)))^m*exp(im*n*x), 0, π)[1] for m in 91:100, n in 91:100])
```

!!! note "Performance"
    The recurrence is faster than generic numerical integration when the intermediate results are useful. In this example, when `m` is small and `n` is large, the function is oscillatory which `QuadGK.jl` has difficulty to handle:
    ```@repl 2
    @btime quadgk(x->(cos(x)/(1+sin(x)))^100*exp(im*100*x), 0, π);
    @btime quadgk(x->(cos(x)/(1+sin(x)))^2*exp(im*100*x), 0, π);
    @btime quadgk(x->exp(im*100*x), 0, π);
    ```
    This is why we don't compare the whole matrix.