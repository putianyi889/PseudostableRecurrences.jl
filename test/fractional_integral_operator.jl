using ClassicalOrthogonalPolynomials
using LinearAlgebra
using BandedMatrices
using SpecialFunctions

OpR(γ,δ,α,β) = Jacobi(γ,δ) \ Jacobi(α,β); # R_{(α,β)}^{(γ,δ)}
OpL(α,β,k,j) = Jacobi(α,β) \ (JacobiWeight(k,j).*Jacobi(α+k,β+j)); # L_{(α+k,β+j)}^{(α,β)}
OpW(β) = Diagonal(β+1:∞); # W^{(β)}
OpInvW(β) = Diagonal(inv.(β+1:∞)); # W^{(β)}^{-1}
OpΛ(r,μ,k) = BandedMatrix(-k=>(gamma.((r+1):μ:∞)./gamma.((r+k*μ+1):μ:∞)))

OpI(α,β,b,p) = p//2^(p-1)*OpL(α,β,0,p)*OpR(α,β+p,α-1,b+p)*OpInvW(b+p-1)*OpR(α,b+p-1,α,β); # I^{(α,β)}_{b,p}
function OpI(α,β,b,p,N) 
    L = p>1 ? *(map(x->x[1:N,1:N],OpL(α,β,0,p).args)...) : OpL(α,β,0,p)[1:N,1:N]
    R1 = b!=β ? *(map(x->x[1:N,1:N],OpR(α,β+p,α-1,b+p).args)...) : OpR(α,β+p,α-1,b+p)[1:N,1:N]
    R2 = b+p-1-β>1 ? *(map(x->x[1:N,1:N],OpR(α,b+p-1,α,β).args)...) : OpR(α,b+p-1,α,β)[1:N,1:N]
    p/2^(p-1)*L*R1*Diagonal(OpInvW(b+p-1).diag[1:N])*R2
end

fracpochhammer(a,b,n) = prod(x/y for (x,y) in zip(range(a,length=n),range(b,length=n)); init=one(promote_type(typeof(a),typeof(b))))
fracpochhammer(a,b,stepa,stepb,n) = prod(x/y for (x,y) in zip(range(a,step=stepa,length=n),range(b,step=stepb,length=n)); init=one(promote_type(typeof(a),typeof(b))))
OpCsingle(α,β,k,n) = (-1)^(n-k)*fracpochhammer(k+β+1,one(β),n-k)*fracpochhammer(n+α+β+1,2*one(β),1,2,k); # [k,n]
function OpCcolumn(α,β,n) # [:,n]
    (α,β)=promote(α,β)
    ret=zeros(typeof(α),n+1)
    ret[1]=(-1)^n*fracpochhammer(β+1,1,n);
    for k=1:n
        ret[k+1]=ret[k]*(n+α+β+k)*(n-k+1)/(-2*k*(k+β))
    end
    ret
end
function OpC(α,β,N)
    (α,β)=promote(α,β)
    ret=zeros(typeof(α),N+1,N+1)
    for n=0:N
        ret[1:n+1,n+1] = OpCcolumn(α,β,n)
    end
    ret
end

function OpI11(α,β,b,p,μ,N)
    k=Int(round(p*μ))
    (α,β,b,p,μ) = promote(α,β,b,p,μ)
    A = OpC(α,β,N+k)
    return A\(2^(μ-k)*OpΛ(b/p,1/p,k)[1:N+k+1,1:N+1]*A[1:N+1,1:N+1])
end