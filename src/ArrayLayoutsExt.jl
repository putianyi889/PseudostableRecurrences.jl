import ArrayLayouts: rowsupport, colsupport

function slicesupport(A::AbstractMatrix, k; dims)
    if dims==1 || dims==(1,)
        (rowsupport(A, k), )
    elseif dims==2 || dims==(2,)
        (colsupport(A, k), )
    else
        throw(ArgumentError("slicesupport(::$(typeof(A)), $k; dims=$dims) is not supported."))
    end
end

function slicesupport(A::AbstractVector, k; dims)
    if dims==1
        ()
    else
        throw(ArgumentError("slicesupport(::$(typeof(A)), $k; dims=$dims) is not supported."))
    end
end