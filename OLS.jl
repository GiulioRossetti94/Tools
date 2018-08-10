function OLS(LHV,RHV,m=0)
    if size(RHV,1) != size(LHV,1)
        sLHV = size(RHV)
        sRHV = size(LHV)
        println("OLS: Dimension mismatch. RHV and LHV must have the same number of rows. Current rows are $sLHV and $sRHV")
    end
    
    T = size(LHV,1)
    N = size(LHV,2)
    K = size(RHV,2)
    bhat = RHV\LHV
    yhat = RHV * bhat
    
    errv = LHV - yhat
    MSE = mean(errv.^2)
    var_y = LHV - ones(T,1) * mean(LHV)
    var_y = mean(var_y.^2)
    
    g = RHV .* errv
    S0 = Newey_West_CM(g,m)
    Exx = - RHV'RHV/T
    V = inv(Exx)'S0 *inv(Exx) 
    sebv = sqrt.(diag(V/T))
    
    R2 = (1 - MSE./var_y)'
    R2adj = (1 - (MSE ./ var_y) *(T - 1)/(T - K))'
    return bhat,sebv,V,R2,R2adj
end
