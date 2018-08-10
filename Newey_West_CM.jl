function Newey_West_CM(g,m=0)
    T = size(g,1)
    m = min(m,T-1) #lags
    
    g1 = g .- mean(g,1)
    
    s = g1'g1/T
    for lags = 1:m
        Gamma = g1[lags+1:T,:]'g1[1:T-lags,:]/T
        s = s + (1 - lags/(m+1))*(Gamma + Gamma')
    end
    
    return s
end
