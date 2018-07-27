
function genFock(rho,h,V)
    F = h
    n = size(h)[1] #read size of basis set
    for i = 1:n
        for j = 1:n
            for k = 1:n
                for l = 1:n
                    F[i,j] += (2*V[i,l,j,k]-V[i,l,k,j])*rho[k,l]
                end
            end
        end
    end
    return F
end

#= ## Need to figure out how to implement lists of matrices in Julia
function diis(errs, Fs)
    n = size(errs)[1]
    B = fill(0., (n,n))
    for i = 1:n
        for j = 1:n
            B[i,j] =  #XXX what does raveling the errors do?
        end
    end
    A[1:n,1:n] = B
    A[n,:] = -1
    A[:,n] = -1
    A[n,n] = 0
    b = fill(0., n+1)
    b[n] = -1
    x = A\b

    w = [1:n]
    F = fill(0., (n,n))
    for i = 1:n
        F += w[i] * Fs[i]
    end
    return F
end
=#

function do_hf(m,h,V)
    n = size(h)[1]
    maxiter = 200
    converged = false
    rho = Diagonal(ones(n))

    F = genFock(rho, h, V)
    println("rho")
    println(rho)
    println("F")
    println(F)
    println()
    for iter = 1:maxiter
        w,C = eig(F)
        Cocc = C[:,1:m]
        println(Cocc)
        rho_new = Cocc*Cocc'
        F_new = genFock(rho_new,h,V)

        err = F_new*rho_new - rho_new*F_new
        if norm(err) < 0.1^8
            converged = true
            println(iter)
            println(norm(err))
            break
        end

        rho = rho_new
        F = F_new
    end

    return (rho, F)

end

n = 10
m = 5
U = 4
h = fill(0., (n,n))
V = fill(0., (n,n,n,n))

#Create h and V for hubbard model
for i = 1:n-1
    h[i,i+1] = -1
end
h[1,n] = -1
h = h + h'

for i = 1:n
    V[i,i,i,i] = U
end

println(h)

rho, F = do_hf(m,h,V)
println("rho")
println(rho)
println("F")
println(F)
