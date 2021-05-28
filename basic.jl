
using Distributed
using Distances
using StatsBase 
using Profile    
using Random
using Distributions
rng = MersenneTwister(1234)
using DelimitedFiles
using Statistics
using LinearAlgebra

softmax(x; dims=1) = exp.(x) ./ sum(exp.(x), dims=dims)

linreg(x, y) = hcat(fill!(similar(x), 1), x) \ y

"""
	f1_f2(X::Array{Int64, 2}, W::Array{Float64, 2}, q::Int64)
	Compute frequencies of single and double sites, f1((i-1)*q+a), f2((i-1)*q+a, (j-1)*q+b).
	Output: f1::Array{Float64, 1} and f2::Array{Float64, 2}.
"""
function f1_f2_c2(X::Array{Int64, 2}, W::Array{Float64, 1}, q::Int64)
	M,L = size(X)
	Meff = sum(W); scale = 1.0 /Meff
	f1 = zeros(Float64, L*q)
	f2 = zeros(Float64, L*q, L*q)
	c2 = zeros(Float64, L*q, L*q)
	# Better to make it to be matrix vector operation.	
	for m in 1:M
		for i in 1:L
			a = X[m,i]+1
			f1[(i-1)*q+a] += W[m] * scale
			for j in (i+1):L
				b = X[m,j]+1
				f2[(i-1)*q+a, (j-1)*q+b] += W[m] * scale
				f2[(j-1)*q+b, (i-1)*q+a] += W[m] * scale
			end
		end
	end
	#c2 = f2 - f1*f1'
		
	for i in 1:L
		for j in (i+1):L
			for a in 1:q
				for b in 1:q
					c2[(i-1)*q+a, (j-1)*q+b] = f2[(i-1)*q+a, (j-1)*q+b] - f1[(i-1)*q+a] * f1[(j-1)*q+b] 
					c2[(j-1)*q+b, (i-1)*q+a] = f2[(j-1)*q+b, (i-1)*q+a] - f1[(i-1)*q+a] * f1[(j-1)*q+b]
				end
			end
		end
	end
	return (Meff, L, f1, f2, c2) 
end

function fitquality(f2_1::Array{Float64,2}, f1_1::Array{Float64,1}, f2_2::Array{Float64,2}, f1_2::Array{Float64,1}, q::Int64; withdiag=false)
    L = Int64(size(f2_1,1)/q)
    C1 = f2_1 - f1_1*f1_1'
    C2 = f2_2 - f1_2*f1_2'

    c1vec = zeros(Float64, Int64(L*(L-1)*q*q/2) )
    c2vec = zeros(Float64, Int64(L*(L-1)*q*q/2) )
    pos = 0
    for i in 1:L
        # id[(i-1)*q .+ (1:q), (i-1)*q .+ (1:q)] .= withdiag
        for j in (i+1):L
            # id[(j-1)*q .+ (1:q),(i-1)*q .+ (1:q)] .= true
            c1vec[pos .+ (1:q^2)] .= vec(C1[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)])
            c2vec[pos .+ (1:q^2)] .= vec(C2[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)])
            pos += q^2
        end
    end

    cc = Statistics.cor(c1vec,c2vec)
    cslope = linreg(c1vec,c2vec)[2]
    froc = LinearAlgebra.norm(c1vec - c2vec)
    cm = Statistics.cor(f1_1,f1_2)
    from = LinearAlgebra.norm(f1_1-f1_2)
    return (cc,cslope,froc,cm,from)
end

function get_independent_site_model_profile(fname,L,rng)
    X=readdlm(fname, Int)
    (M,L2) = size(X)
    for i in 1:L
        shuffle_vec = randperm(rng,M)
        X[:,i] = copy(X[shuffle_vec,i])
    end
    return X
end

function get_independent_site_model_profile(fname::String, L::Int64, rng::MersenneTwister,th::Float64)
    X=readdlm(fname, Int)
    X_temp=[]
    n_id=1

    for n in 1:size(X,1) 
        r=rand()
        if(r<th)
            if(n_id==1)
                X_temp = X[n,:]'
            end
            if(n_id>1)
                X_temp = vcat(X_temp, X[n,:]')
            end
            n_id+=1
        end
    end
    M=size(X_temp,1)
	
    Xcopy = copy(X_temp)
    (M,L2) = size(Xcopy)
    for i in 1:L
        shuffle_vec = randperm(rng,M)
        Xcopy[:,i] = copy(Xcopy[shuffle_vec,i])
    end
    
    Xcopy = Xcopy+ ones(Int, size(Xcopy))
    Xout = zeros(Int, (M, L*q))

    for m in 1:M
        for i in 1:L
            Xout[m, ((i-1)*q+Xcopy[m,i])] = 1
        end
    end
    Xout=Xout'
    return Xout
end

function binary2categorical(n, x)
    q = Int(length(x)//n)
    vec_temp = [argmax( x[((i-1)*q+1):i*q] ) for i in 1:n]
    return vec_temp
end

function categorical2binary(X, q)
	(N, L) = size(X)
	X_binary = zeros(Int64, (N, q*L))
	for n in 1:N
	    for i in 1:L
		X_binary[n, km(i,X[n,i]+1,q)] = 1
	    end
	end
	return X_binary'
end

function write_seq_to_num_file(fname, Data)
    fout_seq = open(fname, "w")
    N,L = size(Data)
    for n in 1:N
        for i in 1:L
            print(fout_seq, Data[n,i], " ")
        end
        println(fout_seq,"")
    end
    close(fout_seq)
end
