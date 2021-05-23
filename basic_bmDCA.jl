
"""
	init_h_J( f1::Array{Float64,1}, q::Int64, L::Int64, Meff::Float64; withdiag=false)
	return (h,J)	
"""
function init_h_J( f1::Array{Float64,1}, q::Int64, L::Int64, Meff::Float64; withdiag=false)
	#J = rand(Float64, (q*L, q*L)) * (1.0/q^2)
	J = rand(Float64, q*L, q*L) * (1.0/q^2)
	J = J+J'	
	#h = log.(f1+0.1/Meff * ones(q*L))
	h = randn(q*L)#log.(f1+0.1/Meff * ones(q*L))
	return (h,J) 
end

"""
	E_i(i::Int64, a::Int64, A::Array{Int64, 1}, J::Array{Float64, 2}, h::Array{Float64, 1})
	return e_i 
"""

function E_i(q::Int64, L::Int, i::Int64, A::Array{Int64, 1}, J::Array{Float64, 2}, h::Array{Float64, 1})
    	a = A[i]+1
 	e_i = - sum( J[km.(i,a,q), km.(1:L, A+ones(Int, L), q)] )	
	e_i += -h[ km(i,a,q) ]
	return e_i 
end 
"""
	dE_i(i::Int64, a1::Int64, a2::Int64, A::Array{Int64, 1}, J::Array{Float64, 2}, h::Array{Float64, 1})
	return de_i # = E_i(a_proposed) - E_i(a 
"""
function dE_i(q::Int64, i::Int64, a1::Int64, a2::Int64, A::Array{Int64, 1}, J::Array{Float64, 2}, h::Array{Float64, 1})
	de_i = 0.0
	for j=1:L
		if j!=i
			b = A[j]+1	
			de_i += -(J[km(i,a2, q), km(j,b, q)] - J[ km(i,a1, q), km(j,b, q) ])
		end
	end
	de_i += -(h[ km(i,a2,q)] - h[ km(i,a1,q)])
	return de_i 
end

"""
	Metropolis_Hastings(i::Int64, A::Array{Int64, 1}, J::Array{Float64, 2}, h::Array{Float64, 1})
	return (accepted, A) 
"""

function Metropolis_Hastings(E_old::Float64, q::Int64, L::Int64, i::Int64, A::Array{Int64, 1}, J::Array{Float64, 2}, h::Array{Float64, 1})
	a = A[i]
    a_proposed = rand(0:(q-1)) # Note entries of X and A are defined as between 0 and 20. 	
    if(a_proposed==a)
        a_proposed = rand(0:(q-1))
    end

    A_prop = copy(A); A_prop[i] = a_proposed
    E_new = E_i(q, L, i, A_prop, J, h)
    dE = E_new - E_old
	w = exp(-dE); 
	r = rand()    
    
    if(w>r)
        return (1, A_prop, E_new) 
    end
    if(w<=r)
        return(0, A, E_old)
    end
end
"""
	Monte_Carlo_sweep(L::Int64, A::Array{Int64, 1}, J::Array{Float64, 2}, h::Array{Float64, 1})
	return n_accepted
"""

function Monte_Carlo_sweep(E_old::Float64, q::Int64, L::Int64, A::Array{Int64,1}, J::Array{Float64,2}, h::Array{Float64,1})
	n_accepted = 0
	for l = 1:L
        i = rand(1:L)
		(accepted, A, E_old) = Metropolis_Hastings(E_old, q, L, i, A, J, h)
		n_accepted += accepted
	end
	return (n_accepted, A, E_old) 
end

function pCDk(X::Array{Int64, 2}, k_max::Int64, M::Int64, q::Int64, L::Int64, J::Array{Float64, 2}, h::Array{Float64, 1}) 
	X_after_transition = zeros(Int64, M, L)
	f1 = zeros(Float64, L*q)
	f2 = zeros(Float64, L*q, L*q)
	scale = 1.0/M	
	A = rand(0:(q-1), L)
    	E_old = E_i(q, L, 1, A, J, h)
	for m=1:M
		for k=1:k_max
			#(n_accepted, A) = Monte_Carlo_sweep(q, L, A, J, h)
			(n_accepted, A, E_old) = Monte_Carlo_sweep(E_old, q, L, A, J, h)
		end
		
		for i in 1:L
			a = X[m,i]+1 #1 to 21
			f1[(i-1)*q+a] += scale
			X_after_transition[m,i] = A[i]	
			for j in (i+1):L
				b = X[m,j] + 1 
				f2[(i-1)*q+a, (j-1)*q+b] += scale
				f2[(j-1)*q+b, (i-1)*q+a] += scale
			end
		end
	end
	return (f1, f2, X_after_transition) 
end

function pCDk_minibatch(X_persistent::Array{Int64, 2}, id_set::Array{Int64, 1}, k_max::Int64, M::Int64, q::Int64, L::Int64, J::Array{Float64, 2}, h::Array{Float64, 1}) 
	X_after_transition = copy(X_persistent) 
	f1 = zeros(Float64, L*q)
	f2 = zeros(Float64, L*q, L*q)
	n_batch = size(id_set, 1)
	scale = 1.0/(n_batch)	
	
	A = rand(0:(q-1), L)
    	E_old = E_i(q, L, 1, A, J, h)
	for n=1:n_batch
		m = id_set[n]	
		A = X_persistent[m, :] 
		for k=1:k_max
			#(n_accepted, A) = Monte_Carlo_sweep(q, L, A, J, h)
			(n_accepted, A, E_old) = Monte_Carlo_sweep(E_old, q, L, A, J, h)
		end
		
		for i in 1:L
			a = X_persistent[m,i]+1 #1 to 21
			f1[(i-1)*q+a] += scale
			X_after_transition[m,i] = A[i]	
			for j in (i+1):L
				b = X_persistent[m,j] + 1 
				f2[(i-1)*q+a, (j-1)*q+b] += scale
				f2[(j-1)*q+b, (i-1)*q+a] += scale
			end
		end
	end
	return (f1, f2, X_after_transition) 
end


function pCDk_weight(X::Array{Int64, 2}, k_max::Int64, M::Int64, W::Array{Float64,1}, q::Int64, L::Int64, J::Array{Float64, 2}, h::Array{Float64, 1}) 
	X_after_transition = zeros(Int64, M, L)
	f1 = zeros(Float64, L*q)
	f2 = zeros(Float64, L*q, L*q)
	M_eff = sum(W)	
	scale = 1.0/M_eff

	A = rand(0:(q-1), L)
	for m=1:M
		for k=1:k_max
			(n_accepted, A) = Monte_Carlo_sweep(q, L, A, J, h)
		end
		
		for i in 1:L
			a = X[m,i]+1 #1 to 21
			f1[(i-1)*q+a] += W[m] * scale
			X_after_transition[m,i] = A[i]	
			for j in (i+1):L
				b = X[m,j] + 1 
				f2[(i-1)*q+a, (j-1)*q+b] += W[m] * scale
				f2[(j-1)*q+b, (i-1)*q+a] += W[m] * scale
			end
		end
	end
	return (f1, f2, X_after_transition) 
end

function gradient_ascent_sparse( lambda_h::Float64, lambda_J::Float64, reg_h::Float64, reg_J::Float64,  f1_1::Array{Float64,1},  f1_2::Array{Float64,1},  f2_2::Array{Float64,2},  f2_msa::Array{Float64,2}, J::Array{Float64, 2},Jfilter::Array{Float64, 2}, h::Array{Float64, 1})
	C2 = f2_2 - f1_2*f1_2'
	C_msa= f2_msa - f1_1*f1_1'
	
	dh = f1_1 - f1_2
	#dJ = (C_msa - C2 ) .* Jfilter 
	dJ = (f2_msa - f2_2 ) .* Jfilter 
	h = h * (1.0 - reg_h) + lambda_h * dh   
	J =  J  + (lambda_J * dJ - lambda_J * J ).* Jfilter  
	
	#c1vec = reshape(C_msa, (L*L*q*q))
	#c2vec = reshape(C2, (L*L*q*q))
	c1vec = zeros(Float64, Int64(L*(L-1)*q*q/2) )
	c2vec = zeros(Float64, Int64(L*(L-1)*q*q/2) )
	pos = 0
	for i = 1:L
	    # id[(i-1)*q .+ (1:q), (i-1)*q .+ (1:q)] .= withdiag
	    for j = (i+1):L
	        # id[(j-1)*q .+ (1:q),(i-1)*q .+ (1:q)] .= true
	        c1vec[pos .+ (1:q^2)] .= vec(C_msa[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)])
	        c2vec[pos .+ (1:q^2)] .= vec(C2[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)])
	        pos += q^2
	    end
	end
	cc = Statistics.cor(c1vec,c2vec)
	cslope = linreg(c1vec,c2vec)[2]
	froc = LinearAlgebra.norm(c1vec - c2vec)

	return (J,h,
		sqrt(sum(dh.^2)), sqrt(sum(dJ.^2)), 
		cc,cslope,froc) 
end


function gradient_ascent( lambda_h::Float64, lambda_J::Float64, reg_h::Float64, reg_J::Float64,  f1_1::Array{Float64,1},  f1_2::Array{Float64,1},  f2_2::Array{Float64,2},  C_msa::Array{Float64,2}, J::Array{Float64, 2}, h::Array{Float64, 1})
	#C1 = f2_1 - f1_1*f1_1' # C1 is only possitive.
	
	C2 = f2_2 - f1_2*f1_2'
	
	dh = f1_1 - f1_2
	dJ = f2_1 - f2_2
	#dJ = C_msa - C2 
	h = h * (1.0 - reg_h) + lambda_h * dh   
	J = J * (1.0 - reg_J) + lambda_J * dJ 
	
	# The computational cost of the gradient_ascent is not significant compare with MC. 
	for i in 1:L
		for j in (1+i):L
			for a in 1:q
				for b in 1:q


				end
			end
		end
	end
	
	#c1vec = reshape(C_msa, (L*L*q*q))
	#c2vec = reshape(C2, (L*L*q*q))
	c1vec = zeros(Float64, Int64(L*(L-1)*q*q/2) )
	c2vec = zeros(Float64, Int64(L*(L-1)*q*q/2) )
	pos = 0
	for i = 1:L
	    # id[(i-1)*q .+ (1:q), (i-1)*q .+ (1:q)] .= withdiag
	    for j = (i+1):L
	        # id[(j-1)*q .+ (1:q),(i-1)*q .+ (1:q)] .= true
	        c1vec[pos .+ (1:q^2)] .= vec(C_msa[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)])
	        c2vec[pos .+ (1:q^2)] .= vec(C2[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)])
	        pos += q^2
	    end
	end
	cc = Statistics.cor(c1vec,c2vec)
	cslope = linreg(c1vec,c2vec)[2]
	froc = LinearAlgebra.norm(c1vec - c2vec)

	return (J,h,
		sqrt(sum(dh.^2)), sqrt(sum(dJ.^2)), 
		cc,cslope,froc) 
end


function get_statistics_BM(L::Int64, q::Int64,  n_sample::Int64,  n_weight::Int64, T_eq::Int64, J::Array{Float64, 2}, h::Array{Float64, 1})
	A = rand(0:(q-1), L)	
	X_output = zeros(Int64, n_sample, L)
    E_old = E_i(q, L, 1, A, J, h)
    E_vec = zeros(T_eq)
    accepted_vec = zeros(T_eq)
	for m=1:T_eq
		(n_accepted, A, E_old) = Monte_Carlo_sweep(E_old, q, L, A, J, h)
		E_vec[m] = E_old
		accepted_vec[m] = n_accepted
	end	
	for m=1:n_sample
		for t=1:n_weight
			(n_accepted, A, E_old) = Monte_Carlo_sweep(E_old, q, L, A, J, h)
		end
		for i=1:L
			X_output[m,i] = A[i]	
		end
	end
	return X_output, E_vec, accepted_vec
end

function output_statistics(t::Int64, L::Int64, q::Int64,  n_sample::Int64, h::Array{Float64, 1}, J::Array{Float64, 2})
	fname_out = "./statistics-t-" *string(t)* ".txt"
	fout = open(fname_out, "w")
	T_eq, T_aut = 1000, 30	
	
	A = rand(0:(q-1), L)	
    	E_old = E_i(q, L, 1, A, J, h)
	for m=1:T_eq
		#(n_accepted, A) = Monte_Carlo_sweep(q, L, A, J, h)
		(n_accepted, A, E_old) = Monte_Carlo_sweep(E_old, q, L, A, J, h)
	end	
	#X_output = zeros(Int64, n_sample, L)	
	
	for m=1:n_sample
		for t=1:T_aut
			#(n_accepted, A) = Monte_Carlo_sweep(q, L, A, J, h)
			(n_accepted, A, E_old) = Monte_Carlo_sweep(E_old, q, L, A, J, h)
		end
		
		for i=1:L
			#X_output[m,i] = A[i]	
			print(fout,A[i], " ")
		end
		println(fout,"")
	end
	close(fout)
end

function output_statistics(fname_out::String, L::Int64, q::Int64,  n_sample::Int64, h::Array{Float64, 1}, J::Array{Float64, 2})
	fout = open(fname_out, "w")
	T_eq, T_aut = 1000, 30	
	
	A = rand(0:(q-1), L)	
    	E_old = E_i(q, L, 1, A, J, h)
	for m=1:T_eq
		#(n_accepted, A) = Monte_Carlo_sweep(q, L, A, J, h)
		(n_accepted, A, E_old) = Monte_Carlo_sweep(E_old, q, L, A, J, h)
	end	
	#X_output = zeros(Int64, n_sample, L)	
	
	for m=1:n_sample
		for t=1:T_aut
			#(n_accepted, A) = Monte_Carlo_sweep(q, L, A, J, h)
			(n_accepted, A, E_old) = Monte_Carlo_sweep(E_old, q, L, A, J, h)
		end
		
		for i=1:L
			#X_output[m,i] = A[i]	
			print(fout,A[i], " ")
		end
		println(fout,"")
	end
	close(fout)
end



function output_paramters(t::Int64, L::Int64, q::Int64, h::Array{Float64, 1}, J::Array{Float64,2})
	fname_out = "./parameters-t"*string(t)*".txt"
	fout = open(fname_out, "w")
	for i=1:L
		for j=(i+1):L
			for a=1:q
				for b=1:q
					println(fout, "J ", i-1, " ", j-1, " ", a-1, " ", b-1, " ", J[(i-1)*q+a, (j-1)*q+b])
				end
			end
		end
	end
	for i=1:L
		for a=1:q
			println(fout, "h ", i-1, " ", a-1, " ", h[(i-1)*q+a])
		end
	end
	
	close(fout)
end

# Run two MCMCs q1 and q2, c :=<q1,q2>
# <q(0),q(t^*)> ~ c
function Test_Autocorrelations(L::Int64, q::Int64, T_eq::Int64, n_max::Int64, J::Array{Float64, 2}, h::Array{Float64, 1})
	av_overlap = 0.0
    @show L, q, T_eq, n_max
    @show size(J)
    @show size(h)
    for n in 1:n_max
	A1 = rand(0:(q-1), L)
        E_old = E_i(q, L, 1, A1, J, h)
        for m=1:T_eq
            i = rand(1:L)
            (n_accepted, A1, E_old) = Metropolis_Hastings(E_old, q, L, i, A1, J, h)
        end        
        
	A2 = rand(0:(q-1), L)
        E_old = E_i(q, L, 1, A2, J, h)
        for m=1:T_eq
            i = rand(1:L)
            (n_accepted, A2, E_old) = Metropolis_Hastings(E_old, q, L, i, A2, J, h)
        end        
	overlap = sum(kr.(A1,A2))
        av_overlap += overlap
    end
    av_overlap /= n_max

    A0_vec = rand(0:(q-1), (n_max, L))
    Auto_corr = zeros(T_eq)
    E_old_vec = zeros(n_max)
    
    #--- This is needed for the recurrsive function ---#
    for n in 1:n_max
        E_old_vec[n] = E_i(q, L, 1, A0_vec[n,:], J, h)
    end
    
    A_vec = copy(A0_vec)
    for n in 1:n_max
    	for t in 1:T_eq
            i = rand(1:L)
	    (n_accepted, A_vec[n,:], E_old_vec[n]) = Metropolis_Hastings(E_old_vec[n], q, L, i, A_vec[n,:], J, h)
            #@show sum(kr.(A_vec[n,:], A0_vec[n,:]))
            Auto_corr[t] += sum(kr.(A0_vec[n,:], A_vec[n,:]))
        end
    end
    Auto_corr ./= n_max
    
    Plots.plot(Auto_corr, label="overlap between <q(t)q(0)>")
    p1 = Plots.plot!(av_overlap*ones(size(Auto_corr)), label="overlap between independent <q1q2>")
    
    return Auto_corr, av_overlap, p1
end

