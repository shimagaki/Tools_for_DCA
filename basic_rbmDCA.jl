function E_i_hidden(q::Int64, i::Int64, P::Int64, L::Int64,  
		    h::Array{Float64, 1}, xi::Array{Float64,2}, J::Array{Float64,2}, Jfiliter::Array{Int64,2},
		    A::Array{Int64, 1}, H::Array{Float64, 1})
	e_i = 0.0
	a = A[i]+1	
	for mu=1:P
		#e_i += - xi[(i-1)*P+mu, a] * H[mu]
		e_i += - xi[mu, km(i,a,q)] * H[mu]
	end	
	for j=1:L
		b = A[j]+1	
		e_i += - J[km(i,a,q), km(j,b,q)] * Jfiliter[km(i,a,q), km(j,b,q)]
	end
	
	e_i += - h[km(i,a,q)]
	return e_i 
end

function E_i_hidden(q::Int64, i::Int64, a::Int64, P::Int64, 
		    h::Array{Float64, 1}, xi::Array{Float64,2}, 
		    H::Array{Float64, 1})
	e_i = 0.0
	
	for mu=1:P
		e_i += - xi[mu, km(i,a,q)] * H[mu]
	end	
	e_i += - h[km(i,a,q)]
	return e_i 
end


function E_i_hidden_diff(i::Int64, a_old::Int64, a_prop::Int64, P::Int64, L::Int64,  
		    h::Array{Float64, 1}, xi::Array{Float64,2}, J::Array{Float64,2}, Jfiliter::Array{Int64,2},
		    A::Array{Int64, 1}, H::Array{Float64, 1})
	e_i = 0.0
	for mu=1:P
		#e_i += - (xi[(i-1)*P+mu, a_prop+1] - xi[(i-1)*P+mu, a_old+1]) * H[mu]
		e_i += - (xi[mu, km(i,a_prop+1,q)] - xi[mu, km(i,a_old+1, q)]) * H[mu]
	end	
	for j=1:L
		if(j!=i)	
			b = A[j]	
			e_i += - (J[(i-1)*q+a_prop+1, (j-1)*q + b+1] * Jfiliter[(i-1)*q+a_prop+1, (j-1)*q + b+1]
				  - J[(i-1)*q+a_old+1, (j-1)*q + b+1] * Jfiliter[(i-1)*q+a_old+1, (j-1)*q + b+1]) 
		end
	end
	e_i += - (h[(i-1)*q+a_prop+1] - h[(i-1)*q+a_old+1]) 
	return e_i 
end

function E_i_hidden_diff(i::Int64, a_old::Int64, a_prop::Int64, P::Int64, L::Int64,  
		    h::Array{Float64, 1}, xi::Array{Float64,2},
		    A::Array{Int64, 1}, H::Array{Float64, 1})
	e_i = 0.0
	for mu=1:P
		e_i += - (xi[mu, km(i,a_propm+1,q)] - xi[mu, km(i,a_old+1, q)]) * H[mu]
	end	
	e_i += - (h[(i-1)*q+a_prop+1] - h[(i-1)*q+a_old+1]) 
	return e_i 
end

function sampling_visible(q::Int64, L::Int64, P::Int64, 
			  H::Array{Float64,1}, 
			  h::Array{Float64,1}, xi::Array{Float64, 2})
	A_return = zeros(Int, L)	
	for i=1:L
		e_i_hidden = zeros(q)
		for a=1:q
			e_i_hidden[a] =  E_i_hidden(q, i, a, P, h, xi, H)
		end
		weight = exp.(-e_i_hidden)
		w = Weights(weight)
		A_return[i] = sample(w) - 1
	end
	return A_return 
end


function sampling_visible_MH(q::Int64, L::Int64, P::Int64, 
			    A::Array{Int64, 1}, H::Array{Float64,1},
			    J::Array{Float64, 2}, Jfiliter::Array{Int64, 2},
			    h::Array{Float64,1}, xi::Array{Float64, 2})
	for i=1:L
		a_old = A[i]	
		a_prop = rand( vcat( 0:(a_old-1), (a_old+1):(q-1) ) )
		
		e_i_hidden =  E_i_hidden_diff(i, a_old, a_prop, P, L, h, xi, J, Jfiliter, A, H)
		r = exp(-e_i_hidden)
		if(rand()<r)
			A[i]=a_prop
		end
	end
	return A
end

function sampling_hidden(P::Int64, L::Int64, 
			 A::Array{Int64,1}, xi::Array{Float64,2 })
	H0 = zeros(P)
	for mu=1:P
		for i=1:L
			H0[mu] += xi[mu, km(i,A[i]+1, q)] / L
		end
	end
	#H0 = ones(P)	
	return H0 +1.0/sqrt(L) * randn(P)
end

function pCDk_rbm_bm(q::Int64, L::Int64, P::Int64, 
	      M::Int64, k_max::Int64, 
	      J::Array{Float64, 2}, Jfiliter::Array{Int64, 2},
	      h::Array{Float64, 1},xi::Array{Float64, 2},  
	      X::Array{Int64, 2}, X_persistent::Array{Int64,2}) 
	
	X_after_transition = zeros(Int64, M, L)
	f1 = zeros(Float64, L*q)
	f2 = zeros(Float64, L*q, L*q)
	psi_data = zeros(Float64, P, q*L)
	psi_model = zeros(Float64, P, q*L)
	scale = 1.0/M
	A_model = zeros(Int64, L); H_model=zeros(P); H_data=zeros(P) #These are necessary since these are local variables and otherwise you cannot use out of for scope
	H_data_mean = zeros(P)	
	H_model_mean = zeros(P)
	for m=1:M
		#positive-term
		H_data = sampling_hidden(P,L, X[m,:],xi)
		#H_data_mean = H_data_mean + H_data
		#negative-term
		A_model = X_persistent[m,:]
		for k=1:k_max
			H_model = copy(sampling_hidden(P,L,A_model,xi)) 
			A_model = copy(sampling_visible_MH(q,L,P, A_model, H_model, J, Jfiliter, h, xi)) 
		end
		#H_model_mean = H_model_mean + H_model	
		
		for i in 1:L
			a = A_model[i]+1
			f1[(i-1)*q+a] += scale
			X_after_transition[m,i] = A_model[i]	
			for mu in 1:P
				psi_data[mu, km(i,X[m,i]+1,q)] += H_data[mu] * scale
				psi_model[mu, km(i,a,q)] += H_model[mu] * scale
			end
			
			for j in (i+1):L
				b = A_model[j]+1 
				f2[(i-1)*q+a, (j-1)*q+b] += scale
				f2[(j-1)*q+b, (i-1)*q+a] += scale
			end
		end
	end
	#@show(H_data_mean/M)
	#@show(H_model_mean/M)
	return (f1, f2, psi_data, psi_model, X_after_transition) 
end

function pCDk_rbm(q::Int64, L::Int64, P::Int64, 
	      M::Int64, k_max::Int64, 
	      h::Array{Float64, 1},xi::Array{Float64, 2},  
	      X::Array{Int64, 2}, X_persistent::Array{Int64,2}) 
	
	X_after_transition = zeros(Int64, M, L)
	f1 = zeros(Float64, L*q)
	f2 = zeros(Float64, L*q, L*q)
	psi_data = zeros(Float64, P, q*L)
	psi_model = zeros(Float64, P, q*L)
	scale = 1.0/M
	A_model = zeros(Int64, L); H_model=zeros(P); H_data=zeros(P) #These are necessary since these are local variables and otherwise you cannot use out of for scope
	H_data_mean = zeros(P)	
	H_model_mean = zeros(P)
	for m=1:M
		#positive-term
		H_data = sampling_hidden(P,L, X[m,:],xi)
		#H_data_mean = H_data_mean + H_data
		#negative-term
		A_model = X_persistent[m,:]
		for k=1:k_max
			H_model = copy(sampling_hidden(P,L,A_model,xi)) 
			A_model = copy(sampling_visible(q,L,P, H_model, h, xi)) 
		end
		#H_model_mean = H_model_mean + H_model	
		
		for i in 1:L
			a = A_model[i]+1
			f1[km(i,a,q)] += scale
			X_after_transition[m,i] = A_model[i]	
			for mu in 1:P
				psi_data[mu, km(i,X[m,i]+1, q)] += H_data[mu] * scale
				psi_model[mu, km(i,a,q)] += H_model[mu] * scale
			end
			
			for j in (i+1):L
				b = A_model[j]+1 
				f2[km(i,a,q), km(j,b,q)] += scale
				f2[km(j,b,q), km(i,a,q)] += scale
			end
		end
	end
	return (f1, f2, psi_data, psi_model, X_after_transition) 
end

function pCDk_rbm_minibatch(q::Int64, L::Int64, P::Int64, 
	      M::Int64, k_max::Int64,
	      id_set::Array{Int64, 1},
	      h::Array{Float64, 1},xi::Array{Float64, 2},  
	      X::Array{Int64, 2}, 
	      X_persistent::Array{Int64,2}) 
	
	X_after_transition = copy(X_persistent) #It should use the X_persistent?  
	f1 = zeros(Float64, L*q)
	f2 = zeros(Float64, L*q, L*q)
	psi_data = zeros(Float64, P, q*L)
	psi_model = zeros(Float64, P, q*L)
	n_batch = size(id_set, 1)	
	scale = 1.0/n_batch
	A_model = zeros(Int64, L); H_model=zeros(P); H_data=zeros(P) #These are necessary since these are local variables and otherwise you cannot use out of for scope
	H_data_mean = zeros(P)	
	H_model_mean = zeros(P) 
	for n=1:n_batch
		m = id_set[n]
		#positive-term
		H_data = sampling_hidden(P,L, X[m,:],xi)
		#H_data_mean = H_data_mean + H_data
		#negative-term
		A_model = X_persistent[m,:]
		for k=1:k_max
			H_model = copy(sampling_hidden(P,L,A_model,xi)) 
			A_model = copy(sampling_visible(q,L,P, H_model, h, xi)) 
		end
		#H_model_mean = H_model_mean + H_model	
		
		for i in 1:L
			a = A_model[i]+1
			f1[(i-1)*q+a] += scale
			X_after_transition[m,i] = A_model[i]	
			for mu in 1:P
				psi_data[mu, km(i, X[m,i]+1, q)] += H_data[mu] * scale
				psi_model[mu, km(i, a, q)] += H_model[mu] * scale
			end
			
			for j in (i+1):L
				b = A_model[j]+1 
				f2[(i-1)*q+a, (j-1)*q+b] += scale
				f2[(j-1)*q+b, (i-1)*q+a] += scale
			end
		end
	end
	
	return (f1, f2, psi_data, psi_model, X_after_transition) 
end

function pCDk_rbm_weight(q::Int64, L::Int64, P::Int64, 
	      M::Int64, k_max::Int64, 
	      w::Array{Float64,1},
	      h::Array{Float64, 1},xi::Array{Float64, 2},  
	      X::Array{Int64, 2}, X_persistent::Array{Int64,2}) 
	
	X_after_transition = zeros(Int64, M, L)
	f1 = zeros(Float64, L*q)
	f2 = zeros(Float64, L*q, L*q)
	psi_data = zeros(Float64, P, q*L)
	psi_model = zeros(Float64, P, q*L)
	
	M_eff=sum(w)
	scale = 1.0/Meff
	
	A_model = zeros(Int64, L); H_model=zeros(P); H_data=zeros(P) #These are necessary since these are local variables and otherwise you cannot use out of for scope
	H_data_mean = zeros(P)
	H_model_mean = zeros(P)
	for m=1:M
		#positive-term
		H_data = sampling_hidden(P,L, X[m,:],xi)
		#H_data_mean = H_data_mean + H_data
		#negative-term
		A_model = X_persistent[m,:]
		for k=1:k_max
			H_model = copy(sampling_hidden(P,L,A_model,xi)) 
			A_model = copy(sampling_visible(q,L,P, H_model, h, xi)) 
		end
		#H_model_mean = H_model_mean + H_model	
		
		for i in 1:L
			a = A_model[i]+1
			f1[(i-1)*q+a] += scale
			X_after_transition[m,i] = A_model[i]	
			for mu in 1:P
				psi_data[mu, km(i,X[m,i]+1,q)] += H_data[mu] * w[m] * scale
				psi_model[mu, km(i,a,q)] += H_model[mu] * w[m] * scale
			end
			
			for j in (i+1):L
				b = A_model[j]+1 
				f2[km(i,a,q), km(j,b,q)] +=  w[m] * scale
				f2[km(j,b,q), km(i,a,q)] +=  w[m] * scale
			end
		end
	end
	return (f1, f2, psi_data, psi_model, X_after_transition) 
end

function pCDk_rbm_weight_minbatch(q::Int64, L::Int64, P::Int64, 
	      M::Int64, k_max::Int64, 
	      id_set::Array{Int64, 1},
	      w::Array{Float64,1},
	      h::Array{Float64, 1},xi::Array{Float64, 2},  
	      X::Array{Int64, 2}, X_persistent::Array{Int64,2}) 
	
	X_after_transition = copy(X) 
	f1 = zeros(Float64, L*q)
	f2 = zeros(Float64, L*q, L*q)
	psi_data = zeros(Float64, P, q*L)
	psi_model = zeros(Float64, P, q*L)
	
	M_eff=sum(w)
	
	n_batch = size(id_set, 1)	
	scale = 1.0/n_batch
	
	A_model = zeros(Int64, L); H_model=zeros(P); H_data=zeros(P) #These are necessary since these are local variables and otherwise you cannot use out of for scope
	H_data_mean = zeros(P)
	H_model_mean = zeros(P)
	for n=1:n_batch
		m = id_set[n]
		#positive-term
		H_data = sampling_hidden(P,L, X[m,:],xi)
		#H_data_mean = H_data_mean + H_data
		#negative-term
		A_model = X_persistent[m,:]
		for k=1:k_max
			H_model = copy(sampling_hidden(P,L,A_model,xi)) 
			A_model = copy(sampling_visible(q,L,P, H_model, h, xi)) 
		end
		#H_model_mean = H_model_mean + H_model	
		
		for i in 1:L
			a = A_model[i]+1
			f1[(i-1)*q+a] += scale
			X_after_transition[m,i] = A_model[i]	
			for mu in 1:P
				psi_data[mu, km(i,X[m,i]+1,q)] += H_data[mu] * w[m] * scale
				psi_model[mu, km(i,a,q)] += H_model[mu] * w[m] * scale
			end
			
			for j in (i+1):L
				b = A_model[j]+1 
				f2[(i-1)*q+a, (j-1)*q+b] +=  w[m] * scale
				f2[(j-1)*q+b, (i-1)*q+a] +=  w[m] * scale
			end
		end
	end
	return (f1, f2, psi_data, psi_model, X_after_transition) 
end

#----- Although it shows clear autocorrelation decay, it cannot work for learning
function PCD_rbm_site_update(q::Int64, L::Int64, P::Int64, 
	      M::Int64, k_max::Int64, 
	      h::Array{Float64, 1},xi::Array{Float64, 2},  
	      X::Array{Int64, 2}, X_persistent::Array{Int64,2}) 
	
	X_after_transition = zeros(Int64, M, L)
	f1 = zeros(Float64, L*q)
	f2 = zeros(Float64, L*q, L*q)
	psi_data = zeros(Float64, P, q*L)
	psi_model = zeros(Float64, P, q*L)
	scale = 1.0/M
	A_model = zeros(Int64, L); H_model=zeros(P); H_data=zeros(P)
	H_data_mean = zeros(P)	
	H_model_mean = zeros(P)
	for m=1:M
		#positive-term
		H_data = sampling_hidden(P,L, X[m,:],xi)
		#negative-term
		A_model = X_persistent[m,:]
        	H_model = copy(sampling_hidden(P,L,A_model,xi))
		for k=1:k_max
		    for i in 1:L
			(A_model, A_i_old) = sampling_visible_i(q,L,P,i,  H_model, A_model, h, xi)
			H_model = sampling_hidden_i(q, P, L, i, A_i_old, copy(H_model), A_model, xi)
		    end
		end
		
		for i in 1:L
			a = A_model[i]+1
			f1[km(i,a,q)] += scale
			X_after_transition[m,i] = A_model[i]	
			for mu in 1:P
				psi_data[mu, km(i,X[m,i]+1, q)] += H_data[mu]
				psi_model[mu, km(i,a,q)] += H_model[mu]
			end
			
			for j in (i+1):L
				b = A_model[j]+1 
				f2[km(i,a,q), km(j,b,q)] += scale
				f2[km(j,b,q), km(i,a,q)] += scale
			end
		end
	end
    psi_data = psi_data * scale
    psi_model = psi_model * scale
	return (f1, f2, psi_data, psi_model, X_after_transition) 
end

function gradient_ascent(q::Int64, L::Int64,P::Int64,  
			 lambda_h::Float64, lambda_xi::Float64, 
			 reg_h::Float64, reg_xi::Float64,  
			 f1_data::Array{Float64,1},  f1_model::Array{Float64,1}, 
			 f2_data::Array{Float64,2}, f2_model::Array{Float64,2},  
			 psi_data::Array{Float64,2},  psi_model::Array{Float64,2}, 
			 h::Array{Float64, 1}, xi::Array{Float64, 2})
	#C1 = f2_1 - f1_1*f1_1' # C1 is only possitive.
	C_data = f2_data - f1_data*f1_data'
	C_model = f2_model - f1_model*f1_model'
	
	#It works. 	
	#reg_h, reg_xi = 1e-3, 1e-3 
	
	dh = f1_data - f1_model
	dxi = psi_data - psi_model
	dh2 = lambda_h*dh-reg_h*h
	dxi2 = lambda_xi*dxi-reg_xi*xi 
	
	h = h * (1.0 - reg_h*lambda_h) + lambda_h * dh   
	xi = xi * (1.0 - reg_xi*lambda_xi) + lambda_xi * dxi
	
	"""
	reg_vector = zeros(P)
	reg_mat = zeros(L*P, q)	
	for mu in 1:P
		for i in 1:L
			for a in 1:q
				reg_vector[mu] += abs(xi[(i-1)*P+mu, a])
			end
		end
	end
	for mu in 1:P
		for i in 1:L
			for a in 1:q
				reg_mat[(i-1)*P+mu, a] += sign(xi[(i-1)*P+mu, a]) * reg_vector[mu]
			end
		end
	end
	xi = xi + lambda_xi * dxi - reg_xi*reg_mat 
	"""

	c1vec = zeros(Float64, Int64(L*(L-1)*q*q/2) )
	c2vec = zeros(Float64, Int64(L*(L-1)*q*q/2) )
	pos = 0 
	for i in 1:L
		for j in (i+1):L
		    c1vec[pos .+ (1:q^2)] .= vec(C_data[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)])
		    c2vec[pos .+ (1:q^2)] .= vec(C_model[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)])
		    pos += q^2
		end
    	end
	
	cc = Statistics.cor(c1vec,c2vec)
	cslope = linreg(c1vec,c2vec)[2]
	froc = LinearAlgebra.norm(c1vec - c2vec)
	return (xi,h,
		sqrt(sum(dh.^2)), sqrt(sum(dxi.^2)), 
		sqrt(sum(dh2.^2)), sqrt(sum(dxi2.^2)), 
	        cc,cslope,froc) 
end

function gradient_ascent_l1(q::Int64, L::Int64,P::Int64,  
			 lambda_h::Float64, lambda_xi::Float64, 
			 reg_h::Float64, reg_xi::Float64,  
			 f1_data::Array{Float64,1},  f1_model::Array{Float64,1}, 
			 f2_data::Array{Float64,2}, f2_model::Array{Float64,2},  
			 psi_data::Array{Float64,2},  psi_model::Array{Float64,2}, 
			 h::Array{Float64, 1}, xi::Array{Float64, 2})
	#C1 = f2_1 - f1_1*f1_1' # C1 is only possitive.
	C_data = f2_data - f1_data*f1_data'
	C_model = f2_model - f1_model*f1_model'
	
	#It works. 	
	#reg_h, reg_xi = 1e-3, 1e-3 
	
	dh = f1_data - f1_model
	dxi = psi_data - psi_model
	dh2 = lambda_h*dh-reg_h*h
	dxi2 = lambda_xi*dxi-reg_xi*xi 
	
	h = h * (1.0 - reg_h*lambda_h) + lambda_h * dh   
	#L1 reg	
	xi = xi - reg_xi*lambda_xi*sign.(xi) + lambda_xi * dxi
	
	"""
	reg_vector = zeros(P)
	reg_mat = zeros(L*P, q)	
	for mu in 1:P
		for i in 1:L
			for a in 1:q
				reg_vector[mu] += abs(xi[(i-1)*P+mu, a])
			end
		end
	end
	for mu in 1:P
		for i in 1:L
			for a in 1:q
				reg_mat[(i-1)*P+mu, a] += sign(xi[(i-1)*P+mu, a]) * reg_vector[mu]
			end
		end
	end
	xi = xi + lambda_xi * dxi - reg_xi*reg_mat 
	"""

	c1vec = zeros(Float64, Int64(L*(L-1)*q*q/2) )
	c2vec = zeros(Float64, Int64(L*(L-1)*q*q/2) )
	pos = 0 
	for i in 1:L
		for j in (i+1):L
		    c1vec[pos .+ (1:q^2)] .= vec(C_data[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)])
		    c2vec[pos .+ (1:q^2)] .= vec(C_model[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)])
		    pos += q^2
		end
    	end
	
	cc = Statistics.cor(c1vec,c2vec)
	cslope = linreg(c1vec,c2vec)[2]
	froc = LinearAlgebra.norm(c1vec - c2vec)
	return (xi,h,
		sqrt(sum(dh.^2)), sqrt(sum(dxi.^2)), 
		sqrt(sum(dh2.^2)), sqrt(sum(dxi2.^2)), 
	        cc,cslope,froc) 
end



function gradient_ascent(q::Int64, L::Int64,P::Int64,  
			 lambda_h::Float64, lambda_xi::Float64, lambda_J::Float64,  
			 reg_h::Float64, reg_xi::Float64, reg_J::Float64, 
			 f1_data::Array{Float64,1},  f1_model::Array{Float64,1}, 
			 f2_data::Array{Float64,2}, f2_model::Array{Float64,2},  
			 psi_data::Array{Float64,2},  psi_model::Array{Float64,2}, 
	      		 J::Array{Float64, 2}, Jfiliter::Array{Int64, 2},
			 h::Array{Float64, 1}, xi::Array{Float64, 2} )
	#C1 = f2_1 - f1_1*f1_1' # C1 is only possitive.
	C_data = f2_data - f1_data*f1_data'
	C_model = f2_model - f1_model*f1_model'
	

	dh = f1_data - f1_model
	dxi = psi_data - psi_model
	dJ = f2_data - f2_model	
	
	dh2 = lambda_h*dh-reg_h*h
	dxi2 = lambda_xi*dxi-reg_xi*xi 
	dJ2 = (lambda_J*dJ - reg_J*J) *Jfiliter

	h = h * (1.0 - reg_h*lambda_h) + lambda_h * dh   
	xi = xi * (1.0 - reg_xi*lambda_xi) + lambda_xi * dxi
	J = (J * (1.0 - reg_J*lambda_J) + lambda_J * dJ) *Jfiliter 	
	
	"""
	reg_vector = zeros(P)
	reg_mat = zeros(L*P, q)	
	for mu in 1:P
		for i in 1:L
			for a in 1:q
				reg_vector[mu] += abs(xi[(i-1)*P+mu, a])
			end
		end
	end
	for mu in 1:P
		for i in 1:L
			for a in 1:q
				reg_mat[(i-1)*P+mu, a] += sign(xi[(i-1)*P+mu, a]) * reg_vector[mu]
			end
		end
	end
	xi = xi + lambda_xi * dxi - reg_xi*reg_mat 
	"""

	c1vec = zeros(Float64, Int64(L*(L-1)*q*q/2) )
	c2vec = zeros(Float64, Int64(L*(L-1)*q*q/2) )
	pos = 0 
	for i in 1:L
		# id[(i-1)*q .+ (1:q), (i-1)*q .+ (1:q)] .= withdiag
		for j in (i+1):L
		    # id[(j-1)*q .+ (1:q),(i-1)*q .+ (1:q)] .= true
		    c1vec[pos .+ (1:q^2)] .= vec(C_data[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)])
		    c2vec[pos .+ (1:q^2)] .= vec(C_model[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)])
		    pos += q^2
		end
    	end
	
	cc = Statistics.cor(c1vec,c2vec)
	cslope = linreg(c1vec,c2vec)[2]
	froc = LinearAlgebra.norm(c1vec - c2vec)
	return (xi,h,
		sqrt(sum(dh.^2)), sqrt(sum(dxi.^2)), sqrt(sum(dJ.^2)),  
		sqrt(sum(dh2.^2)), sqrt(sum(dxi2.^2)), sqrt(sum(dJ2.^2)),  
	        cc,cslope,froc) 
end



function output_paramters(t::Int64, L::Int64, q::Int64,  h::Array{Float64,1},  f1_data::Array{Float64,1},  f1::Array{Float64,1},  xi::Array{Float64,2}, psi_data::Array{Float64,2},  psi_model::Array{Float64,2},  f1_sample::Array{Float64,1},  f2_sample::Array{Float64,2})
	fname_out = "./parameters-t-" *string(t)* "_rbm.txt"
	fout = open(fname_out, "w")
	for i=1:L
		for mu=1:P 
			for a=1:q
				println(fout, "xi ", i-1, " ", mu-1, " ", a-1, " ", xi[mu, km(i,a,q)], " ", psi_data[mu, km(i,a,q)], " ", psi_model[mu, km(i,a,q)])
			end
		end
	end
	for i=1:L
		for a=1:q
			println(fout, "h ", i-1, " ", a-1, " ", h[(i-1)*q+a], " ", f1_data[(i-1)*q+a] ," ",f1[(i-1)*q+a], " ",f1_sample[(i-1)*q+a])
		end
	end
	
	close(fout)
end

function output_statistics(t::Int64, L::Int64, P::Int64, n_sample::Int64, n_weight::Int64,  xi::Array{Float64, 2}, h::Array{Float64, 1}, f1_msa::Array{Float64, 1}, f2_msa::Array{Float64, 2}, c2_msa::Array{Float64, 2})
	fname_out = "./statistics-t-" *string(t)* "_rbm.txt"
	fout = open(fname_out, "w")
	
	fname_out_h = "./statistics-t-" *string(t)* "_rbm_hidden.txt"
	fout_h = open(fname_out_h, "w")
	
	A_model = rand(0:20, L)	
	H_model=zeros(P)
	for m=1:1000
		H_model = sampling_hidden(P,L,A_model,xi)
		A_model = sampling_visible(q,L,P, H_model, h, xi)
	end	
	X_output = zeros(Int64, n_sample, L)	
	
	for m=1:n_sample
		for t=1:n_weight
			H_model = sampling_hidden(P,L,A_model,xi)
			A_model = sampling_visible(q,L,P, H_model, h, xi)
		end
		
		for i=1:L
			X_output[m,i] = A_model[i]	
			print(fout,A_model[i], " ")
		end
		println(fout,"")

		for mu=1:P
			print(fout_h, H_model[mu], " ")
		end
		println(fout_h, "")
	
	end
	close(fout); close(fout_h)
	
	(M_eff, L, f1_out, f2_out, c2_out) = f1_f2_c2(X_output, ones(n_sample),q)	
	fname_out1 = "./frequencies-f1-t-" *string(t)* "_rbm.txt"
	fout1 = open(fname_out1, "w")
	fname_out2 = "./frequencies-f2-t-" *string(t)* "_rbm.txt"
	fout2 = open(fname_out2, "w")
	for i in 1:L
		for a in 1:q
			println(fout1, i, " ", a, " ", f1_msa[(i-1)*q+a], " ", f1_out[(i-1)*q+a])
			for j in (i+1):L
				for b in 1:q
					println(fout2, i, " ", a, " ", j, " ", b, " ", f2_msa[(i-1)*q+a, (j-1)*q+b], " ", f2_out[(i-1)*q+a, (j-1)*q+b], " ", c2_msa[(i-1)*q+a, (j-1)*q+b], " ", c2_out[(i-1)*q+a, (j-1)*q+b])
				end
			end
		end
	end
	close(fout1); close(fout2)
	
	return (f1_out, f2_out)
end

function output_statistics(fname_out::String, L::Int64, P::Int64, n_sample::Int64, n_weight::Int64,  xi::Array{Float64, 2}, h::Array{Float64, 1}, f1_msa::Array{Float64, 1}, f2_msa::Array{Float64, 2}, c2_msa::Array{Float64, 2})
	#fname_out = "./statistics-t-" *string(t)* "_rbm.txt"
	fout = open(fname_out, "w")
	
	A_model = rand(0:20, L)	
	H_model=zeros(P)
	for m=1:1000
		H_model = sampling_hidden(P,L,A_model,xi)
		A_model = sampling_visible(q,L,P, H_model, h, xi)
	end	
	X_output = zeros(Int64, n_sample, L)	
	
	for m=1:n_sample
		for t=1:n_weight
			H_model = sampling_hidden(P,L,A_model,xi)
			A_model = sampling_visible(q,L,P, H_model, h, xi)
		end
		
		for i=1:L
			X_output[m,i] = A_model[i]	
			print(fout,A_model[i], " ")
		end
		println(fout,"")

	
	end
	close(fout);
	
	(M_eff, L, f1_out, f2_out, c2_out) = f1_f2_c2(X_output, ones(n_sample),q)	
	return (f1_out, f2_out)
end



function get_statistics(L::Int64, P::Int64, n_sample::Int64, n_weight::Int64,  xi::Array{Float64, 2}, h::Array{Float64, 1})
	A_model = rand(0:20, L)	
	H_model=zeros(P)
	for m=1:1000
		H_model = sampling_hidden(P,L,A_model,xi)
		A_model = sampling_visible(q,L,P, H_model, h, xi)
	end	
	X_output = zeros(Int64, n_sample, L)	
	
	for m=1:n_sample
		for t=1:n_weight
			H_model = sampling_hidden(P,L,A_model,xi)
			A_model = sampling_visible(q,L,P, H_model, h, xi)
		end
		for i=1:L
			X_output[m,i] = A_model[i]	
		end
	end
	return X_output	
end

function get_statistics2(L::Int64, P::Int64, n_sample_each_chain::Int64, n_chain::Int64, n_weight::Int64,  xi::Array{Float64, 2}, h::Array{Float64, 1})
	A_model = rand(0:20, L)	
	H_model=zeros(P)
	X_output = zeros(Int64, n_sample, L)	
	for k in 1:n_chain	
		for m=1:300
			H_model = sampling_hidden(P,L,A_model,xi)
			A_model = sampling_visible(q,L,P, H_model, h, xi)
		end	
		
		for m=1:n_sample_each_chain
			for t=1:n_weight
				H_model = sampling_hidden(P,L,A_model,xi)
				A_model = sampling_visible(q,L,P, H_model, h, xi)
			end
			for i=1:L
				X_output[km(m,k,n_chain),i] = A_model[i]	
			end
		end
	end	
	return X_output	
end


#T_eq_rbm = 100; T_weight_rbm=3
##X_rbm = get_statistics_RBM(q, L, P, n_sample, T_weight_rbm, T_eq_rbm, 12.5*xi, h);
function get_statistics_RBM(q::Int64, L::Int64, P::Int64, n_sample::Int64, T_weight::Int64, T_eq::Int64, xi::Array{Float64, 2}, h::Array{Float64, 1})
    A_model = rand(0:(q-1), L)	
    for t in 1:T_eq
        H_model = sampling_hidden(P, L, A_model, xi)
        A_model = sampling_visible(q,L,P, H_model, h, xi)
    end
    X_output = zeros(Int64, n_sample, L)
    for n in 1:n_sample
        for t in 1:T_weight
            H_model = sampling_hidden(P, L, A_model, xi)
            A_model = sampling_visible(q,L,P, H_model, h, xi)
        end
        X_output[n,:] = copy(A_model)
    end
    return X_output
end

function get_J_h_from_xi(q::Int64, L::Int64, P::Int64, xi::Array{Float64, 2})
	#NOTE: h = h + h_xi: where h is the original h in RBM.
	J_xi = zeros(q*L,q*L)
	h_xi = zeros(q*L)
	scale = 1.0/L
	for i in 1:L
		for j in i:L
			for a in 1:q
				for b in 1:q
					temp = 0.0
					for m in 1:P	
						temp += xi[m, km(i,a,q)] * xi[m, km(j,b,q)]
					end
					temp = temp * scale
					if(i!=j)	
						J_xi[(i-1)*q+a, (j-1)*q+b] = temp 
						J_xi[(j-1)*q+b, (i-1)*q+a] = temp 
					end
					if(i==j && b==a)	
						h_xi[(i-1)*q+a] = temp 
					end
				end
			end
		end
	end
	return (J_xi, h_xi)
end


# Single site update for the hidden and visible variables. 
function sampling_visible_i(q::Int64, L::Int64, P::Int64, i::Int64,
			  H::Array{Float64,1}, A_original::Array{Int64,1}, 
			  h::Array{Float64,1}, xi::Array{Float64, 2})
	A_return = copy(A_original)
    e_i_hidden = zeros(q)
    for a=1:q
        e_i_hidden[a] =  E_i_hidden(q, i, a, P, h, xi, H)
    end
    weight = exp.(-e_i_hidden)
    w = Weights(weight)
    A_return[i] = sample(w) - 1
	return (A_return, A_original[i]) 
end

# Suppose the change is only site i
function sampling_hidden_i(q::Int64, P::Int64, L::Int64, i::Int64, A_i_old::Int64, H_old::Array{Float64,1},
			 A::Array{Int64,1}, xi::Array{Float64,2 })
	H0 = copy(H_old)
	for mu=1:P
            a = A[i]+1
            b = A_i_old+1
	    H0[mu] += (xi[mu, km(i,a,q)] - xi[mu, km(i,b,q)]) / L 
	end
	return H0 +1.0/sqrt(L) * randn(P)
end


#--- Check Autocorrelation time: This is direct comparison between
#---- the dynamics of the laten variable model and non-latent variable model.
#---- since each update, a single single variable can be updated, therefore we can 
#---- compare with the standard sweep algorithms. 
#n_max=100; T_eq = 1000
#@time (E_vec, overlap_vec, overlap) = Test_Autocorrelations_rbm(q, L, P, n_max, T_eq, xi, h);
function Test_Autocorrelations_rbm_i(q::Int64, L::Int64, P::Int64, n_max::Int64, T_eq::Int64, xi::Array{Float64, 2}, h::Array{Float64, 1})
	av_overlap = 0.0
    (J_HP,h_HP) = convert_xi2J_h(xi, q, L, P)
    h_tot = h_HP + h
    #T_eq_overlap = T_eq; # this is ad hoc

    av_overlap = 0
    for n in 1:n_max
        A_model_1 = rand(0:(q-1), L)	
        H_model_temp = sampling_hidden(P, L, A_model_1, xi)
        for t in 1:T_eq
            #H_model_temp = sampling_hidden(P, L, A_model_1, xi)
            i = rand(1:L)
            (A_model_1, A_i_old) = sampling_visible_i(q,L,P, i, H_model_temp, A_model_1, h, xi)
            H_model_temp = sampling_hidden_i(q, P, L, i, A_i_old, copy(H_model_temp), A_model_1, xi)
            #A_model_1 = sampling_visible(q,L,P, H_model_temp, h, xi)
        end
        A_model_2 = rand(0:(q-1), L)	
        H_model_temp = sampling_hidden(P, L, A_model_2, xi)
        for t in 1:T_eq
            i = rand(1:L)
            (A_model_2, A_i_old) = sampling_visible_i(q,L,P, i, H_model_temp, A_model_2, h, xi)
            H_model_temp = sampling_hidden_i(q, P, L, i, A_i_old, copy(H_model_temp), A_model_2, xi)
            #A_model_2 = sampling_visible(q,L,P, H_model_temp, h, xi)
        end
        overlap = sum(kr.(A_model_1,A_model_2))
        av_overlap += overlap 
    end
    av_overlap /= n_max
       
    E_vec = zeros(T_eq)
    overlap_vec = zeros(T_eq)
    for n in 1:n_max
        A_model = rand(0:(q-1), L)	
        A_model_0 = copy(A_model)
        H_model=zeros(P)    
        H_model = sampling_hidden(P,L,A_model,xi)
        for t=1:T_eq
            i = rand(1:L)
            (A_model, A_i_old) = sampling_visible_i(q,L,P, i, H_model, A_model, h, xi)
            H_model = sampling_hidden_i(q, P, L, i, A_i_old, copy(H_model), A_model, xi)
            E_tot = sum([E_i(q, L, i, A_model, J_HP, h_tot) for i in 1:L])
            E_vec[t] += E_tot
            overlap = sum(kr.(A_model, A_model_0))
            overlap_vec[t] += overlap
        end 
    end
    E_vec = 1.0/n_max .* E_vec
    overlap_vec = 1.0/n_max .* overlap_vec
    Plots.plot(overlap_vec, label="<q(0)q(t)>", color="blue")
    p1 = Plots.plot!(av_overlap*ones(T_eq), label="<q1q2>", color="red")
    return (E_vec, overlap_vec, av_overlap, p1) 
end
