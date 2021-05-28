function read_dist(fname,q,L,id1,id2,id3) 
        X = readdlm(fname)    
        @show (n_max, n_ele) = size(X)
        contact_out = zeros(L,L)

        for n in 1:n_max
            i,j,value  = map(Int64, X[n,id1]),  map(Int64, X[n,id2]), map(Float64, X[n,id3])
            contact_out[i,j] = value
            contact_out[j,i] = value
        end
        return contact_out
end

function read_contact(fname,q,L,id1,id2,id3,dist) 
        X = readdlm(fname)    
        @show (n_max, n_ele) = size(X)
        contact_out = zeros(L,L)

        for n in 1:n_max
            i,j,value  = map(Int64, X[n,id1]),  map(Int64, X[n,id2]), map(Float64, X[n,id3])
            if(value<dist)
                contact_out[i,j] = 1
                contact_out[j,i] = 1
            end
        end
        return contact_out
end

function output_artificial_couplings(q, L, mean_v, var_v, density, cont_mat)
    coupling = zeros(q*L, q*L)
    for i in 1:L
        for j in (i+1):L
            for a in 1:q
            for b in 1:q
              	#Assume nonzero entrories are only diagonal. 
		if(rand()<density)
			value = mean_v*cont_mat[i,j] + var_v*randn()
			
			coupling[km(i,a,q), km(j,b,q)] = value
			coupling[km(j,b,q), km(i,a,q)] = value
            	end
            end
            end
        end
    end
    return coupling
end

#
#Generate local field so that  \sum_a f_i(a) = 1
#    q=# of states
#    L=# of residues
#    prob=frequencies of states (\sum_aa prob[a] = 10) !!! NOT 1. and prob[a] \in Integer
#    eps=pseudo count (>0)
#    Note, frequencies of amino acids are based on the prob that is permuted at each residue positions
#    Note, 
#
function output_artificial_field(rng,q,L,eps)
    #h_i(a) = log(f_i(a) + eps)
    h = zeros(q*L)
    for i in 1:L
	    freq = rand(q)
	    freq /= sum(freq)
        for a in 1:q
		h[km(i,a,q)] = log( (1-eps)*freq[a] + eps/q )
        end
    end
    return h
end
  
function main_artificial_coupling_field(q,L, mean_v, var_v, density, fname_cotact, fname_out)
    fout = open(fname_out, "w")
    cont = read_contact(fname_cotact, q, L, 1,2,4,0.1,6.0);
    #mean_v, var_v, density = 1.5, 0.2, 0.3
    amplitude_h = 10
    J = output_artificial_couplings(q, L, mean_v, var_v, density, cont);
    h = output_artificial_field(rng, q, L, 0.01);
    for i in 1:L
        for j in (i+1):L
            for a in 1:q
                for b in 1:q
                    println(fout, "J ", i-1, " ", j-1, " ", a-1, " ", b-1, " ", J[km(i,a,q), km(j,b,q)])
                end
            end
        end
    end
    for i in 1:L
        for a in 1:q
            println(fout,"h ", i-1, " ", a-1, " ", h[km(i,a,q)])
        end
    end
    close(fout)
end

function write_msa(X, fname)
    fout = open(fname, "w")
    N,L = size(X)
    for n in 1:N
        for i in 1:L
            print(fout, X[n,i], " ")
        end
        println(fout,"")
    end
    close(fout)
end


# n_branch = # of binary branch of MC evlutions
function get_correlated_data(q,L,J,h,T_eq, T_aut, n_branch)
	#---------- Equilibration ---------#
	A = rand(0:(q-1), L); A0 = copy(A)	
	overlap_vec = zeros(T_eq)
	E_old = E_i(q, L, 1, A, J, h)
	for t=1:T_eq
	    i = rand(1:L)
	    (n_accepted, A, E_old) = Metropolis_Hastings(E_old, q, L, i, A, J, h)
		#(n_accepted, A) = Monte_Carlo_adaptive(L,L, A, J, h)# it is setted nMC=L
	    overlap = sum(kr.(A0,A))
	    overlap_vec[t] = overlap
	end
	#----------------------------------#

	Data = zeros(Int64, 2^n_branch, L)
	Data_old = zeros(Int64, 2^n_branch, L)
	Data_old[1,:] = copy(A)

	leaf_id = 1
	for n in 1:n_branch
		#global Data_old, Data	
		leaf_id = 2^(n-1)
		for id in 1:leaf_id		
			#------- Sample1. ------#	
			A = Data_old[id,:]
			E_old = E_i(q, L, 1, A, J, h)
			for t=1:T_aut
			    i = rand(1:L)
			    (n_accepted, A, E_old) = Metropolis_Hastings(E_old, q, L, i, A, J, h)
				#(n_accepted, A) = Monte_Carlo_adaptive(L,nMC, A, J, h)
			end		
			Data[2*id-1,:] = copy(A)
			
			#------- Sample2. ------#	
			A = Data_old[id,:]
			E_old = E_i(q, L, 1, A, J, h)
			for t=1:T_aut
			    i = rand(1:L)
			    (n_accepted, A, E_old) = Metropolis_Hastings(E_old, q, L, i, A, J, h)
			end		
			Data[2*id,:] = copy(A) 
		end
		Data_old = copy(Data)	
	end
	return (Data, overlap_vec)
end

