function get_MSA(fname, q, L, th)
    X=readdlm(fname, Int)
    X_temp=[]
    n_id=1
    @show size(X)
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

    X_temp = X_temp + ones(Int, size(X_temp))
    X = zeros(Int, (M, L*q))

    for m in 1:M
        for i in 1:L
            X[m, ((i-1)*q+X_temp[m,i])] = 1
        end
    end
    X=X'
    return X
end

function get_MSA_fast(fname, q, L, th)
    X=readdlm(fname, Int)
    X_temp = zeros(Int, N, q*L)
    n_id=1
    @show (N,L1) = size(X)
    for n in 1:N
        if(rand()<th)
		X_temp[n_id, km.(1:L, X[n,:] + ones(L),q)] = ones(L) 
            	n_id+=1
        end
    end
    return X_temp[1:n_id,:]'
    #M=size(X_temp,1)

    #X_temp = X_temp + ones(Int, size(X_temp))
    #X = zeros(Int, (M, L*q))

    #for m in 1:M
    #    for i in 1:L
    #        X[m, ((i-1)*q+X_temp[m,i])] = 1
    #    end
    #end
    #X=X'
    #return X
end




function read_Xi(fname,P, L,q) 
	X = readdlm(fname)    
@show (n_max, n_ele) = size(X)
    Xi_vec =  zeros(P,q*L)

    for n in 1:n_max
        if(X[n,1]=="xi")
            i,m,a, value  = map(Int64, X[n,2]+1),  map(Int64, X[n,3]+1), map(Int64, X[n,4]+1),  map(Float64, X[n,5])
	    Xi_vec[m,km(i,a,q)] = value
        end
    end
    return Xi_vec
end

function read_genral(fname,string_key, N1,N2,N0, id_key, id1,id2,id3,id4) 
	X = readdlm(fname)    
	@show (n_max, n_ele) = size(X)
	out_vec =  zeros(N1,N0*N2)

	for n in 1:n_max
	    if(X[n,id_key]==string_key)
		    i,m,a, value  = map(Int64, X[n,id1]+1),  map(Int64, X[n,id2]+1), map(Int64, X[n,id3]+1),  map(Float64, X[n,id4])
		    out_vec[m,km(i,a,q)] = value
	    end
	end
	return out_vec
end




linreg(x, y) = hcat(fill!(similar(x), 1), x) \ y
km(i,a,q) = (i-1)*q+a
kr(i,j) = Int(Bool(i==j))

function get_MSA(fname, fname_w, q, L, th)
    X=readdlm(fname, Int)
    w=readdlm(fname_w)
    X_temp=[]
    n_id=1

    for n in 1:size(X,1) 
        r=rand()
        if(r<th*w[n])
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

    X_temp = X_temp + ones(Int, size(X_temp))
    X = zeros(Int, (M, L*q))

    for m in 1:M
        for i in 1:L
            X[m, ((i-1)*q+X_temp[m,i])] = 1
        end
    end
    X=X'
    return X
end

function get_MSA_weight(fname, fname_w, q, L, th)
    X=readdlm(fname, Int)
    w=readdlm(fname_w)
    X_temp=[]
    n_id=1

    for n in 1:size(X,1) 
        r=rand()
        if(r<th*w[n])
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

    X_temp = X_temp + ones(Int, size(X_temp))
    X = zeros(Int, (M, L*q))

    for m in 1:M
        for i in 1:L
            X[m, ((i-1)*q+X_temp[m,i])] = 1
        end
    end
    X=X'
    return X
end



#M=size(Xs[1,:])
function basic_MSA_analysis(X,q,L)
    
    n_col,M = size(X)
    f1 = sum(X,dims=2)/M;
    @show size(f1)
    f2 = X*X'/M
    @show size(f2)
    c2 = f2 - f1*f1'
    for i in 1:L
        c2[((i-1)*q+1):i*q, ((i-1)*q+1):i*q] = zeros(q,q)
    end
    Frob_c2 = zeros(L,L)
    for i in 1:L
        for j in (i+1):L
            mynorm = norm(c2[((i-1)*q+1):i*q, ((j-1)*q+1):j*q] )
            Frob_c2[i,j] = mynorm; Frob_c2[j, i] = mynorm;
        end
    end
    
    return (f1,f2,c2, Frob_c2)
end

function basic_MSA_analysis(X, w,q,L)
    
    n_col,M = size(X)
    Meff=sum(w)
    scale = 1.0/Meff

    Xweighted1 = zeros(size(X))
    Xweighted2 = zeros(size(X))
    for n in 1:M
        Xweighted1[:,n] = sqrt(w[n])*X[:,n]
        Xweighted2[:,n] = w[n]*X[:,n]
    end
    f1 = sum(Xweighted2,dims=2)/Meff;
    f2 = Xweighted1*Xweighted1'/Meff
    c2 = f2 -f1*f1'
    for i in 1:L
        c2[((i-1)*q+1):i*q, ((i-1)*q+1):i*q] = zeros(q,q)
    end
    Frob_c2 = zeros(L,L)
    for i in 1:L
        for j in (i+1):L
            mynorm = norm(c2[((i-1)*q+1):i*q, ((j-1)*q+1):j*q] )
            Frob_c2[i,j] = mynorm; Frob_c2[j, i] = mynorm;
        end
    end
    
    
    return (f1,f2,c2, Frob_c2)
end

function get_MI(f2, f1, L, q)
    MI = zeros(L,L); mi_local = 0.0
    for i in 1:L
        for j in (i+1):L
            mi_local = 0.0
            for a in 1:q
                for b in 1:q
                    #@show i, j, a, b
                    if(f2[(i-1)*q+a, (j-1)*q+b]> 0)
                        
                        #@show i,j,a,b, f2[(i-1)*q+a, (j-1)*q+b],  f1[(i-1)*q+a],  f1[(j-1)*q+b]
                        mi_local +=  f2[(i-1)*q+a, (j-1)*q+b] * log( f2[(i-1)*q+a, (j-1)*q+b] / ( f1[(i-1)*q+a] *  f1[(j-1)*q+b] ) ) 
                    end                    
                end
            end
            MI[i,j] = mi_local
            MI[j,i] = mi_local            
        end
    end
    return MI
end


function get_MI_cont_noco(MI, Dist, q,L)
    MI_cont=zeros(L*L)
    MI_noco=zeros(L*L)
    n_cont=1;n_noco=1;
    for i in 1:L
        for j in (i+1):L
            if(Dist[i,j]==1)
                MI_cont[n_cont] = MI[i,j]                        
                n_cont += 1
            end

            if(Dist[i,j]==0)
                MI_noco[n_noco] = MI[i,j]
                n_noco += 1
            end
        end
    end
    return MI_cont[1:n_cont], MI_noco[1:n_noco]
end
            
function get_Cvec(X)
    F2 = X*X' / size(X,2)
    mean_X = sum(X,dims=2) / size(X,2)
    C2 = F2 - mean_X*mean_X';
    return C2, F2, mean_X
end

function get_dist(fname_dist, L)
    #L=112#
    #
    param_dist = readdlm(fname_dist);
    Dist = zeros(L,L)
    for n in 1:size(param_dist, 1)
        i,j,value = map(Int64, param_dist[n,1]), map(Int64, param_dist[n,2]), map(Float64, param_dist[n,4])
        if(value<4)
            Dist[i,j] = 1        
            Dist[j,i] = 1
        end
    end
    return Dist
end

function contact_noncontact_cov(C, Dist, L)
    C_cont = zeros(L*q*q*L)
    C_noco = zeros(L*q*q*L)
    n_cont=1; n_noco=1;
    for i in 1:L
        for j in (i+1):L
            
            for a in 1:q
                for b in 1:q
                    if(Dist[i,j]==1)
                        C_cont[n_cont] = C[(i-1)*q+a, (j-1)*q+b]                        
                        n_cont += 1
                    end

                    if(Dist[i,j]==0)
                        C_noco[n_noco] = C[(i-1)*q+a, (j-1)*q+b]
                        n_noco += 1
                    end
                end
            end
        end
    end
    return C_cont[1:n_cont], C_noco[1:n_noco]
end

""" INPUT: fname, q,L,p, id=6: where id is the number of colum of coupling values, the defalut is J_elem
    OUTPUT: (J,xi,h) 
"""
function read_param_in_CA(fname, q,L,P, id=6)
    param = readdlm(fname)
    (n_row,n_col) = size(param)
    J = zeros(q*L, q*L); xi = zeros(P, q*L); h = zeros(q*L)
    for n in 1:n_row
        if(param[n,1]=="J")
            i,j,a,b,v = param[n, 2]+1, param[n, 3]+1, param[n, 4]+1, param[n, 5]+1, param[n,id]
            J[(i-1)*q+a, (j-1)*q+b] = v
            J[(j-1)*q+b,(i-1)*q+a ] = v
        end
        if(param[n,1]=="xi")
            i,m,a,v = param[n, 2]+1, param[n, 3]+1, param[n, 4]+1, param[n, 5]
	    xi[m, km(i,a,q)] = v
        end
        if(param[n,1]=="h")
            i,a,v = param[n, 2]+1, param[n, 3]+1, param[n, 4]
            h[(i-1)*q+a] = v
        end
    end
    return (J,xi,h)
end 

""" INPUT: (xi, q, L, P)
    OUTPUT: (J,h)
"""
function convert_xi2J_h(xi,q,L,P)
    J = zeros(q*L, q*L)
    h = zeros(q*L)
    scale = 0.5/L
    for i in 1:L
        for j in i:L
            for a in 1:q
                for b in 1:q
                    temp = 0
                    for m in 1:P
                        #temp += xi[(i-1)*P+m, a] * xi[(j-1)*P+m, b]
			temp += xi[m, km(i,a,q)] * xi[m, km(j,b,q)] 
                    end
                    temp = temp * scale
                    if(i!=j)
                        J[km(i,a,q), km(j,b,q)] = temp
                        J[km(j,b,q), km(i,a,q)] = temp
                    end
                    if(i==j && a==b)
                        h[km(i,a,q)] = temp
                    end
                end
            end
        end
    end
    return (J,h)
end

""" INPUT: (J, q, L)
    OUTPUT: (F,S), 
    where F is Frobenius norm using Lattice-Gas gauge, and
    S is the APC score.
"""
function J2S_APC(J, q, L)
    F = zeros(L,L); S = zeros(L,L)
    for i in 1:L
        for j in (i+1):L
            c_norm=0.0
            # Lattice-Gas gauge
            for a in 2:q
                for b in 2:q
                   c_norm += J[(i-1)*q+a, (j-1)*q+b]^2
                end
            end
            c_norm = sqrt(c_norm)
            F[i,j] = c_norm; F[j,i] = c_norm            
        end
    end
    F_1site = sum(F, dims=2) / L
    F_norm = sum(F_1site) / L
    S = F - F_1site * F_1site' / F_norm 
    
    return (F,S)
end


""" INPUT: (fname1, fname2, q, L, p1, p2), 
    where p1, and p2 are the number of hidden pattern respectively.
    
"""
function coupling_analysis(fname1,fname2, q, L, p1, p2) 
    (J_p1, xi_p1, h_p1) = read_param_in_CA(fname1, q, L, p1);
    (J_p2, xi_p2, h_p2) = read_param_in_CA(fname2, q, L, p2);
    (J_RBM_p1, h_RBM_p1) = convert_xi2J_h(xi_p1, q, L, p1);
    (J_RBM_p2, h_RBM_p2) = convert_xi2J_h(xi_p2, q, L, p2);
    J_p1_vec = vec(J_p1);
    J_RBM_p1_vec = vec(J_RBM_p1);
    J_p2_vec = vec(J_p2);
    J_RBM_p2_vec = vec(J_RBM_p2);
    
    Plots.scatter(J_p1_vec[:], J_RBM_p1_vec[:], xlabel="Coupling Activation",ylabel="Hopfield Coupling" , color="blue", alpha=0.3, legend=false)
    Plots.savefig("./J-CA_vs_J-HP_p"*string(p1)*".png")
    
    Plots.scatter(J_p2_vec[:], J_RBM_p2_vec[:], xlabel="Coupling Activation",ylabel="Hopfield Coupling" , color="blue", alpha=0.3, legend=false)
    Plots.savefig("./J-CA_vs_J-HP_p"*string(p2)*".png")

    (F_p1, S_p1) = J2S_APC(J_p1, q, L);
    (F_RBM_p1, S_RBM_p1) = J2S_APC(J_RBM_p1, q, L);
    (F_p2, S_p2) = J2S_APC(J_p2, q, L);
    (F_RBM_p2, S_RBM_p2) = J2S_APC(J_RBM_p2, q, L);

    fig1 = Plots.heatmap(F_RBM_p1, size=(600,500), label="Frob")
    fig2 = Plots.heatmap(S_RBM_p1, size=(600,500), label="APC")
    Plots.plot(fig1,fig2, layout=(1,2), size=(1000, 400))
    Plots.savefig("./F_S_heatmap_p"*string(p1)*"_HP.png")

    fig1 = Plots.heatmap(F_RBM_p2, size=(600,500), label="Frob")
    fig2 = Plots.heatmap(S_RBM_p2, size=(600,500), label="APC")
    Plots.plot(fig1,fig2, layout=(1,2), size=(1000, 400))
    Plots.savefig("./F_S_heatmap_p"*string(p2)*"_HP.png")
   
    fig1 = Plots.heatmap(F_p1, size=(600,500), label="Frob")
    fig2 = Plots.heatmap(S_p1, size=(600,500), label="APC")
    Plots.plot(fig1,fig2, layout=(1,2), size=(1000, 400))
    Plots.savefig("./F_S_heatmap_p"*string(p1)*"_CA.png")

    fig1 = Plots.heatmap(F_p2, size=(600,500), label="Frob")
    fig2 = Plots.heatmap(S_p2, size=(600,500), label="APC")
    Plots.plot(fig1,fig2, layout=(1,2), size=(1000, 400))
    Plots.savefig("./F_S_heatmap_p"*string(p2)*"_CA.png") 
    Plots.scatter(J_p1_vec[:], J_RBM_p1_vec[:], xlabel="Coupling Activation",ylabel="Hopfield Coupling" , color="blue", alpha=0.3, legend=false)
    Plots.savefig("./J-CA_vs_J-HP_p"*string(p1)*".png")
    
    Plots.scatter(J_p2_vec[:], J_RBM_p2_vec[:], xlabel="Coupling Activation",ylabel="Hopfield Coupling" , color="blue", alpha=0.3, legend=false)
    Plots.savefig("./J-CA_vs_J-HP_p"*string(p2)*".png")

    @show size(F_p1)
    @show size(F_p2)
    @show size(F_RBM_p1)
    @show size(F_RBM_p2)

    @show size(S_p1)
    @show size(S_p2)
    @show size(S_RBM_p1)
    @show size(S_RBM_p2)

    
    Plots.scatter(vec(F_p1), vec(F_RBM_p1), xlabel="Coupling Activation",ylabel="Hopfield Coupling" , color="blue", alpha=0.3, legend=false)
    Plots.savefig("./F-CA_vs_F-HP_p"*string(p1)*".png")
    
    Plots.scatter(vec(S_p1), vec(S_RBM_p1), xlabel="Coupling Activation",ylabel="Hopfield Coupling" , color="blue", alpha=0.3, legend=false)
    Plots.savefig("./S-CA_vs_S-HP_p"*string(p1)*".png")

    Plots.scatter(vec(F_p2), vec(F_RBM_p2), xlabel="Coupling Activation",ylabel="Hopfield Coupling" , color="blue", alpha=0.3, legend=false)
    Plots.savefig("./F-CA_vs_F-HP_p"*string(p2)*".png")
    
    Plots.scatter(vec(S_p2), vec(S_RBM_p2), xlabel="Coupling Activation",ylabel="Hopfield Coupling" , color="blue", alpha=0.3, legend=false)
    Plots.savefig("./S-CA_vs_S-HP_p"*string(p2)*".png")
    
end

function get_J_opt_Likelihood_Variation(alpha, th, q, L, f1_msa, f2_msa, f2_model)
	scale = 1.0/q^2	
	scale1 = 1.0/q	
	J_list_Int = zeros(Int64, map(Int64, q*q*L*(L-1)/2), 4)	
	J_list_Float = zeros(map(Int64, q*q*L*(L-1)/2), 5)	
	n =1 
	for i in 1:L
		for j in (i+1):L
			delta_l_of_J_block=0.0	
			MI_ele= 0.0	
			for a in 1:q
				for b in 1:q
					f_d = (1-alpha)*f2_msa[(i-1)*q+a, (j-1)*q+b]+alpha*scale 
					f_m =(1-alpha)*f2_model[(i-1)*q+a, (j-1)*q+b]+alpha*scale 
					if(f_d>th && f_m>th)	
						delta_l_of_J_block += f_d * log(f_d / f_m )	
						f_d_a_b =  (1-alpha)*f1_msa[(i-1)*q+a]*f1_msa[(j-1)*q+b]+alpha*scale	
						MI_ele += f_d * log(f_d / (f_d_a_b) )	
					end	
				end
			end
		
		
			for a in 1:q
				for b in 1:q
					f_d = (1-alpha)*f2_msa[(i-1)*q+a, (j-1)*q+b]+alpha*scale 
					f_m =(1-alpha)*f2_model[(i-1)*q+a, (j-1)*q+b]+alpha*scale 
					
					J_elem=0.0; delta_l_of_J_elem=0.0; J_block=0.0	
					if(f_d>th && f_m>th)	
						delta_l_of_J_elem = f_d * log(f_d/f_m) +(1-f_d) * log( (1-f_d) / (1-f_m) )	
						J_elem = log( (f_d * (1-f_m)) / (f_m * (1-f_d)) ) 	
						J_block = log( f_d / f_m ) 	
					end	
					#The likelihood varidation should be discussed only on the couplings that is not introduced. 
					#take into account only zeros couplings. 
				
					J_list_Int[n,1],J_list_Int[n,2] = i,j 
					J_list_Int[n,3],J_list_Int[n,4] = a,b 
					
					J_list_Float[n,1] = J_elem
					J_list_Float[n,2] = delta_l_of_J_elem
					J_list_Float[n,3] = J_block
					J_list_Float[n,4] = delta_l_of_J_block
					J_list_Float[n,5] = MI_ele 
					n += 1	
				end
			end	
				
		end
	end
	
	return (J_list_Int, J_list_Float) 
end

function output_paramters_adding_couplings(fname_out::String, L::Int64, q::Int64, P::Int64, h::Array{Float64, 1}, xi::Array{Float64,2}, J_opt_Int::Array{Int64,2}, J_opt_Float::Array{Float64,2})
	fout = open(fname_out, "w")
	n_max = size(J_opt_Int,1)
	J_out = zeros(q*L, q*L)	
	for n in 1:n_max
		i, j, a, b, = J_opt_Int[n, 1], J_opt_Int[n, 2] , J_opt_Int[n, 3], J_opt_Int[n, 4]
		J_elem, delta_J_elem, J_block, delta_J_block, MI_ele = J_opt_Float[n, 1], J_opt_Float[n, 2], J_opt_Float[n, 3], J_opt_Float[n, 4], J_opt_Float[n, 5] 
		println(fout, "J ", i-1, " ", j-1, " ", a-1, " ", b-1, " ", J_elem, " ", delta_J_elem, " ", J_block, " ", delta_J_block, " ", MI_ele)	
	end
	
	for i=1:L
		for mu=1:P 
			for a=1:q
				println(fout, "xi ", i-1, " ", mu-1, " ", a-1, " ", xi[mu, km(i,a,q)])
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

"""
	f1_f2(X::Array{Int64, 2}, W::Array{Float64, 2}, q::Int64)
	Compute frequencies of single and double sites, f1((i-1)*q+a), f2((i-1)*q+a, (j-1)*q+b).
	Output: f1::Array{Float64, 1} and f2::Array{Float64, 2}.
"""
function f1_f2(X::Array{Int64, 2}, W::Array{Float64, 1}, q::Int64)
	M,L = size(X)
	Meff = sum(W); scale = 1.0 /Meff
	f1 = zeros(Float64, L*q)
	f2 = zeros(Float64, L*q, L*q)
	for m = 1:M
		for i = 1:L
			a = X[m,i]+1	
			f1[km(i,a,q)] += W[m] * scale
			for j = (i+1):L
				b = X[m,j]+1
				f2[km(i,a,q), km(j,b,q)] += W[m] * scale
				f2[km(j,b,q), km(i,a,q)] += W[m] * scale
			end
		end
	end
	
	return (f1, f2) 
end


function f1_f2_onehot(X::Array{Int64, 2}, W::Array{Float64, 1}, q::Int64, L::Int64)
	M,L1 = size(X)
	Meff = sum(W); scale = 1.0 /Meff
	f1 = zeros(Float64, L*q)
	f2 = zeros(Float64, L*q, L*q)
	for m = 1:M
		for i = 1:L
			for a= 1:q
			f1[(i-1)*q+a] += W[m] * scale * X[m,km(i,a,q)]
			for j = (i+1):L
				for b = 1:q
				#b = X[m,j]+1
				f2[(i-1)*q+a, (j-1)*q+b] += W[m] * scale * X[m,km(i,a,q)]*X[m,km(j,b,q)]
				f2[(j-1)*q+b, (i-1)*q+a] += W[m] * scale * X[m,km(i,a,q)]*X[m,km(j,b,q)]
				end
			end
			end
		end
	end
	
	return (f1, f2) 
end



"""
	f3_c3(X::Array{Int64, 2}, W::Array{Float64, 2}, q::Int64, threshold::Float64)
	3pt frequencies and correlations.  
	Output: c3p::Array{corr3p, n_elem}, 
	where corr3p is struct that contains indices of positions 'i,j,k' and states 'a,b,c',  3pt-freq and 3pt-corr.  
"""

"""
    mutable struct corr3p
Three body correlation at positions `i, j, k` for states `a, b, c`, with value `cor`. 
"""
mutable struct corr3p
    i::Int64
    j::Int64
    k::Int64
    a::Int64
    b::Int64
    c::Int64
    frq::Float64
    cor::Float64
end


function f3_c3(X::Array{Int64, 2}, W::Array{Float64, 1}, q::Int64, threshold::Float64)
	M,L = size(X)
	Meff = sum(W); scale = 1.0/Meff
	@time (f1,f2) = f1_f2(X, W, q) 
	
	c3p = Array{corr3p,1}(undef,0)
	println("type of c3p = ", typeof(c3p) )	
	f3 = zeros(Float64, q*q*q) 
	for i in 1:L
		for j in (i+1):L
			for k in (j+1):L
				f3 .= 0.
				for m in 1:M
					a,b,c = X[m,i]+1, X[m,j]+1, X[m,k]+1
                    			f3[(a-1)*q*q + (b-1)*q + c] += W[m] * scale 
				end
				c3 = copy(f3)
				
				for a in 1:q
					for b in 1:q
						for c in 1:q
							c3[(a-1)*q*q + (b-1)*q + c] -= f2[(i-1)*q+a, (j-1)*q+b]*f1[(k-1)*q+c] - f2[(i-1)*q+a, (k-1)*q+c]*f1[(j-1)*q+b] - f2[(j-1)*q+b, (k-1)*q+c]*f1[(i-1)*q+a] + 2*f1[(i-1)*q+a]*f1[(j-1)*q+b]*f1[(k-1)*q+c]
							# Storing only if above threshold
							if abs(c3[(a-1)*q*q + (b-1)*q + c]) >= threshold
								push!(c3p,corr3p(i,j,k,a,b,c,f3[(a-1)*q*q + (b-1)*q + c],c3[(a-1)*q*q + (b-1)*q + c]))
							end
						end
					end
				end
			end
		end
	end
	return c3p
end

function f3_c3_two_set(X1::Array{Int64, 2}, X2::Array{Int64, 2}, q::Int64, L, threshold::Float64)
	@show M1,L1 = size(X1)
	@show M2,L2 = size(X2)
	scale1 = 1.0/M1
	scale2 = 1.0/M2
	@time (f1_1,f2_1) = f1_f2(X1, ones(M1), q) 
	@time (f1_2,f2_2) = f1_f2(X2, ones(M2), q) 
	
	c3p1 = []#Array{corr3p,1}(undef,0)
	c3p2 = []#Array{corr3p,1}(undef,0)
	
	#println("type of c3p = ", typeof(c3p1) )	
	f3_1 = zeros(Float64, q*q*q) 
	f3_2 = zeros(Float64, q*q*q) 
	for i in 1:L
		for j in (i+1):L
			for k in (j+1):L
				f3_1 .= 0.
				f3_2 .= 0.
				for m in 1:M1
					a,b,c = X1[m,i]+1, X1[m,j]+1, X1[m,k]+1
                    			f3_1[(a-1)*q*q + (b-1)*q + c] +=  scale1 
				end
				for m in 1:M2
					a,b,c = X2[m,i]+1, X2[m,j]+1, X2[m,k]+1
                    			f3_2[(a-1)*q*q + (b-1)*q + c] +=  scale2
				end
				c3_1 = copy(f3_1)
				c3_2 = copy(f3_2)
				
				for a in 1:q
					for b in 1:q
						for c in 1:q
							c3_1[(a-1)*q*q + (b-1)*q + c] -= f2_1[(i-1)*q+a, (j-1)*q+b]*f1_1[(k-1)*q+c] - f2_1[(i-1)*q+a, (k-1)*q+c]*f1_1[(j-1)*q+b] - f2_1[(j-1)*q+b, (k-1)*q+c]*f1_1[(i-1)*q+a] + 2*f1_1[(i-1)*q+a]*f1_1[(j-1)*q+b]*f1_1[(k-1)*q+c]
							c3_2[(a-1)*q*q + (b-1)*q + c] -= f2_2[(i-1)*q+a, (j-1)*q+b]*f1_2[(k-1)*q+c] - f2_2[(i-1)*q+a, (k-1)*q+c]*f1_2[(j-1)*q+b] - f2_2[(j-1)*q+b, (k-1)*q+c]*f1_2[(i-1)*q+a] + 2*f1_2[(i-1)*q+a]*f1_2[(j-1)*q+b]*f1_2[(k-1)*q+c]
							# Storing only if above threshold
							if abs(c3_1[(a-1)*q*q + (b-1)*q + c]) >= threshold
								push!(c3p1,c3_1[(a-1)*q*q + (b-1)*q + c])
								push!(c3p2,c3_2[(a-1)*q*q + (b-1)*q + c])
							end
						end
					end
				end
			end
		end
	end
	return (c3p1, c3p2) 
end

function f3_c3_two_set_onehot(X1::Array{Int64, 2}, X2::Array{Int64, 2}, q::Int64, L, threshold::Float64)
	@show M1,L1 = size(X1)
	@show M2,L2 = size(X2)
	scale1 = 1.0/M1
	scale2 = 1.0/M2
	@time (f1_1,f2_1) = f1_f2_onehot(X1, ones(M1), q, L) 
	@time (f1_2,f2_2) = f1_f2_onehot(X2, ones(M2), q, L) 
	
	c3p1 = Array{corr3p,1}(undef,0)
	c3p2 = Array{corr3p,1}(undef,0)
	println("type of c3p = ", typeof(c3p1) )	
	f3_1 = zeros(Float64, q*q*q) 
	f3_2 = zeros(Float64, q*q*q) 
	for i in 1:L
		for j in (i+1):L
			for k in (j+1):L
				f3_1 .= 0.
				f3_2 .= 0.
				for m in 1:M1
					for a in 1:q	
					for b in 1:q	
					for c in 1:q	
					#a,b,c = X1[m,i]+1, X1[m,j]+1, X1[m,k]+1
					f3_1[(a-1)*q*q + (b-1)*q + c] +=  scale1 * X1[m, km(i,a,q)] * X1[m, km(j,b,q)] * X1[m, km(k,c,q)]   
					end	
					end	
					end	
				end
				for m in 1:M2
					#a,b,c = X2[m,i]+1, X2[m,j]+1, X2[m,k]+1
                    			#f3_2[(a-1)*q*q + (b-1)*q + c] +=  scale2
					for a in 1:q	
					for b in 1:q	
					for c in 1:q	
					#a,b,c = X1[m,i]+1, X1[m,j]+1, X1[m,k]+1
					f3_2[(a-1)*q*q + (b-1)*q + c] +=  scale2 * X2[m, km(i,a,q)] * X2[m, km(j,b,q)] * X2[m, km(k,c,q)]   
					end	
					end	
					end	
				end
				c3_1 = copy(f3_1)
				c3_2 = copy(f3_2)
				
				for a in 1:q
					for b in 1:q
						for c in 1:q
							c3_1[(a-1)*q*q + (b-1)*q + c] -= f2_1[(i-1)*q+a, (j-1)*q+b]*f1_1[(k-1)*q+c] - f2_1[(i-1)*q+a, (k-1)*q+c]*f1_1[(j-1)*q+b] - f2_1[(j-1)*q+b, (k-1)*q+c]*f1_1[(i-1)*q+a] + 2*f1_1[(i-1)*q+a]*f1_1[(j-1)*q+b]*f1_1[(k-1)*q+c]
							c3_2[(a-1)*q*q + (b-1)*q + c] -= f2_2[(i-1)*q+a, (j-1)*q+b]*f1_2[(k-1)*q+c] - f2_2[(i-1)*q+a, (k-1)*q+c]*f1_2[(j-1)*q+b] - f2_2[(j-1)*q+b, (k-1)*q+c]*f1_2[(i-1)*q+a] + 2*f1_2[(i-1)*q+a]*f1_2[(j-1)*q+b]*f1_2[(k-1)*q+c]
							# Storing only if above threshold
							if abs(c3_1[(a-1)*q*q + (b-1)*q + c]) >= threshold
								push!(c3p1,corr3p(i,j,k,a,b,c,f3_1[(a-1)*q*q + (b-1)*q + c],c3_1[(a-1)*q*q + (b-1)*q + c]))
								push!(c3p2,corr3p(i,j,k,a,b,c,f3_2[(a-1)*q*q + (b-1)*q + c],c3_2[(a-1)*q*q + (b-1)*q + c]))
							end
						end
					end
				end
			end
		end
	end
	return (c3p1, c3p2) 
end

function save_ppv_fig(fname_key, p_list, typename_key, type_key, image_key)
    fnamekey_add1 =  "/PPV_0.0_"
    fnamekey_add2 =  "/PPV_1.0_"
    mycolor = cgrad([:yellow, :black],size(p_list, 1))
    for i in 1:size(type_key,1)
        m=1
        println(p_list[m], " ", type_key[i])
        if(i!=size(type_key,1))
            fnamekey_add=fnamekey_add1        end
        if(i==size(type_key,1))
            fnamekey_add=fnamekey_add2        end
        x_ppv = readdlm(fname_key*p_list[m]*fnamekey_add*typename_key[i])
        Plots.plot(x_ppv[:,end],color=mycolor[m], markersize=1, markeralpha=1,  xaxis=:log10, label="P="*p_list[m], legend=:bottomleft, legendfontsize=12, tickfontsize=12)
        for m in 2:size(p_list, 1)
            #println(p_list[m])
            x_ppv = readdlm(fname_key*p_list[m]*fnamekey_add*typename_key[i])
            Plots.plot!(x_ppv[:,end],color=mycolor[m], markersize=1, markeralpha=1,  xaxis=:log10, label="P="*p_list[m], legend=:bottomleft, legendfontsize=12, tickfontsize=12)
        end
        x_ppv = readdlm(fname_key*p_list[m]*fnamekey_add1*typename_key[1])
        Plots.plot!(x_ppv[:,end], color="red", markersize=1, markeralpha=1,  xaxis=:log10, label="MI", legend=:bottomleft, legendfontsize=12, tickfontsize=12)
        #println("test", type_key[i])
        Plots.savefig("PPV_"*type_key[i]*image_key*".png")
    end
end

function J_to_Frob(q,L,mat)
    #suppose mat is qL x qL matrix.
    myFrob=zeros(L,L)
    for i in 1:L
    for j in i:L
        temp=0.0
        for a in 1:q
            for b in 1:q
                temp += mat[km(i,a,q), km(j,b,q)]^2
            end
        end
        temp = sqrt(temp)
        myFrob[i,j] = temp
        myFrob[j,i] = temp
    end
    end
    return myFrob
end

# -------------- Usage ---------------#
#fname_in = "../note/MSA_artificial_q4_PF14.txt"
#X = readdlm(fname_in, Int64);
#X = X + 2*ones(Int, size(X));
#fname_out = "MSA_artificial_q4_PF14.fasta"
#num_to_fasta(X, fname_out)
function num_to_fasta(X, fname_out)
    code_alpha = ["-","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"];
    fout = open(fname_out, "w")
    N,L = size(X);
    for n in 1:N
        println(fout, ">Seq-ID-", string(n))
        for i in 1:L
            print(fout, code_alpha[X[n,i]])
        end
        println(fout, "")
    end
    close(fout)
end

# ---------- Skew Normal dist. --------#
# Ref. https://discourse.julialang.org/t/skew-normal-distribution/21549/3
#using StatsFuns
function rand_skewnormal(alpha)
    while true
        z = randn()
        u = rand()
        if u < StatsFuns.normcdf(z*alpha)
            return z
        end
    end
end

#P is dimension 
function get_projection(Mat, P)
    evl, evt = eigen(Mat);
    #Identity_mat = evt*evt'
    #c2_recnst = evt * Diagonal(evl) * evt';
    Project = evt[:, (end-P+1):end]';
    return Project
end


