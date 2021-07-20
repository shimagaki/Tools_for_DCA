
function get_Gamma(f2_original, f1_original, q, L, pseudo_para=0.05)
    f1 = (1. - pseudo_para) * f1_original + pseudo_para / q * ones(q*L) 
    f2 = (1. - pseudo_para) * f2_original + pseudo_para / q/q * ones(q*L, q*L);

    for i in 1:L
        for a in 1:q
            for b in 1:q
                f2[km(i,a,q), km(i,b,q)] = (1.0-pseudo_para) * f2_original[km(i,a,q), km(i,b,q)] + pseudo_para / q * kr(a,b)
            end
        end
    end

    
    Gamma = zeros(q*L, q*L)
    for i in 1:L
        for j in 1:L
            for a in 1:q 
                for b in 1:q 
                    temp = f2[ km(i,a,q), km(j,b ,q) ] - f1[km(i,a,q)] * f1[ km(j,b,q)]
                    Gamma[km(i,a,q), km(j,b,q)] = temp / sqrt(f1[km(i,a,q)] * f1[ km(j,b,q)])
                end
            end
        end
    end
    
    return Gamma
end

""" INPUT: (f2_original, f1_original, q, L, pseudo_para=0.05)
    OUTPUT: (rank, eigenvalues, eigenvectors) # taht is for Gamma matrix.
"""
function get_eigens_of_Gamma(f2_original, f1_original, q, L, pseudo_para=0.05)
    f2 = (1. - pseudo_para) * f2_original + pseudo_para / q/q * ones(q*L, q*L);
    f1 = (1. - pseudo_para) * f1_original + pseudo_para / q * ones(q*L) 
    

    for i in 1:L
        for a in 1:q
            for b in 1:q
                f2[km(i,a,q), km(i,b,q)] = (1.0-pseudo_para) * f2_original[km(i,a,q), km(i,b,q)] + pseudo_para / q * kr(a,b)
            end
        end
    end

    
    Gamma = zeros(q*L, q*L)
    
    for i in 1:L
        for j in 1:L
            for a in 1:q 
                for b in 1:q 
                    temp = f2[ km(i,a,q), km(j,b,q) ] - f1[km(i,a,q)] * f1[ km(j,b,q)]
                    #C[km(i,a,q), km(j,b,q)] = temp
                    Gamma[km(i,a,q), km(j,b,q)] = temp / sqrt(f1[km(i,a,q)] * f1[ km(j,b,q)])
                end
            end
        end
    end
    rank_Gamma = rank(Gamma)
    (evl, evt) = eigen(Gamma)
    
    ll_index = (q*L-rank_Gamma+1):q*L
    ll = evl[ll_index] - ones( size(evl[ll_index]) ) - log.(evl[ll_index]);
    
    
    return (rank_Gamma, ll, evl, evt)
end

""" INPUT: (q, L, P_m, P_p, rank_Gamma, evl, evt, f1), 
    
    OUTPUT: (Xi_p, Xi_m)
"""
function get_Xi_HP_model(q, L, P_p, P_m, rank_Gamma, evl, evt, f1)
    Xi_m = zeros(P_m, q*L)
    scale_m = zeros(P_m)
    for m in 1:P_m
        for i in 1:L
            for a in 1:q
		    scale_m[m] = sqrt(1.0/evl[end-rank_Gamma+m] - 1)
		    Xi_m[m, km(i,a,q)] = sqrt(1.0/evl[end-rank_Gamma+m] - 1) * evt[km(i,a,q), end-rank_Gamma+m] / sqrt(f1[km(i,a,q)]) 
            end
        end
    end
    
    Xi_p = zeros(P_p, q*L);
    scale_p = zeros(P_p)
    for m in 1:P_p
        for i in 1:L
            for a in 1:q
		    scale_p[m] = sqrt(1- 1.0/evl[end+1-m])
		    Xi_p[m, km(i,a,q)] = sqrt(1- 1.0/evl[end+1-m]) * evt[km(i,a,q),end+1-m] / sqrt(f1[km(i,a,q)]) 
            end
        end
    end
    return (Xi_p, Xi_m, scale_p, scale_m)
end

function convert_Xi_HP(Xi_p, Xi_m, P_p, P_m)
	Xi_p_conv = zeros(P_p, q*L)
	Xi_m_conv = zeros(P_m, q*L)
	for m in 1:P_p
	    for i in 1:L
		for a in 1:q
			Xi_p_conv[m, km(i,a,q)] = Xi_p[m, km(i,a,q)]
		end
	    end
	end

	for m in 1:P_m
	    for i in 1:L
		for a in 1:q
			Xi_m_conv[m, km(i,a,q)] = Xi_m[m, km(i,a,q)] 
		end
	    end
	end
	return Xi_p_conv, Xi_m_conv
end

function main_HP(q::Int64, L::Int64, fname_seq::String, file_key::String, fname_fig::String, P_p::Int64, P_m::Int64, th_reading_MSA=1.1, pseudo_para=0.01)
	#q,L = 21, 70
	#fname_seq = "/data/shimagaki/sparse-BM-Analysis/MSA/PF00076_mgap6_compresed_Weff_eq_W.dat"; 
	X = get_MSA(fname_seq,q, L, th_reading_MSA);
	(f1_original,f2_original,c2,Frob) = basic_MSA_analysis(X, q,L);
	(rank_Gamma, ll, eigenvalues, eigenvectors) = get_eigens_of_Gamma(f2_original, f1_original, q, L, pseudo_para);

	#P_p = 400; P_m = 400
	#pseudo_para=0.01
	(Xi_p, Xi_m, scale_p, scale_m)  = get_Xi_HP_model(q, L, P_p, P_m, rank_Gamma, eigenvalues, eigenvectors, f1_original*(1-pseudo_para) +pseudo_para*ones(size(f1_original)));

	f_name_out_para = "./para_xi_Pp"*string(P_p)*"_Pm"*string(P_m)*"_"*file_key*".dat"
	fout = open(f_name_out_para, "w")

	for m in 1:P_p
	    for i in 1:L
		for a in 1:q
		    println(fout, "Xi_p ", m-1, " ", i-1, " ", a-1, " ", Xi_p[m, km(i,a,q)])
		end
	    end
	end
	for m in 1:P_m
	    for i in 1:L
		for a in 1:q
		    println(fout, "Xi_m ", m-1, " ", i-1, " ", a-1, " ", Xi_m[m, km(i,a,q)])
		end
	    end
	end
	close(fout)

	p1 =Plots.plot(eigenvalues, markershape = :circle, markersize=5, markeralpha=0.5, color="blue", xaxis=:log, label="eigenvalues", legend=:topleft, legendfontsize=12)

	ll_index = (q*L-rank_Gamma+1):q*L
	p2 = Plots.plot(eigenvalues[ll_index], ll , markershape = :circle, markersize=5, markeralpha=0.5, color="blue", xaxis=:log, label="lgo-likelihood", xlabel="eingenvalues", legend=:topleft, legendfontsize=12)
	Plots.plot(p1,p2, layout = (2,1), size=(800,700) )
	Plots.savefig(fname_fig)
	return (Xi_p, Xi_m, scale_p, scale_m, p1, p2)
end

function main_HP_with_weight(q::Int64, L::Int64, fname_seq::String, fname_w::String, file_key::String, fname_fig::String, P_p::Int64, P_m::Int64, th_reading_MSA=1.1, pseudo_para=0.01)
	#q,L = 21, 70
	#fname_seq = "/data/shimagaki/sparse-BM-Analysis/MSA/PF00076_mgap6_compresed_Weff_eq_W.dat"; 
	X = get_MSA_weight(fname_seq, fname_w, q, L, th_reading_MSA);
	@show size(X)
	(f1_original,f2_original,c2,Frob) = basic_MSA_analysis(X, q,L);

	(rank_Gamma, ll, eigenvalues, eigenvectors) = get_eigens_of_Gamma(f2_original, f1_original, q, L, pseudo_para);

	#P_p = 400; P_m = 400
	#pseudo_para=0.01
	(Xi_p, Xi_m)  = get_Xi_HP_model(q, L, P_p, P_m, rank_Gamma, eigenvalues, eigenvectors, f1_original*(1-pseudo_para) +pseudo_para*ones(size(f1_original)));

	f_name_out_para = "./para_xi_Pp"*string(P_p)*"_Pm"*string(P_m)*"_"*file_key*".dat"
	fout = open(f_name_out_para, "w")

	for m in 1:P_p
	    for i in 1:L
		for a in 1:q
		    println(fout, "Xi_p ", m-1, " ", i-1, " ", a-1, " ", Xi_p[m, km(i,a,q)])
		end
	    end
	end
	for m in 1:P_m
	    for i in 1:L
		for a in 1:q
		    println(fout, "Xi_m ", m-1, " ", i-1, " ", a-1, " ", Xi_m[m, km(i,a,q)])
		end
	    end
	end
	close(fout)

	p1 =Plots.plot(eigenvalues, markershape = :circle, markersize=5, markeralpha=0.5, color="blue", xaxis=:log, label="eigenvalues", legend=:topleft, legendfontsize=12)

	ll_index = (q*L-rank_Gamma+1):q*L
	p2 = Plots.plot(eigenvalues[ll_index], ll , markershape = :circle, markersize=5, markeralpha=0.5, color="blue", xaxis=:log, label="lgo-likelihood", xlabel="eingenvalues", legend=:topleft, legendfontsize=12)
	Plots.plot(p1,p2, layout = (2,1), size=(800,700) )
	Plots.savefig(fname_fig)
	return (Xi_p, Xi_m, p1, p2)
end

