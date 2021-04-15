using LinearAlgebra

function newton_rapson(f,x0,Na)#::Array{Float64}
	# f = função, x0 = chute, N-dim
	J = ones(Na,Na) # construindo a Jacobiana
	F = ones(Na)::Array{Float64}
	x1 = ones(Na)::Array{Float64}
	dFp = Array{Float64}(undef,Na, Na)# vai guardar as derivadas de F
	#com relação a xi
	dFn = Array{Float64}(undef,Na, Na)
	dx = 0.01
	Id = 1*Matrix(I,Na,Na)#::Array{Float64}
	it = 0 ::Int64
	println("Teste")
	#println([x1,it,norm(F)])
	for it in 0:1000
		F = f(x0) # F no passo atual
		println(x0 + dx*Id[:,1])
		for i in 1:Na
			#println(x0 + dx*Id[:,i])
			dFp[i,:] = f(x0 + dx*Id[:,i])
			dFn[i,:] = f(x0 - dx*Id[:,i])
			#println(dF)
			for j in 1:Na
				J[i,j] = (dFp[i,j]-dFn[i,j])/(2*dx)
			end
		end
		println(J)
		x1 = x0 - inv(J)*F
		
		#println([x1,it,norm(F)])
		if norm(F) < 1.0e-10
			return [x1,it,norm(F)]
		else
			x0 = x1
		end
	end
	return [!true,it,norm(F)]
end
