using LinearAlgebra

# Newton-Raphson method with approximated Jacobian Matrix 
# using finite difference method for derivatives.
# f = function, x0 = initial guess, Na = Number of dimensions (of f and x)
# tol = tolerance (norm(F) < tol ), dx = difference

function newton_rapson(f,x0::Array{Float64} , tol::Float64,dx::Float64)
	
	Na = length(x0)
	J = ones(Na,Na) 			        # Jacobian Matrix 
	F = ones(Na)::Array{Float64} 		# Vector function 
	x1 = ones(Na)::Array{Float64} 		# New function 
	dFp = Array{Float64}(undef,Na, Na) 	# F(x0 + dx) Forward 
	dFn = Array{Float64}(undef,Na, Na)	# F(x0 - dx) Backward
	Id = 1*Matrix(I,Na,Na)
	it = 0 ::Int64
	
	for it in 0:1000
	
		F = f(x0) 		# evaluate and store F(x0) 
		
		for i in 1:Na
			
			dFp[i,:] = f(x0 + dx*Id[:,i]) # evaluate and store F(x0 + dx)
			dFn[i,:] = f(x0 - dx*Id[:,i]) # evaluate and store F(x0 - dx)
			
			for j in 1:Na
				J[i,j] = (dFp[i,j]-dFn[i,j])/(2*dx) # Building the Jacobian Matrix
			end
		end
		
		x1 = x0 - inv(J)*F 		      # Evaluating the x1 from x0
		
		if norm(F) < tol
			return [x1,it,norm(F),true]
				# [ answer , number_iterations , norm_of_F, convergence ]
		else
			x0 = x1
		end
	end
	return [x1,it,norm(F),false]
		# [ answer , number_iterations , norm_of_F, convergence ]
end
