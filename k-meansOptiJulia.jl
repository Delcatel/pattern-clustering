Pkg.add("Ipopt")
using JuMP
#using Clp
using Ipopt


N = 1000
A = randn(N,2)
sigma = 0.5
nbc = 4
H = randn(1,2)
norme = norm(H)
phi = acos(H[1]/norme)
rotation = [[cos(phi),-sin(phi)],[sin(phi),cos(phi)]]

for k =1:N
  a2k1 = A[2*k-1]
	A[2*k-1] = norme*dot(rotation[1],[A[2*k-1],A[2*k]]) + randn()
	A[2*k] = norme*dot(rotation[2],[a2k1,A[2*k]]) + randn()
end

m = Model(solver=IpoptSolver())

@variable(m, theta)
@variable(m, r >= 0 )
@variable(m, C[1:2*nbc] )
@variable(m, T[1:N*nbc] >= 0)
s = 0
#for i = 1:N
	#for j=1:nbc
		#s = s + T[(i-1)*nbc + j]*( (A[2*i - 1] - C[2*j-1])^2 + (A[2*i] - C[2*j] )^2 )
	#end
#end

@setNLObjective(m, Min,sum(sum(T[(i-1)*nbc + j]*( (A[2*i - 1] - C[2*j-1])^2 + (A[2*i] - C[2*j] )^2 ) for j=1:nbc) for i=1:N))
#rot = [[cos(theta),-sin(theta)],[sin(theta),cos(theta)]]
#@NLobjective(m, Min, 0.5*s)

#for i=1:N
	#@constraint(m, T[i*nbc +1] + T[i*nbc + 2] + T[i*nbc + 3] + T[i*nbc + 4] = 1)
#end
@NLconstraint(m, C[3] == C[1] +  cos(theta)*r)
@NLconstraint(m, C[4] == C[2] +  sin(theta)*r)

@NLconstraint(m, C[5] == C[3] +  (-sin(theta))*r)
@NLconstraint(m, C[6] == C[4] +  cos(theta)*r)

@NLconstraint(m, C[7] == C[1] +  (-sin(theta))*r)
@NLconstraint(m, C[8] == C[2] +  cos(theta)*r)
#########
#@constraint(m, C[3] == C[1] +  dot(rot[1],[r,0]))
#@constraint(m, C[4] == C[2] +  dot(rot[2],[r,0]))

#@constraint(m, C[5] == C[3] +  dot(rot[1],[0,r]))
#@constraint(m, C[6] == C[4] +  dot(rot[2],[0,r]))

#@constraint(m, C[7] == C[1] +  dot(rot[1],[0,r]))
#@constraint(m, C[8] == C[2] +  dot(rot[2],[0,r]))

status = solve(m)
println(getValue(r))
println(getValue(theta))
println(getValue(C[1]))
println(getValue(C[2]))
