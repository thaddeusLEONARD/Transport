# Using CBC https://github.com/JuliaOpt/Cbc.jl
using Cbc
# Using Gurobi
using JuMP,JSON
#using GLPKMathProgInterface
#using GLPK
include("instance.jl")

function get_ensemble(inst::Instance)
    n = size(inst.nodes)[1]
    Jl = Int[]
    Jp = Int[]
    Js = Int[]
    Jnonp = Int[]
    JnonpL = Int[]
    for i in 1:n
        if inst.nodes[i].vertex_type!="S"
            append!(Jl,[i])
            if inst.nodes[i].vertex_type!="P"
                append!(JnonpL,[i])
            end
        end
        if inst.nodes[i].vertex_type=="P"
            append!(Jp,[i])
        else
            append!(Jnonp,[i])
        end
        if inst.nodes[i].vertex_type=="S"
            append!(Js,[i])
        end
    end
    return Jl,Js,Jp,Jnonp,JnonpL
end

function get_max(inst::Instance)
    n = size(inst.nodes)[1]
    tmp_max = 0
    for i in 1:n
            if inst.nodes[i].TW_min+inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[n].vertex_idx]>tmp_max
                tmp_max=inst.nodes[i].TW_min+inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[n].vertex_idx]
            end
    end
    return tmp_max
end


# buildTSP
# Given a matrix of city locations, build the TSP
# Inputs:
#   n       Number of cities
#   durations  n-by-n matrix of travel times between cities
# Output:
#   m       JuMP model

function buildTSP(inst::Instance)
    model = Model(with_optimizer(Cbc.Optimizer, logLevel=1))
    #model = Model(with_optimizer(Gurobi.Optimizer, Presolve=0, OutputFlag=0))
    n = size(inst.nodes)[1]
    T = inst.time_horizon
    tmp_max =get_max(inst)
    Q = inst.SV_cap

    Jl,Js,Jp , Jnonp,JnonpL= get_ensemble(inst)

    @variable(model, x[1:n,1:n] ,Bin);
    @variable(model, y[1:n,1:n] ,Bin);
    @variable(model, z[1:n,1:n] ,Bin);
    @variable(model, qL[1:n]>=0 ); #
    @variable(model, qS[1:n]>=0 ); #
    @variable(model, Q>=Cs[1:n]>=0 ); #
    @variable(model, Re[1:n] ,Bin);
    # constraints -- satisfy the demand exactly
    #@constraint(model, Demand[j=1:nf], sum(m[i,j]*x[i] for i=1:n) == D[j]);

    # objective -- minimize total waste
    @objective(model, Min, sum(sum(x[i,j]*(inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx])+y[i,j]*(inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx]*inst.speed_ratio)-z[i,j]*(inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx]*inst.speed_ratio) for i=1:n) for j in 1:n));


    # contraintes de depot entrée sortie 1fois
    #@constraint(model,[j = JnonpL], sum( x[i,j] for i in 1:n j!=i) <= 1);
    #@constraint(model,[i = JnonpL], sum( x[i,j] for j in 1:n j!=i) <= 1);
    @constraint(model,[i = 1:n,j = Js,i!=j], x[i,j] == 0);
    @constraint(model,[j = 1:n,i = Js,i!=j], x[i,j] == 0);
    @constraint(model,[i = 1:n], x[i,i] == 0);
    @constraint(model,[i = 1:n], y[i,i] == 0);
    @constraint(model, x[n,1] == 1);
    @constraint(model, y[n,1] == 1);
    @constraint(model, x[1,n] == 0);
    @constraint(model, y[1,n] == 0);
    @constraint(model,[i = 1:n], sum(x[i,j] for j in 1:n) == sum(x[j,i] for j in 1:n) );
    @constraint(model,[i = 1:n], sum(y[i,j] for j in 1:n) == sum(y[j,i] for j in 1:n) );
    @constraint(model,[i = Jnonp], sum(z[i,j] for j in 1:n) == sum(z[j,i] for j in 1:n));
    #@constraint(model,[i = Jl], sum(x[i,j] for j in 1:n)<=1)
    #@constraint(model,[i = 1:n], sum(y[i,j] for j in 1:n)<=1)

    #@constraint(model,[j = Jnonp], sum(x[i,j] for i in 1:n)+sum(y[i,j] for i in 1:n)-sum(z[i,j] for i in 1:n)   == 1);

    @constraint(model,[j = Js], sum( y[i,j] for i in 1:n j!=i) == 1);

    #@constraint(model,[i = Js], sum( y[i,j] for j in 1:n j!=i) == 1);

    @constraint(model,[i = 1:n,j = 1:n,i!=j], z[i,j] <= (x[i,j]+y[i,j])/2);

    @constraint(model,[i = 1:n], Re[i] <= 1 - (qL[i]-qS[i])/T);
    ## demande de i vers j pour yij
    @constraint(model,[i=1:n,j = 1:n], Cs[j] <= Cs[i]-inst.nodes[j].demand+Q*(1-y[i,j]));
    ## ravitaillement des yij
    @constraint(model,[i=1:n,j = Jp], Cs[j] >= Q*Re[j]);
    @constraint(model,[j = 1:n],  qL[j]>=inst.nodes[j].TW_min );
    @constraint(model,[j = n],  qL[j]<=tmp_max );

    @constraint(model,[j = 1:n,i = 1:n,i!=j],  qL[j]+T*(1-x[i,j])>=qL[i]+inst.service_duration+inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx] );

    @constraint(model,[j = 1:n,i = 1:n,i!=j],  qS[j]+T*(1-y[i,j])>=qS[i]+inst.service_duration+inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx] );

    optimize!(model)

    return cyclex,cycley

end # end buildTSP





# Main program starting here

# PARSE iNstance, parameters and distanceMatrix
instdir = "..\\instances2019\\"
outdir = "..\\output\\"
instNames= ["C1-2-8.txt","C2-2-8.txt","C1-3-10.txt","C1-3-12.txt","C2-3-10.txt","C2-3-12.txt",
      "R1-2-8.txt","R2-2-8.txt","R1-3-10.txt","R1-3-12.txt","R2-3-10.txt","R2-3-12.txt"]

paramf = string(instdir,"parameters.txt")
instf =string(instdir,instNames[1])
matf =string(instdir,"distancematrix98.txt")

inst = parseInstance(paramf,instf,matf)

@time  cyclex,cycley = buildTSP(inst)
println(cyclex)
println(cycley)
# saving some results to be loaded somewhere else:
#saving in a CSV File

CSV.write(string(outdir,"routes.csv"), inst.nodes[cyclex],delim='\t')

# saving in JSON
open(string(outdir,"routes.json"),"w") do f
    JSON.print(f, inst.nodes[cyclex], 4)
end
