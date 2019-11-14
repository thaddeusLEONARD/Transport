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
    J = Int[]
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
            if inst.nodes[i].vertex_type!="D"
                append!(Jnonp,[i])
            end
        end
        if inst.nodes[i].vertex_type=="S"
            append!(Js,[i])
        end
        if inst.nodes[i].vertex_type!="D"
            append!(J,[i])
        end
    end
    return Jl,Js,Jp,Jnonp,JnonpL,J
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

    Jl,Js,Jp , Jnonp,JnonpL, J= get_ensemble(inst)
    println(J)
    @variable(model, x[1:n,1:n] ,Bin);
    @variable(model, y[1:n,1:n] ,Bin);
    @variable(model, z[1:n,1:n] ,Bin);
    #@variable(model, 1>=x[1:n,1:n]>=0);
    #@variable(model, 1>=y[1:n,1:n]>=0);
    #@variable(model, 1>=z[1:n,1:n]>=0);
    @variable(model, qL[1:n]>=0 ); #
    @variable(model, qS[1:n]>=0 ); #
    @variable(model, Q>=Cs[1:n]>=0 ); #
    @variable(model, Re[1:n] ,Bin);
    # constraints -- satisfy the demand exactly
    #@constraint(model, Demand[j=1:nf], sum(m[i,j]*x[i] for i=1:n) == D[j]);

    # objective -- minimize total waste
    @objective(model, Min, sum(sum(x[i,j]*(inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx])+y[i,j]*(inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx]*inst.speed_ratio)-z[i,j]*(inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx]*inst.speed_ratio) for i=1:n) for j in 1:n));


    # contraintes de depot entrée sortie 1fois et flots
    #@constraint(model,[j = JnonpL], sum( x[i,j] for i in 1:n j!=i) <= 1);
    #@constraint(model,[i = JnonpL], sum( x[i,j] for j in 1:n j!=i) <= 1);
    @constraint(model, Cs[1]==Q)
    @constraint(model, qL[n]<=T)
    @constraint(model, qS[n]<=T)

    ## le x peut pas aller dans un Js
    @constraint(model,[i = 1:n,j = Js,i!=j], x[i,j] == 0);

    ## on peut pas aller sur soi meme
    @constraint(model,[i = 1:n], x[i,i] == 0);
    @constraint(model,[i = 1:n], y[i,i] == 0);

    # on sort un fois de 1 et on rentre un fois dans n
    @constraint(model, sum(x[1,j] for j in J) ==1 );
    @constraint(model, sum(y[1,j] for j in J) ==1 );
    @constraint(model, sum(x[j,n] for j in J) ==1 );
    @constraint(model, sum(y[j,n] for j in J) ==1 );
    @constraint(model, sum(z[1,j] for j in J) ==1 );
    @constraint(model, sum(z[j,n] for j in J) ==1 );

    # on rentre jamais en 1 et on sort jamais de n
    @constraint(model, sum(x[n,j] for j in 1+1:n) ==0 );
    @constraint(model, sum(y[n,j] for j in 1+1:n) ==0 );
    @constraint(model, sum(x[j,1] for j in 1:n-1) ==0 );
    @constraint(model, sum(y[j,1] for j in 1:n-1) ==0 );


    # pour tout i,j in J on rentre autant de fois dans un noeud qu'on en sort
    @constraint(model,[i = J], sum(x[i,j] for j in J) == sum(x[j,i] for j in J) );
    @constraint(model,[i = J], sum(y[i,j] for j in J) == sum(y[j,i] for j in J) );
    @constraint(model,[i = Jnonp], sum(z[i,j] for j in J) == sum(z[j,i] for j in J));


    @constraint(model,[i = 1:n], sum(x[i,j] for j in 1:n)<=1)
    @constraint(model,[i = 1:n], sum(y[i,j] for j in 1:n)<=1)
    @constraint(model,[j = 1:n], sum(x[i,j] for i in 1:n)<=1)
    @constraint(model,[j = 1:n], sum(y[i,j] for i in 1:n)<=1)


    ## tu rentres exactement une fois dans chaque Jnonp
    @constraint(model,[j = Jnonp], sum(x[i,j] for i in 1:n j!=i)+sum(y[i,j] for i in 1:n)-sum(z[i,j] for i in 1:n j!=i)   == 1);


    #def de z
    @constraint(model,[i = 1:n,j = 1:n,i!=j], z[i,j] <= (x[i,j]+y[i,j])/2);

    #def de rechargement
    @constraint(model,[i = Jp], Re[i] <= 1 - (qS[i]-qL[i])/T);

    ## demande de i vers j pour yij
    @constraint(model,[i=1:n,j = Jnonp], Cs[j] <= Cs[i]-inst.nodes[j].demand*(1-x[i,j])+Q*(1-y[i,j]));

    ## ravitaillement des yij
    @constraint(model,[i=1:n,j = Jp], Cs[j] >= Q*Re[j]);

    #time windows
    @constraint(model,[j = 1:n,i = 1:n,i!=j],  qL[j]+T*(1-x[i,j])>=qL[i]+inst.service_duration+inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx] );
    @constraint(model,[j = 1:n,i = 1:n,i!=j],  qS[j]+T*(1-y[i,j])>=qS[i]+inst.service_duration+inst.dist_matrix[inst.nodes[i].vertex_idx,inst.nodes[j].vertex_idx]);
    #=@constraint(model,[i = 1:n],  qL[i]>=inst.service_duration+inst.nodes[i].TW_min)
    @constraint(model,[i = 1:n],  qS[i]>=inst.service_duration+inst.nodes[i].TW_min)
    @constraint(model,[i = 1:n],  qL[i]<=inst.nodes[i].TW_max)
    @constraint(model,[i = 1:n],  qS[i]<=inst.nodes[i].TW_max)=#


    optimize!(model)
    x_val = JuMP.value.(x)
    for i=1:n
        for j=1:n
            if x_val[i,j]>0
                println(" i,j : ",i,",",j," val_x : ",x_val[i,j])
            end
        end

    end
    x_val = JuMP.value.(y)
    cs = JuMP.value.(Cs)
    for i=1:n
        for j=1:n
            if x_val[i,j]>0
                println(" i,j : ",i,",",j," val_y : ",x_val[i,j])
            end
        end
        println(" i : ",i," cs : ",cs[i])
    end
    x_val = JuMP.value.(z)

    for i=1:n
        for j=1:n
            if x_val[i,j]>0
                println(" i,j : ",i,",",j," val_z : ",x_val[i,j])
            end
        end

    end
    qs = JuMP.value.(qS)
    ql = JuMP.value.(qL)
    for i=1:n

        println(" i : ",i," qs : ",qs[i])
        println(" i : ",i," ql : ",ql[i])

    end

    return

end # end buildTSP




function main()
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
end
