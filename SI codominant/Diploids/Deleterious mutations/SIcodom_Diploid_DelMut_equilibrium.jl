using Distributions
using StatsBase


## MUTATIONS

#Mutations on the genome 
function mutation_genome!(Npop, list_ind, U)
    for i in 1:Npop
        for j in 2:3
            nbMut = rand(Poisson(U)) #number of mutations
            posMut = rand(nbMut)  #positions of the mutations (creates a vector of length nbMut with random values in [0,1))
            chr = list_ind[i][j]
            if nbMut > 0
                append!(chr,posMut)
                sort!(chr)
            end
        end
    end
end


#Mutations Si to Sj 
function mutation_SI!(Npop, list_ind, U_SI, k)
    nbMut = rand(Binomial(2*Npop,U_SI)) #number of mutations among all S alleles
    for _ in 1:nbMut
        indMut = rand(1:Npop) #choosing the individual
        #choosing the allele of that individual that will mutate
        chrMut = rand(1:2)
        allMut = list_ind[indMut][1][chrMut] #S allele that will mutate
        if allMut != -1 #mutates only if the allele is SI
            mutant_S = rand(0:k-1) #Choosing a random S allele to replace the existing one
            while mutant_S == allMut #Making sure that the new S allele is different from the old one
                mutant_S = rand(0:k-1)
            end
            list_ind[indMut][1][chrMut] = mutant_S #Replacing with mutant S allele
        end
    end
end


#Mutations Si to SC
function mutation_SC!(Npop, list_ind, U_SC)
    nbMut = rand(Binomial(2*Npop,U_SC)) #number of mutations among all S alleles
    for _ in 1:nbMut
        indMut = rand(1:Npop) #choosing the individual
        #choosing the allele of that individual that will mutate
        chrMut = rand(1:2)
        allMut = list_ind[indMut][1][chrMut]
        if allMut != -1 #mutates only if the allele is SI
            list_ind[indMut][1][chrMut] = -1 #Replacing with mutant SC allele
        end
    end
end
            








## RECOMBINATION 

function recombination(par, L)
    #number and positions of the crossover
    nb_crossovers = rand(Poisson(L))
    pos_crossovers = rand(nb_crossovers)
    sort!(pos_crossovers)

    #recombined chromosomes
    chr1 = par[2]
    chr2 = par[3]
    Sallele1 = par[1][1]
    Sallele2 = par[1][2]
    #shuffling the two chromosomes
    if rand() < 0.5
        chr1,chr2 = chr2,chr1
        Sallele1,Sallele2 = Sallele2,Sallele1
    end

    gamete::Vector{Float64} = []
    sizehint!(gamete, floor((length(chr1)+length(chr2))/2)) #Suggest that gamete reserve some capacity in advance. This can improve performance.
    Sallele = Sallele1

    first_possible_index_chr1 = 1
    first_possible_index_chr2 = 1
    
    for i in 1:nb_crossovers + 1
        #Looking between two crossovers
        if nb_crossovers == 0
            j=0
            k=1
        elseif i == 1
            j = 0
            k = pos_crossovers[i]
        elseif i == nb_crossovers +1
            j = pos_crossovers[i-1]
            k = 1
        else
            j = pos_crossovers[i-1]
            k = pos_crossovers[i]
        end

        #Take on the first chromosome
        if i % 2 == 0
            if j <= 0.5 && k >=0.5
                Sallele = Sallele1
            end
            for l in first_possible_index_chr1:length(chr1)
                if chr1[l] >= k
                    first_possible_index_chr1 = l
                    break
                end
                if chr1[l] >= j
                    push!(gamete,chr1[l])
                end
            end
        end

        #Take on the second chromosome
        if i%2 != 0
            if j <= 0.5 && k >=0.5
                Sallele = Sallele2
            end
            for l in first_possible_index_chr2:length(chr2)
                if chr2[l] >= k
                    first_possible_index_chr2 = l
                    break
                end
                if chr2[l] >= j
                    push!(gamete,chr2[l])
                end
            end
        end
        
    end

    return (Sallele, gamete)
end








## FITNESS 

#Wi = (1 − hs)**Nhe (1 − s)**Nho, where Nhe and Nho are the number of mutations in the heterozygous and homozygous state within its genome


#returns the fitness of one individual
function fitness_ind(ind, h, s)
    #chromosomes
    c1 = ind[2]
    c2 = ind[3]
    s1 = length(c1)
    s2 = length(c2)

  #Couting the number of heterozygous and homozygous mutations
    Nho = 0

    i = 1
    j = 1
    while i <= s1 && j <= s2
        if c1[i] == c2[j]
            Nho += 1
            i += 1
        elseif c1[i] < c2[j]
            i += 1
        else
            j += 1
        end
    end

    Nhe = s1+s2-2*Nho

    W = (1-h*s)^Nhe*(1-s)^Nho #Computing the fitness of the individual
    return W
end


#returns the list of fitness for all individuals
function set_fitness!(list_ind, list_W, h, s)
    for i in 1:length(list_ind)
        list_W[i] = fitness_ind(list_ind[i],h,s)
    end
end







## SELFING RATE (OK)

#SC dominant in the pollen
function selfing_rate(Npop, list_ind, par, alpha, list_W, dominance)
    Salleles = list_ind[par][1] #list of the parent's S-alleles

    SR = 0.0
    #Special cases
    if Salleles[1] != -1 && Salleles[2] != -1 #the individual is SI
        return 0.0
    elseif Salleles[1] == -1 && Salleles[2] == -1 ##Two SC alleles, compatible with everyone
        gamma_ii = 2
        gamma_ij = 2
        sum_val = gamma_ij*(sum(list_W)-list_W[par])
        SR = (alpha*list_W[par]*gamma_ii)/(alpha*list_W[par]*gamma_ii + ((1-alpha)/(Npop-1))*sum_val)
    
        
    #Case 1 SC
    else

        #SC DOMINANT
        if dominance ==1
            gamma_ii = 2
    
            Si = Salleles[1]
            if Salleles[1] == -1
                Si = Salleles[2]
            end
    
            sum_val = 0
            for j in 1:Npop
                gamma_ij = 0
                if j != par
                    Sj = list_ind[j][1]
                    if Sj[1] == -1 || Sj[2] == -1 #At least 1SC allele in the pollen, compatible because SC dominant in the pollen
                        gamma_ij = 2
                    else
                        if Sj[1] == Si || Sj[2] == Si #At least 1 Sallele is shared with the mother, both alleles are SI cause SSI
                            gamma_ij = 0
                        else # No S alleles shared with the mother, both alleles are compatible
                            gamma_ij = 2
                        end
                    end
                end
                sum_val += list_W[j]*gamma_ij
            end      
    
            SR = (alpha*list_W[par]*gamma_ii)/(alpha*list_W[par]*gamma_ii + ((1-alpha)/(Npop-1))*sum_val)

        elseif dominance == 0
            return 0.0      
        end
        
    end

    return SR
end






## REPRODUCTION 

function reproduction!(Npop, list_ind, list_W, list_off, list_W_off, h,s, alpha, L, phase, dominance)
    fitnessMax = maximum(list_W)

    Npop_off = 0

    while Npop_off != Npop

        # Selecting a mother
        par1 = rand(1:Npop) #Selecting a parent randomly
        fitnessPar1 = list_W[par1]

        while fitnessPar1/fitnessMax < rand() #A new maternal parent is selected if the fitness of the first selected parent is too low
            par1 = rand(1:Npop)
            fitnessPar1 = list_W[par1]
        end
        

        Salleles_par1 = list_ind[par1][1]


        #Selfing or outcrossing ?
        selfRate = selfing_rate(Npop, list_ind, par1, alpha, list_W, dominance)

        #OUTCROSSING
        if selfRate <=  rand()
            OK_father = false
            while OK_father == false
                #Selecting a father
                par2 = rand(1:Npop) #A second parent is selected
                fitnessPar2 = list_W[par2]
                Salleles_par2 = list_ind[par2][1]
                
                #Checking the father
                if par2 == par1 # Checking that the father is different from the mother
                    continue
                end
                if fitnessPar2/fitnessMax < rand() # Checking its fitness
                    continue
                end
                
                # Checking the compatibility at the S-locus
                
                #SC DOMINANT
                if dominance ==1
                    if Salleles_par2[1] == -1 || Salleles_par2[2] == -1 # At least one allele in the pollen is SC
                        OK_father = true
                        break
                    end
                    if Salleles_par1[1] == Salleles_par2[1] || Salleles_par1[2] == Salleles_par2[1] # OK first S allele of the pollen
                        continue
                    end
                    if Salleles_par1[1] != Salleles_par2[2] && Salleles_par1[2] != Salleles_par2[2] #OK second S allele of the pollen
                        OK_father = true
                    end

                #SC RECESSIVE OR CODOMINANT
                elseif dominance == 0
                    if Salleles_par2[1] == -1 && Salleles_par2[2] == -1 # Both alleles in the pollen is SC
                        OK_father = true
                        break
                    else #1 or 0 SC alleles
                        if Salleles_par2[1] == -1 #One SC allele
                            if Salleles_par2[2] != Salleles_par1[1] && Salleles_par2[2] != Salleles_par1[2] #Check if the second is compatible
                                OK_father = true
                                break
                            end
                        elseif Salleles_par2[2] == -1
                            if Salleles_par2[1] != Salleles_par1[1] && Salleles_par2[1] != Salleles_par1[2] #Check if the second is compatible
                                OK_father = true
                                break
                            end
                        else #2 SI alleles
                            if Salleles_par1[1] != Salleles_par2[1] && Salleles_par1[1] != Salleles_par2[2] && Salleles_par1[2] != Salleles_par2[1] && Salleles_par1[2] != Salleles_par2[2]  
                                OK_father = true
                                break
                            end
                        end
                    end
                end
                    
            end
            
        #SELFING
        #All pollen grains are compatible when SC
        else
            par2 = par1
        end

        #Recombination
        Sallele1, chr1 = recombination(list_ind[par1], L)
        Sallele2, chr2 = recombination(list_ind[par2], L)
        #Creating a new individual for the next generation
        Npop_off += 1
        offspring = [[Sallele1, Sallele2], chr1, chr2]
        list_off[Npop_off] = offspring
        list_W_off[Npop_off] = fitness_ind(offspring,h, s)
    end
end












#INBREDING DEPRESSION

function inbreedingDepression(Npop, list_ind, list_W, h,s, L)
    #Sampling 100 individuals from the population
    nbSample = 100
    if nbSample > Npop
        nbSample = Npop
    end
    
    fitness_allof = 0
    fitness_autof = 0

    #outcrossing
    for _ in 1:nbSample
        par1 = rand(1:Npop)
        par2 = rand(1:Npop)
        while par1 == par2
            par2 = rand(1:Npop)
        end

        #Creating a new offspring
        Sallele1, chr1 = recombination(list_ind[par1], L)
        Sallele2, chr2 = recombination(list_ind[par2], L)
        offspring = [[Sallele1, Sallele2], chr1, chr2]

        fitness_allof += fitness_ind(offspring, h,s)
    end

    #selfing
    for _ in 1:nbSample
        par1 = rand(1:Npop)
        par2 = par1

        #Creating a new offspring
        Sallele1, chr1 = recombination(list_ind[par1], L)
        Sallele2, chr2 = recombination(list_ind[par2], L)
        offspring = [[Sallele1, Sallele2], chr1, chr2]

        fitness_autof += fitness_ind(offspring, h,s)
    end

    inbreeding_depression = 1 - fitness_autof/fitness_allof
    return inbreeding_depression
end






function dichotomie(U,Umin,Umax, inDep, delta)
    diff = inDep-delta
    if diff > 0 
        Umax = U
    else
        Umin = U
    end
    U = (Umin+Umax)/2

    return U, Umin, Umax
end



## MAIN 

#Phases:
# - Only mutations Si to Sj (2000 generations)
# - Then mutations on the genome (2000 generations)
# - Then mutations Si to SC (max 250000 generations, and stop when invasion of SC)

function simulation(nbSim, Npop, h,s, alpha, U, U_SI, U_SC, L, dominance, delta)
    #Creating output lists (for all runs)
    Tot_nbSalleles::Vector{Float64} = []
    Tot_ID::Vector{Float64} = []
    Tot_nbSC::Vector{Float64} = []
    Tot_freqSC::Vector{Float64} = []
    Tot_U::Vector{Float64} = []
    
    threshold_mutation_filter = Npop*100000.0

    Umin = 0
    Umax = 6
    U = (Umin+Umax)/2
    
    inDep = 2
    phase2_not_launched = true
    
    while phase2_not_launched == true    
        nbSalleles = 0
        nbSCalleles = 0
        freqSC = 0.0
        ID = 0
        
        it = 0
        println("U test=" , U)
        # Creating intermediate lists
        list_ID::Vector{Float64}, list_nbSalleles::Vector{Int} = [],[]
        list_nbSC::Vector{Int} = []


        #Initialisation
        nbS = 100 #number of possible S alleles
        list_ind::Vector{Vector{Vector{Float64}}} = []
        for i in 1:Npop
            ind::Vector{Vector{Float64}} = [[0.0,0.0], [],[]]
            push!(list_ind,ind)
            list_ind[i][1][1] = rand(0:nbS-1)
            list_ind[i][1][2] = rand(0:nbS-1) 
        end
        list_W = Vector{Float64}(undef, Npop)
        set_fitness!(list_ind, list_W, h, s)
        
        list_off::Vector{Vector{Vector{Float64}}} = []
        for _ in 1:Npop
            ind::Vector{Vector{Float64}} = [[0.0,0.0], [],[]]
            push!(list_off,ind)
        end
        list_W_off = Vector{Float64}(undef, Npop)


        generation = 0
        phase = 0
        
        max_gen = 500000

        #Phase 0: Starting with mutations on the S locus only (SI to SI)
        #Phase 1: Introduction of deleterious mutations on the genome
        #Phase 2: Introduction of S locus mutations to SC alleles
        
        while phase <= 2 
            generation += 1
            if generation%100==0
                println("generation=", generation)
            end


            #Creating the next generation
            reproduction!(Npop, list_ind, list_W, list_off, list_W_off, h,s, alpha, L, phase, dominance)
            
            #Mutations
            mutation_SI!(Npop, list_off, U_SI, nbS) #on the S locus: SI to SI
            if phase >= 1
                mutation_genome!(Npop, list_off, U) #on the genome
            end
            if phase >= 2
                mutation_SC!(Npop, list_off, U_SC) #on the S locus: SI to SC
            end


            #Computing Inbreeding depression and number of S alleles
            if phase == 1 && generation > 2000 - 50
                inDep = inbreedingDepression(Npop, list_ind, list_W, h,s, L)
                if abs(inDep - delta) > 1e-2 && Umax-Umin > 1e-2
                    println("ID=",inDep)
                    #Dichotomie
                    U, Umin, Umax = dichotomie(U,Umin,Umax, inDep, delta)
                    break
                end
                push!(list_ID,inDep)
                list_alleles = []
                for i in 1:Npop
                    S1, S2 = list_ind[i][1][1], list_ind[i][1][2]
                    append!(list_alleles,S1)
                    append!(list_alleles,S2)
                end
                nS = length(unique(list_alleles))
                push!(list_nbSalleles,nS)
            end

            
            list_ind, list_off = list_off, list_ind
            list_W, list_W_off = list_W_off, list_W

            
            #Filtering mutations
            num_mut = 0
            for i in 1:Npop
                num_mut +=  length(list_ind[i][2])+length(list_ind[i][3])
            end
            if (phase >= 1 && generation > 1 && generation%500==0) || num_mut > threshold_mutation_filter
                
                genome_mutation_counts = Dict()
                for i in 1:Npop
                    for mutation in list_ind[i][2]
                        count = get!(genome_mutation_counts, mutation, 0)
                        genome_mutation_counts[mutation] = count + 1
                    end
                    for mutation in list_ind[i][3]
                        count = get!(genome_mutation_counts, mutation, 0)
                        genome_mutation_counts[mutation] = count + 1
                    end
                end
                
                # Keep only mutations present fewer than 2*Npop times in the population
                for i in 1:Npop
                    list_ind[i][2] = [m for m in list_ind[i][2] if genome_mutation_counts[m] < 2*Npop]
                    list_ind[i][3] = [m for m in list_ind[i][3] if genome_mutation_counts[m] < 2*Npop]
                end
                num_mut = sum(length(list_ind[i][2])+length(list_ind[i][3]) for i in 1:Npop)
                threshold_mutation_filter = 2*num_mut
            end
   
            
            if phase == 0 && generation > 2000
                println("Going to phase 1")
                phase+= 1
                generation = 0
            end
            if phase == 1 && generation > 2000
                println("Going to phase 2")
                ID = sum(list_ID)/length(list_ID) #Inbreeding depression averaged over the last 50 generations
                nbSalleles = sum(list_nbSalleles)/length(list_nbSalleles) #number of S alleles before introducing the SC allele averaged over the last 50 generations
                phase+= 1
                phase2_not_launched = false
                generation = 0
            end
            
            #Computing the number and frequency of SC alleles + Stopping criteria (Stop after 50,000 generations if SC is fixed in the population)
            if phase == 2 && generation%100==0 && generation > 50000
                nbSC = 0
                for i in 1:Npop
                    S1, S2 = list_ind[i][1][1], list_ind[i][1][2]
                    if S1 < 0
                        nbSC +=1
                    end
                    if S2 < 0
                        nbSC+=1
                    end
                end
                
                push!(list_nbSC,nbSC)

                #Stopping criteria
                if nbSC == 2*Npop #SC is fixed in the population
                    nbSCalleles = nbSC
                    freqSC = 1
                    phase +=1 #Stop the simulation
                end
            
                if generation > max_gen-150
                    nbSCalleles = sum(list_nbSC[end-299:end])/300 #average of the nbSalleles over the last 30000 generations
                    freqSC = nbSCalleles/(2*Npop)
                    phase += 1 #Stops the simulation
                end
            end
                
        end


        #Complete the list with what to return
        if phase2_not_launched == false
            push!(Tot_nbSalleles,nbSalleles)
            push!(Tot_ID,ID)
            push!(Tot_nbSC,nbSCalleles)
            push!(Tot_freqSC,freqSC)
            push!(Tot_U, U)
            println("Done !")
        end
        
    end
    
    return (Tot_nbSalleles, Tot_ID, Tot_nbSC, Tot_freqSC, Tot_U)
end












#PARAMETERS

nbSim = 1

Npop = 1000
L = 10
h = 0.25
s= 0.05
alpha = 1 #Tested from 0 to 1
delta = 0.4 #Tested from 0 to 1
U_SC = 1e-4
U_SI = 1e-5
U = 5 #To change in order to control Inbreeding depression
dom = 0 #dominance of SC allele (0 if recessive or codominant, 1 if dominant  over the other SI alleles)



Tot_nbSalleles, Tot_ID, Tot_nbSC, Tot_freqSC, Tot_U = simulation(nbSim, Npop, h,s, alpha, U, U_SI, U_SC, L, dom, delta)


