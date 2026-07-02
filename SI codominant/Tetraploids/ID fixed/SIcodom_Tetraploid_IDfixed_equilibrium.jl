using Distributions
using StatsBase
using Random



## MUTATIONS 


#Mutations Si to Sj 
function mutation_SI!(Npop, list_ind, U_SI, k)
    nbMut = rand(Binomial(4*Npop,U_SI)) #number of mutations among all S alleles
    for _ in 1:nbMut
        indMut = rand(1:Npop) #choosing the individual
        #choosing the allele of that individual that will mutate
        chrMut = rand(1:4)
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
    nbMut = rand(Binomial(4*Npop,U_SC)) #number of mutations among all S alleles
    for _ in 1:nbMut
        indMut = rand(1:Npop) #choosing the individual
        #choosing the allele of that individual that will mutate
        chrMut = rand(1:4)
        allMut = list_ind[indMut][1][chrMut]
        if allMut != -1 #mutates only if the allele is SI
            list_ind[indMut][1][chrMut] = -1 #Replacing with mutant SC allele
        end
    end
end


# RECOMBINATION 

function recombination(par, L)
    #number and positions of the crossover
    nb_crossovers = rand(Poisson(L))
    pos_crossovers = rand(nb_crossovers)
    sort!(pos_crossovers)

    #recombined chromosomes
    
    list_chr = [par[2],par[3], par[4], par[5]]
    list_Salleles = [par[1][1], par[1][2], par[1][3], par[1][4]]
    first_possible_index = [1,1,1,1]
    
    #Shuffling the chromosomes
    chr1 = list_chr[1]
    chr2 = list_chr[2]
    chr3 = list_chr[3]
    chr4 = list_chr[4]
    Sallele1 = list_Salleles[1]
    Sallele2 = list_Salleles[2]
    

    gamete_chr1::Vector{Float64} = []
    sizehint!(gamete_chr1, floor((length(chr1)+length(chr2)+length(chr3)+length(chr4))/4)) #Suggest that gamete reserve some capacity in advance. This can improve performance.
    gamete_chr2::Vector{Float64} = []
    sizehint!(gamete_chr2, floor((length(chr1)+length(chr2)+length(chr3)+length(chr4))/4)) 
    Sallele_chr1 = list_Salleles[1]
    Sallele_chr2 = list_Salleles[2]
    
    L = [1,2,3,4]
    for i in 1:nb_crossovers + 1
        
        #Shuffling the chromosomes
        shuffle!(L)
        chr1 = list_chr[L[1]]
        chr2 = list_chr[L[2]]
        Sallele1 = list_Salleles[L[1]]
        Sallele2 = list_Salleles[L[2]]
        
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


        if j <= 0.5 && k >=0.5
            Sallele_chr1 = Sallele1
            Sallele_chr2 = Sallele2
        end
        for l in first_possible_index[L[1]]:length(chr1)
            if chr1[l] >= k
                first_possible_index[L[1]] = l
                break
            end
            if chr1[l] >= j
                push!(gamete_chr1,chr1[l])
            end
        end
        for l in first_possible_index[L[2]]:length(chr2)
            if chr2[l] >= k
                first_possible_index[L[2]] = l
                break
            end
            if chr2[l] >= j
                push!(gamete_chr2,chr2[l])
            end
        end
 
        
    end

    return (Sallele_chr1, Sallele_chr2, gamete_chr1, gamete_chr2) #Changed for tetraploids
end




## SELFING RATE 

function selfing_rate(Npop, list_ind, par, alpha, list_W, dominance)
    Salleles = list_ind[par][1] #list of the parent's S-alleles

    SR = 0.0
    #Special cases
    if Salleles[1] != -1 && Salleles[2] != -1 && Salleles[3] != -1 && Salleles[4] != -1 #the individual is SI
        return 0.0
    elseif Salleles[1] == -1 && Salleles[2] == -1 && Salleles[3] == -1 && Salleles[4] == -1 ##Four SC alleles, compatible with everyone
        gamma_ii = 2
        gamma_ij = 2
        sum_val = gamma_ij*(sum(list_W)-list_W[par])
        SR = (alpha*list_W[par]*gamma_ii)/(alpha*list_W[par]*gamma_ii + ((1-alpha)/(Npop-1))*sum_val)
    
        
    #Case 1,2 or 3 SC
    else

        #SC DOMINANT (if SC codominant or recessive, SR=0)
        if dominance ==1
            gamma_ii = 2

            Salleles_SI = filter(!=(-1), Salleles) #We keep only the SI alleles of the list Salleles
    
            sum_val = 0
            for j in 1:Npop
                gamma_ij = 0 
                
                if j != par
                    Sj = list_ind[j][1] #S alleles in the pollen
                    gamma_ij = 2 #When at least 1 SC in the pollen, or when all SI alleles from the pollen are different from the SI alleles of the mother
                    if Sj[1] != -1 && Sj[2] != -1 && Sj[3] != -1 && Sj[4] != -1  #No SC allele in the pollen
                        for k in 1:length(Salleles_SI) #If 1 SI allele in the pollen is shared with the mother -> Not compatible 
                            if Sj[1] == Salleles_SI[k] || Sj[2] == Salleles_SI[k] || Sj[3] == Salleles_SI[k] || Sj[4] == Salleles_SI[k]
                                gamma_ij = 0
                                break
                            end
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

function reproduction!(Npop, list_ind, list_W, list_off, list_W_off, alpha, delta, L, phase, dominance)
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
        fit_off = 1.0
        OK_father = false

        #OUTCROSSING
        if selfRate <=  rand()
            OK_father = false
            maxAttempts = 50
            loop = 0
            while OK_father == false && loop < maxAttempts
                loop +=1
                
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
                    if Salleles_par2[1] == -1 || Salleles_par2[2] == -1 || Salleles_par2[3] == -1 || Salleles_par2[4] == -1 # At least one allele in the pollen is SC
                        OK_father = true
                        break
                    end
                    if Salleles_par1[1] == Salleles_par2[1] || Salleles_par1[2] == Salleles_par2[1] || Salleles_par1[3] == Salleles_par2[1] || Salleles_par1[4] == Salleles_par2[1] # OK first S allele of the pollen
                        continue
                    end
                    if Salleles_par1[1] == Salleles_par2[2] || Salleles_par1[2] == Salleles_par2[2] || Salleles_par1[3] == Salleles_par2[2] || Salleles_par1[4] == Salleles_par2[2] # OK second S allele of the pollen
                        continue
                    end
                    if Salleles_par1[1] == Salleles_par2[3] || Salleles_par1[2] == Salleles_par2[3] || Salleles_par1[3] == Salleles_par2[3] || Salleles_par1[4] == Salleles_par2[3] # OK third S allele of the pollen
                        continue
                    end
                    if Salleles_par1[1] != Salleles_par2[4] && Salleles_par1[2] != Salleles_par2[4] && Salleles_par1[3] != Salleles_par2[4] && Salleles_par1[4] != Salleles_par2[4] #OK fourth S allele of the pollen
                        OK_father = true
                    end

                
                #SC RECESSIVE OR CODOMINANT
                elseif dominance == 0
                    if Salleles_par1[1] == -1 && Salleles_par1[2] == -1 && Salleles_par1[3] == -1 && Salleles_par1[4] == -1 #All alleles in the mother are SC
                        OK_father = true
                        break
                    end
                    if Salleles_par2[1] == -1 && Salleles_par2[2] == -1 && Salleles_par2[3] == -1 && Salleles_par2[4] == -1 # All alleles in the pollen are SC
                        OK_father = true
                        break
                    else 
                        Salleles_par2_SI = filter(!=(-1), Salleles_par2) #Contains only SI alleles of the pollen
                        nb_incomp = 0 #will contain the number of SI alleles shared by the pollen and the mother
                        for i in 1:length(Salleles_par2_SI)
                            if Salleles_par2_SI[i] == Salleles_par1[1] || Salleles_par2_SI[i] == Salleles_par1[2] || Salleles_par2_SI[i] == Salleles_par1[3] || Salleles_par2_SI[i] == Salleles_par1[4]   #The i-th SI allele of the pollen is the same as 1 S-allele of the mother
                                nb_incomp +=1
                            end
                        end 
                        if nb_incomp == 0
                            OK_father = true
                            break
                        end
                    end
                end
                    
            end
            
        #SELFING
        #All pollen grains are compatible when SC
        else
            par2 = par1
            fit_off = 1-delta
            OK_father = true
        end

        
        #Recombination
        if OK_father == true
                Sallele1, Sallele2, chr1, chr2 = recombination(list_ind[par1], L)
                Sallele3, Sallele4, chr3, chr4 = recombination(list_ind[par2], L)
                #Creating a new individual for the next generation
                Npop_off += 1
                offspring = [[Sallele1, Sallele2, Sallele3, Sallele4], chr1, chr2, chr3, chr4]
                list_off[Npop_off] = offspring
                list_W_off[Npop_off] = fit_off
        end
    end
end





## MAIN 

#Phases:
# - Only mutations Si to Sj (2000 generations)
# - Then mutations on the genome (2000 generations)
# - Then mutations Si to SC (max 200000 generations, and stop when invasion of SC)

function simulation(nbSim, Npop, alpha, delta, U_SI, U_SC, L, dominance)
    #Creating output lists (for all runs)
    Tot_nbSalleles::Vector{Float64} = []
    Tot_ID::Vector{Float64} = []
    Tot_nbSC::Vector{Float64} = []
    Tot_freqSC::Vector{Float64} = []

    for simu in 1:nbSim
        nbSalleles = 0
        nbSCalleles = 0
        freqSC = 0.0
        ID = 0
        
        it = 0
        println("nbSim=" , simu)
        #Creating intermediate lists
        list_ID::Vector{Float64}, list_nbSalleles::Vector{Int} = [],[]
        list_nbSC::Vector{Int} = []


        #Initialisation
        nbS = 100 #number of possible S alleles
        list_ind::Vector{Vector{Vector{Float64}}} = []
        for i in 1:Npop
            ind::Vector{Vector{Float64}} = [[0.0,0.0,0.0,0.0], [],[],[],[]] #4 empty lists to represent 4 chromosomes
            push!(list_ind,ind)
            #Initialising the 4 S-alleles
            list_ind[i][1][1] = rand(0:nbS-1)
            list_ind[i][1][2] = rand(0:nbS-1)
            list_ind[i][1][3] = rand(0:nbS-1) 
            list_ind[i][1][4] = rand(0:nbS-1) 
        end
        list_W = [1.0 for _ in 1:Npop]
        
        list_off::Vector{Vector{Vector{Float64}}} = []
        for _ in 1:Npop
            ind::Vector{Vector{Float64}} = [[0.0,0.0,0.0,0.0], [],[],[],[]]
            push!(list_off,ind)
        end
        list_W_off = Vector{Float64}(undef, Npop)


        generation = 0
        phase = 0
        
        max_gen = 250000

        #Phase 0: Starting with mutations on the S locus only (SI to SI)
        #Phase 1: Introduction of deleterious mutations on the genome
        #Phase 2: Introduction of S locus mutations to SC alleles
        
        while phase <= 2 
            generation += 1
            if generation%1000==0
                println("generation=", generation)
            end


            #Creating the next generation
            reproduction!(Npop, list_ind, list_W, list_off, list_W_off, alpha, delta, L, phase, dominance)
            
            #Mutations
            mutation_SI!(Npop, list_off, U_SI, nbS) #on the S locus: SI to SI

            if phase >= 2
                mutation_SC!(Npop, list_off, U_SC) #on the S locus: SI to SC
            end

            
            #Computing Inbreeding depression and number of S alleles
            if phase == 1 && generation > 2000 - 50
                list_alleles = []
                for i in 1:Npop
                    S1, S2, S3, S4 = list_ind[i][1][1], list_ind[i][1][2], list_ind[i][1][3], list_ind[i][1][4]
                    append!(list_alleles,S1)
                    append!(list_alleles,S2)
                    append!(list_alleles,S3)
                    append!(list_alleles,S4)
                end
                nS = length(unique(list_alleles))
                push!(list_nbSalleles,nS)
            end

            
            list_ind, list_off = list_off, list_ind
            list_W, list_W_off = list_W_off, list_W 
           
            
            if phase == 0 && generation > 2000
                println("Going to phase 1")
                phase+= 1
                generation = 0
            end
            if phase == 1 && generation > 2000
                println("Going to phase 2")
                nbSalleles = sum(list_nbSalleles)/length(list_nbSalleles) #number of S alleles before introducing the SC allele averaged over the last 50 generations
                phase+= 1
                generation = 0
            end
            
            #Computing the number and frequency of SC alleles + Stopping criteria (Stop after 50,000 generations if SC is fixed in the population)
            if phase == 2 && generation%100==0 && generation > 50000
                nbSC = 0
                for i in 1:Npop
                    S1, S2, S3, S4 = list_ind[i][1][1], list_ind[i][1][2], list_ind[i][1][3], list_ind[i][1][4]
                    if S1 < 0
                        nbSC +=1
                    end
                    if S2 < 0
                        nbSC+=1
                    end
                    if S3 < 0
                        nbSC+=1
                    end
                    if S4 < 0
                        nbSC+=1
                    end
                end
                
                push!(list_nbSC,nbSC)

                #Stopping criteria
                if nbSC == 4*Npop #SC is fixed in the population
                    nbSCalleles = nbSC
                    freqSC = 1
                    phase +=1 #Stop the simulation
                end
            
                if generation > max_gen-150
                    nbSCalleles = sum(list_nbSC[end-299:end])/300 #average of the nbSalleles over the last 30000 generations
                    freqSC = nbSCalleles/(4*Npop)
                    phase += 1 #Stops the simulation
                end
            end
                
        end


        #Complete the list with what to return
        push!(Tot_nbSalleles,nbSalleles)
        push!(Tot_ID,ID)
        push!(Tot_nbSC,nbSCalleles)
        push!(Tot_freqSC,freqSC)
        println("Done !")


    end
    
    return (Tot_nbSalleles, Tot_ID, Tot_nbSC, Tot_freqSC)
end










#PARAMETERS

nbSim = 10

Npop = 1000
L = 10
alpha = 0.1 #Tested between 0 and 1
U_SC = 1e-4
U_SI = 1e-5
dom = 0 #dominance of SC allele (0 if recessive or codominant, 1 if dominant  over the other SI alleles)
delta = 0.3 #Tested between 0 and 1



Tot_nbSalleles, Tot_ID, Tot_nbSC, Tot_freqSC = simulation(nbSim, Npop, alpha, delta, U_SI, U_SC, L, dom)

