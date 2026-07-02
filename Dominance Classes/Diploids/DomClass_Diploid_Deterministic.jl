##COMPATIBILITY AND PROBABILITY

function compatibility(parM, parF, list_class)
    SM1, SM2 = parM[1], parM[2]   # S alleles of the mother
    SF1, SF2 = parF[1], parF[2]   # S alleles of the father

    class_SF1 = list_class[SF1]   # Dominance class of the S alleles of the father 
    class_SF2 = list_class[SF2]

    # We artificially consider that all S alleles in the father are one of the most dominant one
    if class_SF1 > class_SF2
        SF2 = SF1
    elseif class_SF1 < class_SF2
        SF1 = SF2
    end

    if SM1 == 1 && SM2 == 1   # Mother is SC SC
        return 1

    elseif SM1 == 1           # Mother is SC Si
        if SM2 != SF1 && SM2 != SF2
            return 1
        end

    elseif SM2 == 1           # Mother is Si SC 
        if SM1 != SF1 && SM1 != SF2
            return 1
        end

    else                      # Mother is Si Sj 
        if SM1 != SF1 && SM1 != SF2 && SM2 != SF1 && SM2 != SF2
            return 1
        end
    end

    return 0
end


function probability(parM, parF, offspring)
    SM1, SM2 = parM[1], parM[2] #S alleles of the mother
    SF1, SF2 = parF[1], parF[2] #S alleles of the father
    Soff1, Soff2 = offspring[1], offspring[2] #S alleles of the offspring
    prob = 0
    for S1 in (SM1, SM2)
        for S2 in (SF1, SF2)
            if S1 == Soff1 && S2 == Soff2
                prob += 1
            end
        end
    end
    return prob/4
end




##SELFING RATE

function get_Wbar(tabFreq, delta)
    n = size(tabFreq[2],1)
    Wself =(1-delta)* sum(tabFreq[1][i,j] for i in 1:n, j in 1:n) #sum of fitness of individuals obtained through selfing
    Wout = sum(tabFreq[2][i,j] for i in 1:n, j in 1:n) #same for outcrossing
    return Wself + Wout
end


function selfingRate(par, Wbar, delta, alpha, list_class)
    S1, S2 = par[1], par[2] #S alleles of the parent
    c1, c2 = list_class[S1], list_class[S2] #Dominance classes for both S alleles
    classSC = list_class[1]

    if !(S1 == 1 || S2 == 1)
        return 0
    end

    theta = 0

    #Case SCSC
    if S1 == 1 && S2 == 1
        theta = 1
    #Case SCSi
    elseif S1 == 1
        if classSC > c2
            theta = 1
        end
    elseif S2 == 1
        if classSC > c1
            theta = 1
        end
    end

    W = 1
    if par[3] == 1 #The parent was obtained through selfing
        W = 1-delta
    end

    if theta == 0 || W == 0
        SR = 0
    else
        SR = (W*theta*alpha)/(W*theta*alpha+(1-alpha)*Wbar)
    end
    
    return SR
end


#FREQUENCIES AT THE NEXT GENERATION

#Return the entire tab of frequencies for all genotypes
#For each pair of parents, we look at what offspring they produce, and add that to the offspring's frequency
function freqNext_oc_fast(tabFreq, delta, alpha, list_class)
    x = tabFreq
    n = size(x[2],1)
    x_next_oc = zeros(n,n)

    Wbar = get_Wbar(x, delta)
    
    
    #For all maternal parent
    tab = zeros(n,n)
    for iM in 1:n
        for jM in 1:n
            for oM in (1,2)
                parM = (iM,jM,oM)
                aM = selfingRate(parM, Wbar, delta, alpha, list_class)
                if oM == 1
                    WM = 1-delta
                else
                    WM = 1
                end
                fitMom = x[oM][iM,jM]*(1-aM)*WM

                fill!(tab, 0.0)
                denom = 0.0

                #For all paternal parent
                for iF in 1:n
                    for jF in 1:n
                        for oF in (1,2)
                            xf = x[oF][iF,jF]
                            if xf == 0.0
                                continue
                            end
                            parF = (iF,jF,oF)
                            if oF == 1
                                WF = 1-delta
                            else
                                WF = 1
                            end
                                
                            compa_MF = compatibility(parM, parF, list_class)
                            if compa_MF == 0.0
                                continue
                            end
                            denom += WF*x[oF][iF,jF]*compa_MF
                            #All possible offspring from these parents
                            if iM == jM #homozygous
                                i_vals = (iM,)
                            else
                                i_vals = (iM,jM)
                            end
                            
                            if iF == jF
                                j_vals = (iF,)
                            else
                                j_vals = (iF,jF)
                            end
                            for i in i_vals
                                for j in j_vals
                                    offspring = (i,j)
                                    probij_MF = probability(parM, parF, offspring)
                                    tab[i,j] += WF*x[oF][iF,jF]*compa_MF*probij_MF
                                end
                            end
                        end
                    end
                end
                x_next_oc .+= fitMom.*(tab ./denom)
            end
        end
    end
    x_next_oc = x_next_oc./Wbar

    return x_next_oc
end


function freqNext_self_fast(tabFreq, delta, alpha, list_class)
    x = tabFreq
    n = size(x[2],1)
    x_next_self = zeros(n,n) 

    Wbar = get_Wbar(x, delta)

    for iM in 1:n
        for jM in 1:n
            for oM in (1,2)
                if oM == 1
                    WM = 1-delta
                else
                    WM = 1
                end
                parM = (iM,jM,oM)
                aM = selfingRate(parM, Wbar, delta, alpha, list_class)
                compa_M = compatibility(parM, parM, list_class)
                if compa_M == 0.0
                    continue
                end

                if iM == jM #homozygous
                    i_vals = (iM,)
                else
                    i_vals = (iM,jM)
                end
                for i in i_vals
                    for j in i_vals
                        offspring = (i,j)
                        probij_M = probability(parM, parM, offspring)
                        x_next_self[i,j] += aM*WM*x[oM][iM,jM]*compa_M*probij_M

                    end
                end
            end
        end
    end

    x_next_self = x_next_self./Wbar

    return x_next_self
end



## RECURSION

function get_freq(x, class_SC)
    n = size(x[2],1)
    list_class = [class_SC; list_class0]
    nbClass = length(unique(list_class))
    freq_classes = zeros(nbClass + 1) 
    nbInClasses = zeros(nbClass + 1) 
    freq_SC = 0
    x_self = x[1]
    x_out = x[2]
    freq_SC = sum(x_self[1,i]+x_self[i,1]+ x_out[1,i]+ x_out[i,1] for i in 1:n)/2
    freq_classes[1] = freq_SC
    
    for i in 1:n
        for j in 1:n
            if i!=1
                freq_classes[list_class[i]+1] += x_self[i,j]/2
                freq_classes[list_class[i]+1] += x_out[i,j]/2
            end
            if j!= 1
                freq_classes[list_class[j]+1] += x_self[i,j]/2
                freq_classes[list_class[j]+1] += x_out[i,j]/2
            end
        end
    end

    return freq_classes
end




function recursion_fast(delta, alpha, class_SC)

    list_class = [class_SC;list_class0]
    nbClass = length(unique(list_class))
    freq_classes = zeros(nbClass + 1) 
    #Initialization
    #All SI alleles have the same frequency, and SC allele is absent
    x = [fill(1/(2*nbSalleles^2), nbSalleles+1, nbSalleles+1) for _ in 1:2]
    for i in (1,2)
        for j in 1:nbSalleles +1
            x[i][1,j] = 0
            x[i][j,1] = 0
        end
    end

    s = sum(x[1]) + sum(x[2])


    #Phase 1: convergence of the SI alleles
    generation = 0
    phase = 0
    while generation < 1000
        generation +=1
        
        x_next = [zeros(nbSalleles+1, nbSalleles+1) for _ in 1:2]
        Wbar = get_Wbar(x, delta)
        
        x_next[1] = freqNext_self_fast(x, delta, alpha, list_class)
        x_next[2] = freqNext_oc_fast(x, delta, alpha, list_class)
        
        s = sum(x_next[1]) + sum(x_next[2])

        maxdiff = 0.0
        for i in 1:nbSalleles +1
            for j in 1:nbSalleles + 1
                x_next[1][i,j] /= s
                x_next[2][i,j] /= s
                maxdiff = max(maxdiff, abs(x_next[1][i,j]-x[1][i,j]), abs(x_next[2][i,j]-x[2][i,j]))
            end
        end
        x = x_next
        if maxdiff < 1e-5
            j = rand(2: nbSalleles+1)
            x[2][j,1] = 1/1000
            s = sum(x_next[1]) + sum(x_next[2])
            for i in 1:nbSalleles +1
                for j in 1:nbSalleles + 1
                  x_next[1][i,j] /= s
                  x_next[2][i,j] /= s
                end
            end
            freq_classes = get_freq(x, class_SC)
            generation = 0
            break
        end
    end


    while generation < 5000
        
        generation +=1
        
        x_next = [zeros(nbSalleles+1, nbSalleles+1) for _ in 1:2]
        Wbar = get_Wbar(x, delta)
        
        x_next[1] = freqNext_self_fast(x, delta, alpha, list_class)
        x_next[2] = freqNext_oc_fast(x, delta, alpha, list_class)

        s = sum(x_next[1]) + sum(x_next[2])
        
        maxdiff = 0.0
        for i in 1:nbSalleles +1
            for j in 1:nbSalleles + 1
              x_next[1][i,j] /= s
              x_next[2][i,j] /= s
              maxdiff = max(maxdiff, abs(x_next[1][i,j]-x[1][i,j]), abs(x_next[2][i,j]-x[2][i,j]))
            end
        end
        x = x_next
    end

    freq_classes = get_freq(x, class_SC)

    return x, freq_classes
end


function test_all(nbSalleles, class_SC, list_class0)
    list_freqSC = []
    list_x = []
    for j in 0:10
        println("alpha=",j*0.1)
        alpha = j*0.1
        list_alpha = []
        list_tot = []
        for i in 0:100
            delta = i*0.01
            x, freq_classes = recursion_fast(delta, alpha, class_SC)
            push!(list_alpha, freq_classes[1])
            push!(list_tot, x)
        end
        push!(list_freqSC, list_alpha)
        push!(list_x, list_tot)
    end
    return list_freqSC, list_x
end






##PARAMETERS

nbSalleles = 15
class_SC= 1 #Tested classes: 1,2,3,4,5
list_class0 = [1,2,2,3,3,3,3,3,4,4,4,4,4,4,4,4]



list_freqSC, list_x = test_all(nbSalleles, class_SC, list_class0)


    
