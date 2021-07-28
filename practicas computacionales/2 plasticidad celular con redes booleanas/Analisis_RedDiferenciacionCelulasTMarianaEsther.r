# Título: Plasticidad celular con redes booleanas 



# Nombre: Biología de sistemas



# Fecha: febrero 2021

####################################################################################



# cargar librería

library(BoolNet) 
library(igraph)


# cargar toda la red

net <- loadNetwork("RedDiferenciacionCelulasTMarianaEsther.txt") 

net



# mostrar la red

plotNetworkWiring(net)
plotNetworkWiring(net,layout=layout.circle, main="circle")
# for more information:
#https://www.r-graph-gallery.com/247-network-chart-layouts.html


# obtener atractores

# esto da mas de 320 atractores.. 

attr <- getAttractors(net)

attr



# mostrar atractores

#plotAttractors(attr)



##### Pregunta 1: Ambiente basal ####################

# fijar las condiciones ambientales

net_no_inputs <- fixGenes(net, c("IFNGe", "IL2e", "IL4e", "IL21e", "TGFBe", "IL10e"), c(0,0,0,0,0,0))



# obtener los atractores

attr_no_inputs <- getAttractors(net_no_inputs) 



#TBET IFNG GATA3 IL2 IL4 PORGT IL21 FOXP3 TGFB IL10 BCL6 IL9 IFNGe IL2e IL4e IL21e TGFBe IL10e



number_attractors = max(attr_no_inputs$stateInfo$attractorAssignment)

number_attractors



# Contar el número de atractores con Gata3 = 0

number_attractors = max(attr_no_inputs$stateInfo$attractorAssignment)

counter_success = 0;



for(AttNumber in 1:number_attractors) {

        

	AA = getAttractorSequence(attr_no_inputs,AttNumber)

	A = AA[ ,c(3)] # we care only about the third element, which is Gata3

        

	if (sum(A==1) > 0){ # algun elemento es 1

	        

           print("oh no, gata3 is on! :(")

	        

          } else {

                  

           counter_success=counter_success + 1

           

         }

}



number_attractors



counter_success



# cuantificar el tamaño total de las cuencas de atracción con Gata3 = 0



BasinSize_GataOff = 0

BasinSize_GataOn = 0



counter_success = 0



for(AttNumber in 1:number_attractors) {

        

	AA = getAttractorSequence(attr_no_inputs,AttNumber)

	SizeBasinCurrentAttractor = attr_no_inputs$attractors[[AttNumber]]$basinSize

	#### PARTE DE CODIGO ADECUADA PARA PERMITIR EVALUAR ATRACTORES CÍCILICOS

	A = AA[,3] # we care only about the third element, which is Gata3

        if (sum(A==1)){ 

                

           print("oh no, gata3 is on! :(")

           BasinSize_GataOn=BasinSize_GataOn+SizeBasinCurrentAttractor

           

        } else {

                

           counter_success=counter_success+1

           BasinSize_GataOff=BasinSize_GataOff+SizeBasinCurrentAttractor

        }

}



BasinSize_GataOn =  BasinSize_GataOn/2^12

BasinSize_GataOn



BasinSize_GataOff = BasinSize_GataOff/2^12

BasinSize_GataOff





##### Pregunta 2: Ambientes pro -inflamatorios ####################



#lo mismo que arriba - pero con IL4 y o IL2 prendidos





# uncomment if you want to assess the effect of IL4 only:

# net_pro_inflammatory_env <- fixGenes(net, c("IFNGe", "IL2e", "IL4e", "IL21e", "TGFBe", "IL10e"), c(0,0,1,0,0,0))

# uncomment if you want to assess the effect of IL2 only:

# net_pro_inflammatory_env <- fixGenes(net, c("IFNGe", "IL2e", "IL4e", "IL21e", "TGFBe", "IL10e"), c(0,1,0,0,0,0))

# uncomment if you want to assess the effect of IL2 and IL4 combined:

net_pro_inflammatory_env <- fixGenes(net, c("IFNGe", "IL2e", "IL4e", "IL21e", "TGFBe", "IL10e"), c(0,1,1,0,0,0))

attr_pro_inflammatory_env <- getAttractors(net_pro_inflammatory_env) 

number_attractors=max(attr_pro_inflammatory_env$stateInfo$attractorAssignment)



BasinSize_GataOff_pro_inflammatory_env = 0

BasinSize_GataOn_pro_inflammatory_env = 0



counter_success_pro_inflammatory_env = 0



for(AttNumber in 1:number_attractors) {

        

	AA = getAttractorSequence(attr_pro_inflammatory_env,AttNumber)

	SizeBasinCurrentAttractor = attr_pro_inflammatory_env$attractors[[AttNumber]]$basinSize

	A = AA[1,c(3)] # we care only about the thid element, which is Gata3

	#print(AA)

        #print(A)

        if (A==1){ 

                

           print("oh yeah, gata3 is on! :)")

           BasinSize_GataOn_pro_inflammatory_env =BasinSize_GataOn_pro_inflammatory_env +SizeBasinCurrentAttractor

        

        } else {

                

            BasinSize_GataOff_pro_inflammatory_env =BasinSize_GataOff_pro_inflammatory_env +SizeBasinCurrentAttractor

        }

}



BasinSize_GataOn_pro_inflammatory_env = BasinSize_GataOn_pro_inflammatory_env/2^12

# 2^(N-K) K: número de "genes" fijos (N=18, K=6)

BasinSize_GataOn_pro_inflammatory_env 



BasinSize_GataOff_pro_inflammatory_env = BasinSize_GataOff_pro_inflammatory_env/2^12

BasinSize_GataOff_pro_inflammatory_env 





############ Tercer pregunta : cómo pasar de un atractor a otro, empujando con condiciones iniciales #####



# Definir una condición inicial apropiada (\verb|ci=getAttractorSequence(attr_no_inputs,1)|)\verb|ci[,c(14,15)]=1)

# Obtener la dinámica de esta condición inicial a su atractor (\verb|getPathToAttractor(net_pro_inflammatory_env, ci)|)

# Checar si el atractor al que converge tiene Gata3 = 1 (\verb|jj=length(path[,3])|), \verb|path[jj,3]==1|)

# Checar si este atractor es uno de los atractores cuando IL2e y IL4e están prendidos (pregunta 2)

# attr_pro_inflammatory_env

# Ver en qué paso entró a la vasija de atracción de este atractor





time_to_attr = matrix(0,1, 13) # preallocate; sabemos que son  13  los atractores en condiciones sin inflamación

attractor_off_converges_to_attractor_on = 0 # set the counter 

for (attrNum in 1:13) {

        

        print(attrNum)

        

        # 1. Definir una condición inicial apropiada:

        print('we are at atractor:')

        ci = getAttractorSequence(attr_no_inputs,attrNum) # elejir un atractor del 1 al 13 de esta subred

        ci

        ci[,c(14,15)] = 1 # fijar IL2e y IL4e=1

        

        # 2. Obtener la dinámica de esta condición inicial a su atractor 

        print('the path to the attractor is:')

        path = getPathToAttractor(net_pro_inflammatory_env, ci)

        

        # 3 checar cuántos pasos al atractor

        jj = length(path[,3])

        jj

        

        # 4 Checar si el atractor al que converge tiene Gata3 = 1 

        path[jj,3] == 1

        path[jj,3]

        

        # 5. Checar si este atractor es uno de los atractores cuando IL2e y IL4e están prendidos (pregunta 2)

        number_attractors=max(attr_pro_inflammatory_env$stateInfo$attractorAssignment)

        path[jj,] == 1

        unknown_attractor_number = -1

        att_number = 1

        

        while(unknown_attractor_number == -1 & att_number < number_attractors + 1) {

                

                comparison= path[jj,]==getAttractorSequence(attr_pro_inflammatory_env,att_number)

                

                if (length(which((comparison)==FALSE))==0) {

                        

                	print("yeah - under IL2 and IL4 forcing the system converges to attractor:")

                	print( att_number) 

                	unknown_attractor_number=att_number

                	attractor_off_converges_to_attractor_on=attractor_off_converges_to_attractor_on + 1

                	

        	} else {

                        att_number=att_number+1

                

                }

        }

        

        # 5. En qué paso entró a la cuenca de atracción de este atractor?

        att_basin=getBasinOfAttraction(attr_pro_inflammatory_env, att_number)

        number_elements_in_basin=length(att_basin$initialState.TBET)

        max_steps=length(path[,3])

        

        # vamos a barrer primero por cada uno de los pasos

        inside_basin = -1

        step = 1

        basin_element = 1

        

        while (inside_basin == -1 & step < max_steps+1) {

                

                while (inside_basin==-1 & basin_element<number_elements_in_basin+1) {

                        comparison2=(path[step,]==c(att_basin[basin_element,]$initialState.TBET, att_basin[basin_element,]$initialState.IFNG, att_basin[basin_element,]$initialState.GATA3, att_basin[basin_element,]$initialState.IL2, att_basin[basin_element,]$initialState.IL4, att_basin[basin_element,]$initialState.PORGT, att_basin[basin_element,]$initialState.IL21, att_basin[basin_element,]$initialState.FOXP3, att_basin[basin_element,]$initialState.TGFB, att_basin[basin_element,]$initialState.IL10, att_basin[basin_element,]$initialState.BCL6, att_basin[basin_element,]$initialState.IL9, att_basin[basin_element,]$initialState.IFNGe, att_basin[basin_element,]$initialState.IL2e, att_basin[basin_element,]$initialState.IL4e, att_basin[basin_element,]$initialState.IL21e, att_basin[basin_element,]$initialState.TGFBe, att_basin[basin_element,]$initialState.IL10e))

                        if  (length(which((comparison2)==FALSE))==0) {

                                print("yeah!!! time of the transient found!!")

                                print(step-1)

                                print("basin element")

                                print(basin_element)

                                inside_basin=1

                        } else {

                        (basin_element= basin_element+1)}

                }

                step=step+1

        }

        time_to_attr[attrNum]=step-2

}







mean(time_to_attr)

attractor_off_converges_to_attractor_on/13