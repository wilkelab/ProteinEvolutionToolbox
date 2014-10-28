##This code renumbers pdb protein file structure (ATOM data accessible through "pdb_name$atom", which is a matrix table with all the ATOM data info), if inserts are present (numbering of residues has "A","B","C"... in the "inserts" column of "pdb_name$atom". This step is required when using DSSP function because DSSP function takes in the numbering from "resno" column of pdb ATOM data "pdb_name$atom" table, and does not consider the "inserts" column from the same table. That way the "resno" vector has repeating numbering, like this resno <- c(1,2,3,4,4,5,6,7,8,8,9).	
require(io3d)

#label identities of pbd sturctures present int he alignment "aln". 
ids <- aln$id

#for loop runs through each pdb structure present in the alignment and checks weather the structure has inserts. If inserts are present the renumbering happens in the second for loop in this for loop. 
for (i in 1:length(ids)){

	#get raw (unaltered pdb file acquired through "get.pdb()" function) pdb into a file "raw_pdbs". 
	get.pdb(ids[i], path = "raw_pdbs/")
        
        #read in and lable just downloaded structure into R. 
	pdb <- read.pdb(paste("raw_pdbs/",ids[i],sep = ""))
    
        #index out and label the "insert" column from "pdb$atom" table.
	insert <- pdb$atom[,"insert"]
    
        #if insert(s) present, renumber.
	if (!all(is.na(insert))){
        
        	#index out and label all the places (or rows) in the insert column where inserts are present, by using "which()" function. This function returns a vector with indecies of the TRUE expression in the argument. Here this expression finds all the places where "insert" column does not have "NA".
		inds_inserts_all <- which(!is.na(insert))

        	#find TRUE/FALSE vector where indecies of inserts in the "inds_inserts_all" interrupt and are not consequtive anymore. 
		#ex: if inds_inserts_all <- c(3,4,8,9),then a <- c(TRUE,FALSE,TRUE)
        	a <- inds_inserts_all[2:length(inds_inserts_all)] == inds_inserts_all[1:length(inds_inserts_all)-1]+1

        	#make numeric and label residue numbering vector from "pdb$atom" ATOM table, to adjust for new numbering.
 		nums <- as.numeric(pdb$atom[,"resno"])
        
       		#if only one insert present (found by checking weather vector "a" has only TRUE and there should be only one TRUE), renumber (add 1) starting at the insert position and until the end of "nums" vector.
  		#ex: nums <- c(1,2,3,3,4,5), renums <- c(3,4,5), new renums <- (4,5,6), new nums <- c(1,2,3,4,5,6).
       		if (all(a)==TRUE){
                        
			#find the insert position start from all "inds_inserts_all" where all the indecies of inserts are present. In the case of one insert this vector would have the length of one, and contain one index value of the insert.
            		start_num = inds_inserts_all[1]
            
			#use the starting number to index out all the numbers needed to be renumbered in the "nums" vector
            		renums = nums[start_num:length(nums)]
    
                        #renumber this values by adding 1.
            		renums = renums+1
                       
                        #assign the renumbered vector to the original number vector "nums" starting at the insert position and until the end of the vector.
            		nums[start_num:length(nums)] = renums
        	}

                #else, if more than one insert is present 
        	else{
			
			#add FALSE value at the beginning of the "a" vector, because the "a" vector is cut by first value when the consequetiveness is compared to find the vector itself. This first FALSE value will indicate the interuption in the indecies, like all other FALSES in this vector do.  
            		a <- append(a,FALSE,0)
            
			#index out and label all the indecies where renumbering needs to start from all the inserts indecies vector "inds_inserts_all" using the locations of negated or the opposite of "a" vector. The negated or opposite vector of "a" ("!a") has now FALSE values for TRUE, and TRUE values for FALSE.
           	 	inds_inserts_start <- inds_inserts_all[!a]
        
            		#for loop renumberes the "nums" vector starting at each insert postion from "inds_inserts_start" until the following insert position. 
              		for (j in 1:(length(inds_inserts_start)-1)){
				
				#index out of nums all the numbers to be renumbered, starting at the insert position present until the following insert position. Both positions are acquired through "inds_inserts_start" vector.  
                		renums = nums[inds_inserts_start[j]:(inds_inserts_start[j+1]-1)]
			
				#here renumbering depends on the number of inserts. If two inserts are present, numbers between first and second insert are need to be increased by one, while the numbers between second and third inserts need to be increased by 2. Value "j" counts the loops, and by how much each interval need to be adjusted by (i.e. each loop adjust the intervals by the number of that loop)
                		renums = renums+j
				
				#assign the renumbered vector to the original number vector "nums" starting at the insert position and until the following insert position.  
                		nums[inds_inserts_start[j]:(inds_inserts_start[j+1]-1)] = renums
			}
            	}

        #write the new renumbered pdb
	write.pdb(pdb=pdb,xyz = pdb$xyz, resno = nums, resid = pdb$atom[,"resid"], eleno = pdb$atom[,"eleno"], elety = pdb$atom[,"elety"], chain = pdb$atom[,"chain"], o = pdb$atom[,"o"], b = pdb$atom[,"b"], file = paste("raw_pdbs/",ids[i],sep=""))

	#split the whole pdb struture into chain pdb structures
	#split.pdb(paste("raw_pdbs/",ids[i],sep=""),path = "raw_pdbs/split_chain")

	}
	#if no inserts are present in the structure, move on to the next structure in the alignment.
	else next 
}
