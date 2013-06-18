##This code calculates RSA for all residues of protein structures in the alignment using the DSSP function, and normalizing values from M. Tien, A. G. Meyer, S. J. Spielman, C. O. Wilke. Maximum allowed solvent accessibilites of residues in proteins. 

#read in the final alingment for a particular structure. Here it's set to Hemagglutinin (chain A) 1RD8_A.pdb
aln <- read.fasta("Phy/1RD8_A.fa")

#label the protein structures present in that alignment
ids <- aln$id

#dictionary with  Normalization Values
NormValLst <-list("A"=129 ,"R"=274,"N"=195,"D"=193,"C"=158,"E"=223,"Q"=224,"G"=104,"H"=209,"I"=197,"L"=201,"K"=237,"M"=218,"F"=239,"P"=159,"S"=151,"T"=172,"W"=282,"Y"=263,"V"=174)

#recreate an RSA matrix from the alignment matrix. Both matrices contain same gap positions, as well as identical dimensions.
RSA <- matrix(NA,length(aln$ali[,1]),length(aln$ali[1,]))

#NOTE: For some protein pdb structures (specifically Hemmagglutinin in our case) the numbering is not exactly consecutive, like this 1,2,3,..n. Instead it uses this type of numbering 1,2,2A,3,3A,4..n. The reason for this is that some residues in the structure were discovered later, and therefore were inserted into this structure later. The function DSSP does not work with the second type of numbering, hence the pdb structures need to be rewritten. This is done in a separate program RenumPDB.r.  

#get SA and RSA for each protein, while stripping the protein structure of other organic molecules.  
for (i in 1:length(ids)){
        #read in a pdb chain structure 
	raw_pdb <- read.pdb(paste("raw_pdbs/split_chain/",ids[i],sep=""))
	
        #keep only the protein structure, by locating the residuals of that protein through presence of C-alpha in the ATOM data "pdb_name$atom". 

	#select atoms that are only C-alpha atoms ("calpha"), no random ATOMs such as HOH, NAG, PO4 etc.
        #inds_CA is a vector of indecies (locations) of C-alphas in this pdb structure in the ATOM data.  
	inds_CA <- atom.select(raw_pdb,"calpha")
	
	#get residual names from C-alpha indices, by indexing those names out of pdb ATOM data, column 4, where all the residual names are present
	resid <- raw_pdb$atom[inds_CA$atom,][,4]

	#get the indicies of residuals from residual names. 
	inds_resid <- atom.select(raw_pdb, resid = unique(resid))

	#get new ATOM data using residual names indices, to keep only the residual atom data 
	atom <- raw_pdb$atom[inds_resid$atom,]

	#get the XYZ data from the same residual names indices used in previous step. Again only residual XYZ data is kept.
	xyz <- raw_pdb$xyz[inds_resid$xyz]

	#rewrite the pdb file with new ATOM data.
	write.pdb(pdb=raw_pdb,xyz = xyz, resno = atom[,6], resid = atom[,4], eleno = atom[,1], elety = atom[,2], chain = atom[,5], o = atom[,11], b = atom[,12], file = paste("raw_pdbs/split_chain/clean_pdbs/",ids[i],sep=""))
	
	#read in the new ATOM data pdb structure.
	pdb <- read.pdb(paste("raw_pdbs/split_chain/clean_pdbs/",ids[i],sep =""))

	#select pdb structure sequence.
	seq <- seq.pdb(pdb)

	#run DSSP on this structure.
	DSSP <- dssp(pdb)

	#get the SA from this function. "$acc" accesses SA from the returned results of DSSP function. 
        SA <- DSSP$acc
	k = 1 #assign k for indexing purpose in the SA value vector that is independent from all other counting in the loops.

	#for loop runs through the each i,j position in the RSA matrix and fills them in with RSA values for a correspoding residue position in the alignment matrix "aln". If gap is present in the alignment matrix, RSA is assigned "NA" value.  
	for (j in 1:length(RSA[1,])){
                
		#if a gap is not present in the alignment continue with calculating RSA. 
		if (aln$ali[i,j] != "-"){

                        #get the normalization value from the Normalization Value Dictionary based on the residual identity.
                        #return "NULL" if the residual is not present in this Dictionary
			NormVal <- NormValLst[[seq[k]]]
			
			#check if the residual exists based on the value returned. Assign NA to the RSA value if the residual does not exist 
			if (is.null(NormVal)) RSA_val = NA

			#if residual exists assign the RSA value equaling to SA divided by Normalization Value
			else RSA_val = SA[k]/NormVal

		k = k+1 #adjust k to indicate the further move along the SA vector.
               
                #Place the RSA value into the RSA matrix based on the residual location in the alignment matrix
		RSA[i,j] = RSA_val
		}	
                
                #if the gap is present move to the next row location in the alignment
		else next
	}
}



					





			
	
