#Andrew Dhawan
#Sept 2017
#Contains the code to generate a linear model of miRNA predicting mRNA signature scores on a dataset
#Uses 10x cross-validation/penalised linear regression

library(RankProd)
library(reshape2)
library(penalized)

get_coefficients_pre_filter <- function(cancer_type,scores){
	#this function actually does the linear modelling, whilst pre-filtering out predictors that we know will function poorly

	#load in the miRNA data
	miRNA_data <- all_miRNA_datasets[[cancer_type]]

	#take only common subset of miRNA and scores
	common_colNames <- intersect(colnames(miRNA_data),names(scores))

	#take just the common pieces
	miRNA_data <- miRNA_data[,common_colNames]
	scores <- scores[common_colNames]

	#z-transform the scores
	scores <- as.numeric(scores) - mean(as.numeric(scores))/sd(as.numeric(scores))

	#expression filter for miRNA
	expression_threshold <- 0.80 # means that at least 20% of samples must have a nonzero value of the mRNA
	miRNA_data <-miRNA_data[which((rowSums(miRNA_data==0)) < ((1-expression_threshold) * length(colnames(miRNA_data)))),]

	#z-transform the miRNA data
	for (j in 1:length(rownames(miRNA_data))){
		miRNA_data[j,] <- (as.numeric(miRNA_data[j,]) - mean(as.numeric(miRNA_data[j,])))/sd(as.numeric(miRNA_data[j,]))
	}
	
	#----penalised linear regression----------
	#first we need to subset the data into folds
	new_df <- na.omit(t(rbind(scores,miRNA_data))) #create a combined dataframe
	colnames(new_df) <- c('scores',rownames(miRNA_data))
	folds <- 10 #10 fold cross-validation

	nrows_combined_df <- 1:dim(new_df)[1]

	best_overall_error <- 99999999 #we are optimising based on this
	
	for (i in 0:(folds-1)){
		new_df_subset <- as.data.frame(new_df[!(nrows_combined_df%%folds==i),]) #takes out the 1/nth row of the data set
		
		#train the univariate model, and use these significant predictors as inputs to the penalized model

		linear_models_miRNA <- matrix(NA,nrow=length(rownames(miRNA_data)),ncol=1)
		row.names(linear_models_miRNA) <- rownames(miRNA_data)
		for (j in 1:length(rownames(miRNA_data))){
			univariate_data <- as.data.frame(cbind(new_df_subset[,1],new_df_subset[,(j+1)]))

			colnames(univariate_data) <- c('sig_score','miRNA')
			univariate_model <- lm(formula = sig_score ~ miRNA,data = univariate_data)
	 		linear_models_miRNA[j] <- (summary(univariate_model)$coefficients)[2,4]
		}

		#significant miRNAs are those with p < 0.2:
		significant_miRNAs <- rownames(linear_models_miRNA)[which(linear_models_miRNA < 0.2 & !is.nan(linear_models_miRNA))]
		
		#penalised linear regression
		
		lambda_2_values <- c(0, 0.01, 0.1,1,10,100) #these are the range of lambda parameters we will try in the regression, as we optimise L1
		max_likelihood <- -9999999999 #want to maximise this

		for (lambda2_val in lambda_2_values){
			cross_val_model <- optL1(response = new_df_subset[,1],penalized = new_df_subset[,significant_miRNAs], lambda2 = lambda2_val,data=as.data.frame(new_df_subset),model="linear",fold=10,trace=F)#,trace=F,maxiter=1000,tol=.Machine$double.eps^0.23)
			if ((cross_val_model$fullfit)@loglik > max_likelihood){
				best_model <<- cross_val_model
				best_lambda <- lambda2_val
			}
		}

		#now that we know the best model, let's test it on the other 1/n of the data, and record the error
		unused_df <- as.data.frame(new_df[(nrows_combined_df%%folds==i),])
		current_predictions <- predict(best_model$fullfit, penalized=unused_df[,significant_miRNAs],data=unused_df)

		cur_error <- norm((as.numeric(unused_df[,1]) - as.numeric(current_predictions)),type="2")

		if (cur_error < best_overall_error){
			best_overall_error <- cur_error
			best_overall_model <- best_model
			best_overall_lambda <- best_lambda
		}
	}

	miRNA_names_reported <- intersect(names(coef(best_overall_model$fullfit)), rownames(miRNA_data))

	#return the coefficients
	coef(best_overall_model$fullfit)[miRNA_names_reported]
}

#-----this is for generating the random datasets-----------

all_cancer_types <- c('type1','type2','type3')
all_signature_names <- c('sig1','sig2')
all_signature_genes <- list()

num_samples <- 100
num_mRNA_genes <- 1000
num_miRNA_genes <- 50
sig_length <- 20

for(sig_name in all_signature_names){
	all_signature_genes[[sig_name]]<- paste0('mRNA',sample(x = 1:num_mRNA_genes,size = sig_length,replace=F))
}

all_mRNA_datasets <- list()
all_miRNA_datasets <- list()
all_miRNA <- c() #this is to store all miRNA that are reported in all datasets
for(cancer_type in all_cancer_types){
	all_mRNA_datasets[[cancer_type]] <- replicate(num_samples, rnorm(num_mRNA_genes)) 
	colnames(all_mRNA_datasets[[cancer_type]]) <- paste0('sample',1:num_samples)
	rownames(all_mRNA_datasets[[cancer_type]]) <- paste0('mRNA',1:num_mRNA_genes)

	all_miRNA_datasets[[cancer_type]] <- replicate(num_samples, rnorm(num_miRNA_genes)) 
	colnames(all_miRNA_datasets[[cancer_type]]) <- paste0('sample',1:num_samples)
	rownames(all_miRNA_datasets[[cancer_type]]) <- paste0('miR',1:num_miRNA_genes)
	all_miRNA <- unique(c(all_miRNA,rownames(all_miRNA_datasets[[cancer_type]])))
}

#-----------------------------------------------------------------

all_coeffs <- list(); #this will store all the miRNA coefficients we identify in the model

all_rank_product_matrices <- list(); #this will store all the rank product matrices

for(sig_name in all_signature_names){

	all_coeffs[[sig_name]] <- matrix(0,nrow=length(all_miRNA),ncol=length(all_cancer_types))
	
	row.names(all_coeffs[[sig_name]]) <- all_miRNA
	colnames(all_coeffs[[sig_name]]) <- all_cancer_types

	for (cancer_type in all_cancer_types){

		genes_present <- intersect(rownames(all_mRNA_datasets[[cancer_type]]),all_signature_genes[[sig_name]])
		
		#compute and score the scores
		scores <- apply(all_mRNA_datasets[[cancer_type]][genes_present,], 2, function(x) median(x,na.rm=T))

		#cross-validated linear model
		coeffs <- get_coefficients_pre_filter(cancer_type,scores)

		#store the miRNA results
		all_coeffs[[sig_name]][names(coeffs),cancer_type] <- coeffs
	}
}

#-----next, ocne we have all of the results from the penalised linear models, we will aggregate these with the rank product across cancer types
rank_prod_tables <- list();
RP_out_values <- list();

for (sig_name in all_signature_names){
	all_coeffs[[sig_name]] <- all_coeffs[[sig_name]][which(rowSums(all_coeffs[[sig_name]]==0) < length(colnames(all_coeffs[[sig_name]]))),] # we only want the negative coefficients

	RP.out <- RP(all_coeffs[[sig_name]],rep(1,length(all_cancer_types)))

	RP_out_values[[sig_name]] <- RP.out
	
	rank_prod_tables[[sig_name]] <- topGene(RP.out,cutoff = 0.05,method="pfp",gene.names=rownames(all_coeffs[[sig_name]]))
}

#then save the outputs
save(file='rank_prod_output_pre_filtered.rda',RP_out_values)
save(file='rank_prod_tables_out_pre_filtered.rda',rank_prod_tables)
