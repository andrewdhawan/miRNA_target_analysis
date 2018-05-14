#Andrew Dhawan
#Nov 2016
#Contains the code to generate the RR plots for miRNA regulation of mRNA


#-----this is for generating the random datasets-----------

cancer_types <- c('type1','type2','type3')

num_samples <- 100
num_mRNA_genes <- 1000
num_miRNA_genes <- 50

all_mRNA_datasets <- list()
all_miRNA_datasets <- list()
all_miRNA <- c() #this is to store all miRNA that are reported in all datasets
all_mRNA <- c() #this is to store all miRNA that are reported in all datasets

for(cancer_type in all_cancer_types){
	all_mRNA_datasets[[cancer_type]] <- replicate(num_samples, rnorm(num_mRNA_genes)) 
	colnames(all_mRNA_datasets[[cancer_type]]) <- paste0('sample',1:num_samples)
	rownames(all_mRNA_datasets[[cancer_type]]) <- paste0('mRNA',1:num_mRNA_genes)
	all_mRNA <- unique(c(all_mRNA,rownames(all_mRNA_datasets[[cancer_type]])))

	all_miRNA_datasets[[cancer_type]] <- replicate(num_samples, rnorm(num_miRNA_genes)) 
	colnames(all_miRNA_datasets[[cancer_type]]) <- paste0('sample',1:num_samples)
	rownames(all_miRNA_datasets[[cancer_type]]) <- paste0('miR',1:num_miRNA_genes)
	all_miRNA <- unique(c(all_miRNA,rownames(all_miRNA_datasets[[cancer_type]])))
}

#we also need in this case, the matrix of how the miRNA are predicted to regulate the mRNA
#to do this, let's define an adjacency matrix called full_pred_matrix for the graph of these interactions
#in this case, it will be random, but in real life, we fill it with a matrix from target prediction algorithms

full_pred_matrix <- as.matrix(1 * (replicate(length(all_mRNA), rnorm(length(all_miRNA))) > 0))
rownames(full_pred_matrix) <- all_miRNA
colnames(full_pred_matrix) <- all_mRNA

#-----------------------------------------------------------------

for (cancer_type in cancer_types){
	dir.create(cancer_type) #make the output directory
	mRNA_expr_matrix <- all_mRNA_datasets[[cancer_type]]
	miRNA_matrix <- all_miRNA_datasets[[cancer_type]]
	#take common samples
	common_colNames <- intersect(colnames(mRNA_expr_matrix),colnames(miRNA_matrix))

	#take just the common pieces
	mRNA_expr_matrix <- mRNA_expr_matrix[,common_colNames]
	miRNA_matrix <- miRNA_matrix[,common_colNames]

	#before we check the correlation, we will filter out the mRNAs which are barely expressed among the samples, since they will give errors in correlation
	expression_threshold <- 0.10 # means that at least 10% of samples must have a nonzero value of the mRNA
	mRNA_expr_matrix <-mRNA_expr_matrix[which((rowSums(mRNA_expr_matrix==0)) < ((1-expression_threshold) * length(colnames(mRNA_expr_matrix)))),]

	#compute the correlations of all miRNA with all mRNA
	full_corr_matrix <- cor(t(mRNA_expr_matrix[,common_colNames]),t(miRNA_matrix[,common_colNames]),use="complete.obs",method="spearman")

	correlation_cutoff <- seq(-1,0,0.05) #this is going to be the x axis of the RR plots

	miRNA_RR_graphs <- matrix(0,nrow=length(colnames(full_corr_matrix)),ncol=length(correlation_cutoff))
	row.names(miRNA_RR_graphs) <- colnames(full_corr_matrix) #these are the names of the miRNAs

	#now we compute the "relative risk" for each miRNA
	for (i in 1:length(colnames(full_corr_matrix))){
		for (j in 1:length(correlation_cutoff)){
			nume <- sum(((full_corr_matrix[,i] <= correlation_cutoff[j]) * 1 )  * ((full_pred_matrix[i,]>0)*1),na.rm=TRUE )
			nume <- nume / sum((full_corr_matrix[,i] <= correlation_cutoff[j])*1,na.rm=TRUE)

			deno <- sum(((full_corr_matrix[,i] >= correlation_cutoff[j]) *1) * ((full_pred_matrix[i,]>0)*1),na.rm=TRUE )
			deno <- deno / sum((full_corr_matrix[,i] >= correlation_cutoff[j])*1,na.rm=TRUE)

			miRNA_RR_graphs[i,j] <- (nume/deno)
		}
	}

	N_random_sampling <- 100 #this is the number of bootstrap resamples we will do
	stored_sampling_results <- array(0,c(N_random_sampling,length(correlation_cutoff),length(colnames(full_corr_matrix))))

	for (k in 1:length(colnames(full_corr_matrix))){
		N_transcripts_to_take <- sum(full_pred_matrix[k,]) #this should be the same as the number of transcripts predicted as targets for the miRNA

		for (i in 1:N_random_sampling){
			samp <- sample(1:length(rownames(full_corr_matrix)),N_transcripts_to_take,replace=TRUE)
			for (j in 1:length(correlation_cutoff)){
				nume <- sum(((full_corr_matrix[samp,k] <= correlation_cutoff[j]) * 1 ) ,na.rm=TRUE)
				nume <- nume / sum((full_corr_matrix[,k] <= correlation_cutoff[j])*1,na.rm=TRUE)

				deno <- sum(((full_corr_matrix[samp,k] >= correlation_cutoff[j]) *1) ,na.rm=TRUE)
				deno <- deno / sum((full_corr_matrix[,k] >= correlation_cutoff[j])*1,na.rm=TRUE)

				stored_sampling_results[i,j,k] <- (nume/deno)
			}
		}
	}

	#store all these values as we will need these for plotting
	sampling_res_mean <- array(0,c(length(correlation_cutoff),length(colnames(full_corr_matrix))))
	sampling_res_sd <- array(0,c(length(correlation_cutoff),length(colnames(full_corr_matrix))))
	sampling_res_error_high <- array(0,c(length(correlation_cutoff),length(colnames(full_corr_matrix))))
	sampling_res_error_low <- array(0,c(length(correlation_cutoff),length(colnames(full_corr_matrix))))
	
	for (k in 1:length(colnames(full_corr_matrix))){
		for (j in 1:length(correlation_cutoff)){
			sampling_res_mean[j,k] <- mean(stored_sampling_results[,j,k])
			sampling_res_sd[j,k] <- sd(stored_sampling_results[,j,k])
			sampling_res_error_high[j,k] <- quantile(na.omit(stored_sampling_results[,j,k]),0.95)
			sampling_res_error_low[j,k] <- quantile(na.omit(stored_sampling_results[,j,k]),0.05)
		}
	}

	save(miRNA_RR_graphs,sampling_res_mean,sampling_res_sd,sampling_res_error_high,sampling_res_error_low,correlation_cutoff,file=paste(cancer_type,'/miRNA_RR_graphs.rda',sep=''))

	#plot all the values below
	for(i in 1:length(rownames(miRNA_matrix))){
		which_miRNA <- i

		start_ind <-  min(which(!is.nan(miRNA_RR_graphs[which_miRNA,])))
		end_ind <- max( which(!is.nan(miRNA_RR_graphs[which_miRNA,])))
		if(is.finite(start_ind + end_ind)){
			error_high <- sampling_res_error_high[start_ind:end_ind,which_miRNA]
			error_low <- sampling_res_error_low[start_ind:end_ind,which_miRNA]
			dev.new()
			plot(correlation_cutoff[start_ind:end_ind],miRNA_RR_graphs[which_miRNA,start_ind:end_ind],ylim=c(0,max(max(miRNA_RR_graphs[which_miRNA,start_ind:end_ind]),max(error_high,na.rm=T))),
				xlab='Correlation with miRNA expression',ylab='RR of anti-correlated targets')
			lines(correlation_cutoff[start_ind:end_ind],miRNA_RR_graphs[which_miRNA,start_ind:end_ind])
			lines(correlation_cutoff[start_ind:end_ind], error_low,lty=2) 
			lines(correlation_cutoff[start_ind:end_ind], error_high,lty=2)
			lines(correlation_cutoff[start_ind:end_ind],sampling_res_mean[start_ind:end_ind,which_miRNA],lwd=2,lty=2)

			title(rownames(miRNA_RR_graphs)[which_miRNA])
			dev.copy(pdf,file=paste0(cancer_type,'/',rownames(miRNA_matrix)[which_miRNA],'_RR_plot.pdf'),width=5,height=5)
			dev.off()
			graphics.off()
		}
	}
}
