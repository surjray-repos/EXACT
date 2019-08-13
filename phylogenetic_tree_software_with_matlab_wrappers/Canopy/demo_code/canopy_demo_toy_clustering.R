#######################################################
#######################################################
#######                                         #######
#######            Toy: try it yourself         #######
#######                                         #######
#######################################################
#######################################################
library(Canopy)
data(toy)
projectname = 'toy'
R = toy$R; X = toy$X; WM = toy$WM; Wm = toy$Wm
epsilonM = toy$epsilonM; epsilonm = toy$epsilonm; Y = toy$Y


## HERE WE OVERIDE THIS DATA WITH SOME DATA READ FROM SOME FILE
R <- read.delim("R_file.tsv",header=T);
R <- as.matrix(as.data.frame(R))
X <- read.delim("X_file.tsv",header=T);
X <- as.matrix(as.data.frame(X))
#####################################################
#We can use the WM and Wm provided by Canopy since our method is CNA agnostic
WM <- read.delim("Wm_file_Case_Conflict.tsv",header=T);
WM <- as.matrix(as.data.frame(WM))
Wm <- read.delim("Wm_file.tsv",header=T);
Wm <- as.matrix(as.data.frame(Wm))
#####################################################
Y <- read.delim("Y_file.tsv",header=T);
Y <- as.matrix(as.data.frame(Y))
#####################################################

#######################################################
#######################################################
#######                                         #######
#######               SNA clustering            #######
#######                                         #######
#######################################################
#######################################################
cluster_parameters <- read.delim("cluster_parameters_file.tsv", header=F);
cluster_parameters <- as.matrix(as.data.frame(cluster_parameters))

num_cluster = cluster_parameters[[1]]:cluster_parameters[[2]];
#num_cluster=2:5 # Range of number of clusters to run
num_run=10 # How many EM runs per clustering step for each mutation cluster wave
canopy.cluster=canopy.cluster(R = R,
                              X = X,
                              num_cluster = num_cluster,
                              num_run = num_run)

# BIC to determine the optimal number of mutation clusters
bic_output=canopy.cluster$bic_output
plot(num_cluster,bic_output,xlab='Number of mutation clusters',ylab='BIC',type='b',main='BIC for model selection')
abline(v=num_cluster[which.max(bic_output)],lty=2)
# Visualization of clustering result
Mu=canopy.cluster$Mu # VAF centroid for each cluster
Tau=canopy.cluster$Tau  # Prior for mutation cluster, with a K+1 component
sna_cluster=canopy.cluster$sna_cluster # cluster identity for each mutation
colc=c('green4','red3','royalblue1','darkorange1','royalblue4',
       'mediumvioletred','seagreen4','olivedrab4','steelblue4','lavenderblush4')
pchc=c(17,0,1,15,3,16,4,8,2,16)
plot((R/X)[,1],(R/X)[,2],xlab='Sample1 VAF',ylab='Sample2 VAF',col=colc[sna_cluster],pch=pchc[sna_cluster],ylim=c(0,max(R/X)),xlim=c(0,max(R/X)))
library(scatterplot3d)
scatterplot3d((R/X)[,1],(R/X)[,2],(R/X)[,3],xlim=c(0,max(R/X)),ylim=c(0,max(R/X)),zlim=c(0,max(R/X)),color=colc[sna_cluster],pch=pchc[sna_cluster],
              xlab='Sample1 VAF',ylab='Sample2 VAF',zlab='Sample3 VAF')
dev.copy(png, 'toy3_clustering.png')
dev.off()
write.table(canopy.cluster$sna_cluster, "output_sna_cluster_file.tsv", sep="\t") 

#######################################################
#######################################################
#######                                         #######
#######               MCMC sampling             #######
#######                                         #######
#######################################################
#######################################################
parameters <- read.delim("parameters_file.tsv",header=F);
parameters <- as.matrix(as.data.frame(parameters))

K = parameters[[3]]:parameters[[4]]; numchain = parameters[[5]]
sampchain = canopy.sample.cluster(R = R, X = X, WM = WM, Wm = Wm, epsilonM = epsilonM, 
                          epsilonm = epsilonm, C = NULL, sna_cluster = sna_cluster, Y = Y, K = K, 
                          numchain = numchain, max.simrun = parameters[[6]],
                          min.simrun = parameters[[7]], writeskip = parameters[[8]],
                          projectname = projectname, cell.line = FALSE,
                          plot.likelihood = TRUE)
#sampchain = canopy.sample.cluster.nocna(R = R, X = X, sna_cluster = sna_cluster, K = K, 
#                          numchain = numchain, max.simrun = parameters[[6]],
#                          min.simrun = parameters[[7]], writeskip = parameters[[8]],
#                          projectname = projectname, cell.line = FALSE,
#                          plot.likelihood = TRUE)
save.image(file = paste(projectname, '_postmcmc_image.rda',sep=''),
           compress = 'xz')


#######################################################
#######################################################
#######                                         #######
#######   BIC to determine number of subclones  #######
#######                                         #######
#######################################################
#######################################################
#library(Canopy)
#data(toy)
projectname='toy'
load(paste(projectname, '_postmcmc_image.rda', sep=''))
burnin = parameters[[1]]
thin = parameters[[2]]
# If pdf = TRUE, a pdf will be generated.
bic = canopy.BIC(sampchain = sampchain, projectname = projectname, K = K,
                 numchain = numchain, burnin = burnin, thin = thin, pdf = TRUE)
optK = K[which.max(bic)]


#######################################################
#######################################################
#######                                         #######
#######         posterior tree evaluation       #######
#######                                         #######
#######################################################
#######################################################
post = canopy.post(sampchain = sampchain, projectname = projectname, K = K,
                   numchain = numchain, burnin = burnin, thin = thin, 
                   optK = optK, post.config.cutoff = 0.05)
samptreethin = post[[1]]   # list of all post-burnin and thinning trees
samptreethin.lik = post[[2]]   # likelihoods of trees in samptree
config = post[[3]]
config.summary = post[[4]]
print(config.summary)
# first column: tree configuration
# second column: posterior configuration probability in the entire tree space
# third column: posterior configuration likelihood in the subtree space
# note: if modes of posterior probabilities aren't obvious, run sampling longer.


#######################################################
#######################################################
#######                                         #######
#######          Tree output and plot           #######
#######                                         #######
#######################################################
#######################################################
# choose the configuration with the highest posterior likelihood
config.i = config.summary[which.max(config.summary[,3]),1]
cat('Configuration', config.i, 'has the highest posterior likelihood.\n')
output.tree = canopy.output(post, config.i, C=NULL)
pdf.name = paste(projectname, '_config_highest_likelihood.pdf', sep='')
txt.name = paste(projectname, '_config_highest_likelihood_clusters.txt', sep='')
canopy.plottree(output.tree, pdf = TRUE, pdf.name = pdf.name)
canopy.plottree(output.tree, txt = TRUE, txt.name = txt.name)
canopy.plottree(output.tree, pdf = FALSE)

#Read the output.tree object in R and write output.tree$P to get M probability values
#Read output.tree$Z to get the tree matrix

write.table(output.tree$P, "output_P_file.tsv", sep="\t") 
write.table(output.tree$Z, "output_Z_file.tsv", sep="\t") 
