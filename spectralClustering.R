####################################################
######EDA Visualization & Clustering Functions######
####################################################

###############################
######Spectral Clustering######
###############################

library(Rcpp)
library(inline)
library(parallel)
library(foreach)
library(doSNOW)


####Compute the Graph Laplacian Matrix L###
###########################################

###Similarity Matrix###
#######################

###Pure R Implementation###
###########################
#L2 norm between two vectors
l2_norm <- function(xi,xj){
  ans <- sum((xi - xj)^2)
  return(ans)
}

#L2 norm LAPACK version
l2_normLAPACK <- function(xi,xj){
  ans <- norm(as.matrix(xi-xj),type="F")
  return(ans)
}

#Gaussian RBF kernel (xi and xj are rows of a matrix/dataframe)
rbf_kern <- function(xi,xj,sigma){
  ans <- exp(-l2_norm(xi=xi,xj=xj) / (2*(sigma^2)))
  return(ans)
}

#Wrapper for Gaussian RBF kernel to compute the rows of a matrix/dataframe
rbf_wrap <- function(i,j,data,sigma){
  xi <- data[i,]
  xj <- data[j,]
  ans <- rbf_kern(xi=xi,xj=xj,sigma=sigma)
  ans <- as.vector(unname(ans))
  return(ans)
}

#Vectorize the function to be compatible with outer for fast processing
rbf_wrapVec <- Vectorize(rbf_wrap,vec=c("i","j"))

#Similarity Graph S
sim_matrix <- function(data,sigma=1){
  n <- nrow(data)
  ans <- outer(1:n,1:n,rbf_wrapVec,data=data,sigma=sigma)
  return(ans)
}


###C++ Implementation###
########################
#Rcpp L2 norm between two vectors
cppFunction("double l2_normC(NumericVector xi, NumericVector xj){
  double ans;
  ans = sqrt(sum(pow((xi - xj),2)));
  return(ans);
}")

#Rcpp L2 norm on low syntactic sugar diet
cppFunction("double l2_normNoSugar(NumericVector xi, NumericVector xj){
  int n = xi.size();
  double norm_int;
  for(int i=0; i <= n; ++i){
    norm_int += pow((xi[i] - xj[i]),2);
  }
  double ans;
  ans = pow(norm_int,.5);
  return(ans);
}")

#Rcpp Gaussian RBF kernel
cppFunction("double rbf_kernC(NumericVector xi, NumericVector xj, double sigma) {
  double norm;
  norm = sum(pow((xi - xj),2));
  double ans;
  ans = norm / (2*pow(sigma,2));
  return(exp(-(ans)));
}")

#Rcpp Gaussian RBF kernel on low syntactic sugar diet
cppFunction("double rbf_kernNoSugar(NumericVector xi, NumericVector xj, double sigma){
  int n = xi.size();
  double norm;
  for(int i=0; i < n; ++i){
    norm +=  pow((xi[i] - xj[i]),2);
  }
  double ans;
  ans = norm / (2*pow(sigma,2));
  return(exp(-(ans)));
}")

#Rcpp Wrapper for Gaussian RBF kernel to compute the rows of a matrix/dataframe
rbf_wrapC <- function(i,j,data,sigma){
  xi <- as.numeric(data[i,])
  xj <- as.numeric(data[j,])
  ans <- rbf_kernNoSugar(xi=xi,xj=xj,sigma=sigma)
  return(ans)
}

#Rcpp Vectorize the function to be compatible with outer for fast processing
rbf_wrapVecC <- Vectorize(rbf_wrapC,vec=c("i","j"))

#Rcpp Similarity Graph S
sim_matrixC <- function(data,sigma=1){
  n <- nrow(data)
  ans <- outer(1:n,1:n,rbf_wrapVecC,data=data,sigma=sigma)
  return(ans)
}

###100% C++
cppFunction("NumericMatrix gaussian_kern(NumericMatrix data, NumericVector gridA, NumericVector gridB, double sigma){
            int n_row = pow(gridA.size(),0.5);
            int n_iter = gridA.size();
            NumericMatrix output(n_row,n_row);
            for(int k=0; k<n_iter; ++k){
                int i = gridA[k];
                int j = gridB[k];
                NumericVector xi = data.row(i);
                NumericVector xj = data.row(j);
                double norm;
                norm = sum(pow((xi - xj),2));
                double ans;
                ans = norm / (2*pow(sigma,2));
                output(i,j) = exp(-(ans));
            }
            return(output);
}")

sim_matrixCpp <- function(data,sigma=1){
  grid <- expand.grid((1:nrow(data))-1,(1:nrow(data))-1)
  mat_out <- gaussian_kern(data=data,gridA=grid$Var1,gridB=grid$Var2,sigma=sigma)
  return(mat_out)
}









#benchmarking
#2 vectors
x1 <- 1:1e1
y1 <- x1 + rnorm(n=length(x1),mean=0,sd=.25)
x2 <- 1:1e2
y2 <- x2 + rnorm(n=length(x2),mean=0,sd=.25)
x3 <- 1:1e3
y3 <- x3 + rnorm(n=length(x3),mean=0,sd=.25)
x4 <- 1:1e4
y4 <- x4 + rnorm(n=length(x4),mean=0,sd=.25)
x5 <- 1:1e5
y5 <- x5 + rnorm(n=length(x5),mean=0,sd=.25)
x6 <- 1:1e6
y6 <- x6 + rnorm(n=length(x6),mean=0,sd=.25)

kern1 <- microbenchmark(rbf_kern(x1,y1,1),rbf_kernC(x1,y1,1),rbf_kernNoSugar(x1,y1,1))
kern2 <- microbenchmark(rbf_kern(x2,y2,1),rbf_kernC(x2,y2,1),rbf_kernNoSugar(x2,y2,1))
kern3 <- microbenchmark(rbf_kern(x3,y3,1),rbf_kernC(x3,y3,1),rbf_kernNoSugar(x3,y3,1))
kern4 <- microbenchmark(rbf_kern(x4,y4,1),rbf_kernC(x4,y4,1),rbf_kernNoSugar(x4,y4,1))
kern5 <- microbenchmark(rbf_kern(x5,y5,1),rbf_kernC(x5,y5,1),rbf_kernNoSugar(x5,y5,1))
kern6 <- microbenchmark(rbf_kern(x6,y6,1),rbf_kernC(x6,y6,1),rbf_kernNoSugar(x6,y6,1))

kernBench_df <- data.frame(time=kern1$time[order(kern1$expr)],call=rep(c("PureR","SugarC","PureC"),each=100),length=rep(1,300))
kernBench_df <- rbind(kernBench_df,data.frame(time=kern2$time[order(kern2$expr)],call=rep(c("PureR","SugarC","PureC"),each=100),length=rep(2,300)))
kernBench_df <- rbind(kernBench_df,data.frame(time=kern3$time[order(kern3$expr)],call=rep(c("PureR","SugarC","PureC"),each=100),length=rep(3,300)))
kernBench_df <- rbind(kernBench_df,data.frame(time=kern4$time[order(kern4$expr)],call=rep(c("PureR","SugarC","PureC"),each=100),length=rep(4,300)))
kernBench_df <- rbind(kernBench_df,data.frame(time=kern5$time[order(kern5$expr)],call=rep(c("PureR","SugarC","PureC"),each=100),length=rep(5,300)))
kernBench_df <- rbind(kernBench_df,data.frame(time=kern6$time[order(kern6$expr)],call=rep(c("PureR","SugarC","PureC"),each=100),length=rep(6,300)))

ggplot(data=kernBench_df) +
  geom_boxplot(aes(x=call,y=time,fill=call),position="dodge") +
  facet_grid(.~length,labeller=labeller(length=c("1"="1e1","2"="1e2","3"="1e3","4"="1e4","5"="1e5","6"="1e6"))) +
  guides(fill=FALSE) +
  labs(y="Time (microseconds)") +
  theme_bw() +
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=13.5),axis.text.x=element_text(vjust=.75,angle=45,size=11.5)) +
  scale_y_log10()

#matrix
mat5 <- as.matrix(rbind(sapply(rep(seq(1,5,by=1),each=1),function(x) {rnorm(n=5,mean=x,sd=.25)})))
mat25 <- as.matrix(rbind(sapply(rep(seq(1,25,by=1),each=1),function(x) {rnorm(n=25,mean=x,sd=.25)})))
mat50 <- as.matrix(rbind(sapply(rep(seq(1,50,by=1),each=1),function(x) {rnorm(n=50,mean=x,sd=.25)})))
mat100 <- as.matrix(rbind(sapply(rep(seq(1,100,by=1),each=1),function(x) {rnorm(n=100,mean=x,sd=.25)})))

mat5bench <- microbenchmark(sim_matrix(mat5,1),sim_matrixC(mat5,1),sim_matrixCpp(mat5,1))
mat25bench <- microbenchmark(sim_matrix(mat25,1),sim_matrixC(mat25,1),sim_matrixCpp(mat25,1))
mat50bench <- microbenchmark(sim_matrix(mat50,1),sim_matrixC(mat50,1),sim_matrixCpp(mat50,1))
mat100bench <- microbenchmark(sim_matrix(mat100,1),sim_matrixC(mat100,1),sim_matrixCpp(mat100,1))

kernMatBench_df <- data.frame(time=mat5bench$time[order(mat5bench$expr)],call=rep(c("PureR","MixedC","PureC"),each=100),size=rep(5,300))
kernMatBench_df <- rbind(kernMatBench_df,data.frame(time=mat25bench$time[order(mat25bench$expr)],call=rep(c("PureR","MixedC","PureC"),each=100),size=rep(25,300)))
kernMatBench_df <- rbind(kernMatBench_df,data.frame(time=mat50bench$time[order(mat50bench$expr)],call=rep(c("PureR","MixedC","PureC"),each=100),size=rep(50,300)))
kernMatBench_df <- rbind(kernMatBench_df,data.frame(time=mat100bench$time[order(mat100bench$expr)],call=rep(c("PureR","MixedC","PureC"),each=100),size=rep(100,300)))

ggplot(data=kernMatBench_df) +
  geom_boxplot(aes(x=call,y=time,fill=call),position="dodge") +
  facet_grid(.~size) +
  guides(fill=FALSE) +
  labs(y="Time (microseconds)") +
  theme_bw() +
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=13.5),axis.text.x=element_text(vjust=.75,angle=45,size=11.5)) +
  scale_y_log10()



#Next need to compute the Affinity Matrix A
#make affinity matrix with KNN filter applied to S
aff_matrix <- function(sim_mat,neighbor=10){
  n <- dim(sim_mat)[1]
  ans <- matrix(0,ncol=n,nrow=n)
  for(i in 1:n){
    n_sim <- sort(sim_mat[i,],decreasing=TRUE)[1:neighbor]
    n_index <- which(sim_mat[i,] %in% n_sim)
    ans[i,n_index] <- sim_mat[i,n_index]
    ans[n_index,i] <- sim_mat[i,n_index]
  }
  return(ans)
}

#Next need to compute the Degree Matrix D
#make degree matrix, a diagonal matrix where element is each vertex's degree of connectedness
deg_matrixR <- function(aff_mat){
  d_mat <- matrix(0,ncol=ncol(aff_mat),nrow=nrow(aff_mat))
  d_degree <- rowSums(aff_mat)
  diag(d_mat) <- d_degree
  return(d_mat)
}

#Compute the Graph Laplacian Matrix
graph_laplacian <- function(degreeM,affinityM,type=NULL){
  if(is.null(type)){
    stop("You need to tell me what Graph Laplacian you want! (unnorm, simple, norm, general)")
  }
  if(!(type %in% c("unnorm","simple","norm","general"))){
    stop("You need to tell me what Graph Laplacian you want! (unnorm, simple, norm, general)")
  }
  if(type=="unnorm"){
    l_mat <- degreeM - affinityM
    message("Giving you the unnormalized Laplacian!")
  }
  if(type=="simple"){
    l_mat <- diag(nrow(degreeM)) - solve(degreeM) %*% affinityM
    message("Giving you the simple (random walk) Laplacian!")
  }
  if(type=="norm"){
    l_unnorm <- degreeM - affinityM
    "%^%" <- function(M, power){with(eigen(M), vectors %*% (values^power * solve(vectors)))}
    l_mat <- (degreeM %^% (-1/2)) %*% l_unnorm %*% (degreeM %^% (-1/2))
    message("Giving you the normalized (symmetric) Laplacian!")
  }
  if(type=="general"){
    l_unnorm <- degreeM - affinityM
    l_mat <- solve(degreeM) %*% l_unnorm
    message("Giving you the generalized Laplacian!")
  }
  return(l_mat)
}



  

#Now we want to preform spectral decomposition of the Graph Laplacian
spec_decomp <- function(L,k){
  
  #spectral decomposition of the Graph Laplacian
  spec_decomp <- eigen(L,symmetric=TRUE)
  eigenval <- spec_decomp$values
  n_eigen <- length(eigenval)
  
  #plot 15 smallest eigenvalues
  eigenval_vis <- data.frame(cbind(eigenval,index=1:length(eigenval)))
  require(ggplot2)
  eig_plot <- ggplot(data=eigenval_vis[(n-15):n,],aes(x=index,y=eigenval)) +
    geom_point(colour="chartreuse3",size=4) +
    scale_x_reverse() + 
    labs(x="Index",y="Eigenvalue (15 smallest)") +
    theme_bw()
  print(eig_plot)
  
  #spectral clustering on data projected onto space spanned bb k smallest eigenvectors
  k_eigenvec <- spec_decomp$vector[,(n_eigen-k+1):n_eigen]
  eigenvec_vis <- data.frame(cbind(e1=k_eigenvec[,3],e2=k_eigenvec[,2],e3=k_eigenvec[,1],index=1:n_eigen))
  #plot against the 3rd (y axis) and 2nd (x axis) smallest eigenvectors
  eigenspace_plot <- ggplot(data=eigenvec_vis,aes(x=e2,y=e3,colour=y_dat)) +
    geom_point(size=4) +
    scale_color_brewer(palette="Paired") +
    guides(colour=FALSE) +
    labs(x="2nd Smallest Eigenvector",y="1st Smallest Eigenvector") +
    theme_bw()
  print(eigenspace_plot)
  
  #plot data projected onto each eigenvector
  eigenvec_plot1 <- ggplot(data=eigenvec_vis) +
    geom_point(aes(x=index,y=e1,colour=y_dat)) +
    scale_color_brewer(palette="Paired") +
    guides(colour=FALSE) +
    labs(x="Index",y="3rd Smallest Eigenvector") +
    theme_bw()
  eigenvec_plot2 <- ggplot(data=eigenvec_vis) +
    geom_point(aes(x=index,y=e2,colour=y_dat)) +
    scale_color_brewer(palette="Paired") +
    guides(colour=FALSE) +
    labs(x="Index",y="2nd Smallest Eigenvector") +
    theme_bw()
  eigenvec_plot3 <- ggplot(data=eigenvec_vis) +
    geom_point(aes(x=index,y=e3,colour=y_dat)) +
    scale_color_brewer(palette="Paired") +
    guides(colour=FALSE) +
    labs(x="Index",y="1st Smallest Eigenvector") +
    theme_bw()
}



#spectral decomposition of the Laplacian
p1_eigen <- eigen(p1_lmat,symmetric=TRUE)
p1_eigenval <- p1_eigen$values

#plot 15 smallest eigenvalues
p1_plot2dat <- data.frame(cbind(eig_vec=p1_eigenval,index=1:length(p1_eigenval)))
p1_plot2 <- ggplot(data=p1_plot2dat[435:450,],aes(x=index,y=eig_vec)) + geom_point(colour="chartreuse3",size=4) + labs(x="Number",y="Eigenvalue",title="15 Smallest Eigenvalues") + scale_x_reverse() + theme_bw()

#spectral clustering with k=3 smallest eigenvectors 
p1_v <- p1_eigen$vectors[,448:450]
#plot spectral clustering
p1_plot3dat <- data.frame(cbind(e_1=p1_v[,2],e_2=p1_v[,1]))
p1_plot3 <- ggplot(data=p1_plot3dat,aes(x=e_1,y=e_2)) + geom_point(colour=p1_color,size=4) + labs(x="Second smallest eigenvector",y="Third smallest eigenvector",title="Spectral Clustering (k=10)") + theme_bw()

#plot eigenvectors
p1_plot4dat <- data.frame(cbind(p1_plot3dat,index=1:450))
p4_plot4a <- ggplot(data=p1_plot4dat) + geom_point(aes(x=index,y=e_1),colour=p1_color) + labs(y="2nd smallest eigenvector",x="Index") + theme_bw()
p4_plot4b <- ggplot(data=p1_plot4dat) + geom_point(aes(x=index,y=e_2),colour=p1_color) + labs(y="3rd smallest eigenvector",x="Index") + theme_bw()
grid.arrange(p4_plot4a,p4_plot4b,ncol=2)

#k-means on transformed data
p1_kmean <- kmeans(p1_v[,2:1],centers=3)
p1_clus_col <- NULL
for(i in 1:length(p1_kmean$cluster)){
  if(p1_kmean$cluster[i] == 1){
    p1_clus_col[i] <- "tomato"
  }
  else if(p1_kmean$cluster[i] == 2){
    p1_clus_col[i] <- "steelblue"
  }
  else if(p1_kmean$cluster[i] == 3){
    p1_clus_col[i] <- "chartreuse3"
  }
}
p1_plot5 <- ggplot(data=p1_data,aes(x=X,y=Y)) + geom_point(colour=p1_clus_col) + labs(title="Spectral Clustering (k=10)") + theme_bw()
