# library(igraph)     ###
# library(ggplot2)  #
# library(gplots)   #
# library(corrplot) #
# library(grDevices)
# library(ggcorrplot) #
# library(cowplot)    #
# library(gridExtra)
# library(ggpubr)
# library(reshape2)  ###

get_int_pairs <- function(g, nbias, multiway, diff.cor.vars=NULL, int.partner.list=NULL){

  kvec <- igraph::degree(g)

  # which variables are not connected in network
  idx.not.connected <- which(kvec == 0)

  # which variables are connected in network
  if(length(idx.not.connected) > 0){                              # if any are not connected
    idx.connected <- seq(1,length(kvec),by=1)[-idx.not.connected]
  }else{                                                          # if all are connected
    idx.connected <- seq(1,length(kvec),by=1)
  }

  mycomp <- igraph::components(g,mode="strong")
  idx.max <- which.max(mycomp$csize)

  comp.vertices <- which(as.numeric(mycomp$membership)==idx.max)

  diff.cor.vars <- diff.cor.vars
  int.partner.list <- int.partner.list

  # functional variables are randomly defined from connected variables
  if(is.null(diff.cor.vars)){

    if(length(comp.vertices) >= nbias){                             # if number connected at least nbias

      diff.cor.vars <- sample(comp.vertices,size=nbias,replace=F)

    }else{                                                          # if too few are connected

      n.add.vars <- nbias - length(comp.vertices)
      idx.add <- which(c(idx.connected %in% comp.vertices)==F)
      diff.cor.vars <- c(comp.vertices, sample(idx.add,size=n.add.vars,replace=F)) # add some from unconnected

    }

  }

  ## make make noisy covariance matrix form structure of adjacency matrix
  Adj <- as.matrix(igraph::get.adjacency(g))

  if(!is.null(diff.cor.vars) && is.null(int.partner.list)){

    connected.list <- list()
    int.partner.list <- list()
    for(i in 1:length(diff.cor.vars)){

      connected.list[[i]] <- which(abs(Adj[diff.cor.vars[i],] - 1) < 1e-9)

      int.possible <- intersect(connected.list[[i]],diff.cor.vars[-i])

      if(length(int.possible) >= (multiway-1)){
        mysamp <- sample(int.possible, size=(multiway-1), replace=F)
      }else{
        mysamp <- int.possible
      }
      int.partner.list[[i]] <- mysamp
    }

  }

  list(deg.vec = kvec, A.mat = Adj, sig.vars = diff.cor.vars, interacting.partners=int.partner.list)
}

################################################################################################
## parameters
#
# g                - (igraph) igraph network or something that could be coerced to be an igraph object
# num.variables    - (numeric) number of variables in simulation
# hi.cor.tmp       - (numeric) baseline correlation between variables connected in network
# lo.cor.tmp       - (numeric) baseline correlation between non-connected variables in network
# diff.cor.vars    - (numeric) functional variables - i.e. those with differential correlation
# int.partner.list - (list) interacting partners of functional variables
# plot.graph       - (logical) will make network plot from input igraph object
# make.diff.cors   - (logical) used to generate differential correlation, controlled by setting
#                    different values of hi.cor.tmp for different groups (i.e. cases and controls)
# nbias            - (numeric) number of functional variables
# use.Rcpp         - (logical) ideally this will speed up eigendecomposition, but it actually seems to be slower
# multiway         - (numeric) controls how many pairwise correlations will be different between groups
################################################################################################
generate_structured_corrmat <- function(g=NULL,
                                        num.variables=100,
                                        hi.cor.tmp=0.8,
                                        lo.cor.tmp=0.2,
                                        hi.cor.fixed=0.8,
                                        diff.cor.vars=NULL,
                                        int.partner.list=NULL,
                                        corr.structure=NULL,
                                        plot.graph=F,
                                        make.diff.cors=F,
                                        nbias=1, use.Rcpp=F,
                                        multiway=2){
  if(use.Rcpp){
    Rcpp::cppFunction(depends="RcppArmadillo",
                      'arma::vec getEigenValues(arma::mat M) {
                      return arma::eig_sym(M);
  }')
  }
  if(abs(as.integer(num.variables) - num.variables) > 1e-9){
    stop("generate_structured_corrmat: num.variables should be a positive integer")
  }

  if(is.null(corr.structure)){
    all.stuff.for.corrs <- get_int_pairs(g, nbias, multiway, diff.cor.vars=diff.cor.vars, int.partner.list=int.partner.list)
  }else{
    all.stuff.for.corrs <- corr.structure
  }

  p <- num.variables

  kvec <- all.stuff.for.corrs$deg.vec

  # make plot of random graph
  if(plot.graph==T){
    plot(g,vertex.size=1,vertex.label.color=kvec)
  }

  if(is.null(diff.cor.vars) && is.null(int.partner.list)){

    diff.cor.vars <- all.stuff.for.corrs$sig.vars
    int.partner.list <- all.stuff.for.corrs$interacting.partners

  }

  ## make make noisy covariance matrix form structure of adjacency matrix
  Adj <- all.stuff.for.corrs$A.mat

  # if simply generating a correlation matrix without between-group differential correlation
  if(make.diff.cors==F){

    adj.upper.tri <- Adj[upper.tri(Adj)]                                                         # upper triangle of adjacency
    hi.cors <- (hi.cor.tmp + rnorm(length(adj.upper.tri), sd = 0.1)) * adj.upper.tri             # connected (high) correlations
    adj.upper.tri.adinv <- -(adj.upper.tri - 1)                                                  # additive inverse of adjacency upper triangle
    lo.cors <- (lo.cor.tmp + rnorm(length(adj.upper.tri.adinv), sd = 0.1)) * adj.upper.tri.adinv # non-connected (low) correlations

    upper.mat <- matrix(0, nrow=dim(Adj)[1], ncol=dim(Adj)[1]) # initialize connected correlations matrix
    lower.mat <- matrix(0, nrow=dim(Adj)[1], ncol=dim(Adj)[1]) # initialize non-connected correlations matrix

    upper.mat[upper.tri(upper.mat)] <- hi.cors                 # load upper triangle of connected matrix
    lower.mat[upper.tri(lower.mat)] <- lo.cors                 # load upper triangle of non-connected matrix
    full.mat <- upper.mat + lower.mat                          # add connected/non-connected to get full matrix
    full.mat <- full.mat + t(full.mat)                         # make symmetric
    diag(full.mat) <- 1                                        # 1's on the diagonal

    # if generating correlation matrix and between-group differential correlation
  }else{

    # find connections of functional attributes, store connections in list,
    # and create binary matrix for functional and connections only (Subset of Adj)
    #
    adj.tmp1 <- Adj*0                                                      # initialize binary matrix
    for(i in 1:length(diff.cor.vars)){                                     # for each functional,

      adj.tmp1[diff.cor.vars[i], int.partner.list[[i]]] <- 1 # populate Adjacency for functional vars and connections from graph
      adj.tmp1[int.partner.list[[i]], diff.cor.vars[i]] <- 1 # symmetry

    }

    nonfunc.adj <- Adj - adj.tmp1  # binary matrix for all non-functional connections (Subset of Adj)

    adj.noncon.adinv <- -(Adj - 1) # additive inverse of Adjacency for all non-connections

    adj.func.upper.tri <- adj.tmp1[upper.tri(adj.tmp1)]                                                                                                       # functional connections only (upper triangle)
    func.upper.tri <- (hi.cor.tmp + rnorm(length(adj.func.upper.tri), sd = 0.1)) * adj.func.upper.tri                                                         # functional connections correlations matrix
    noncon.upper.tri <- (lo.cor.tmp + rnorm(length(adj.noncon.adinv[upper.tri(adj.noncon.adinv)]), sd = 0.1)) * adj.noncon.adinv[upper.tri(adj.noncon.adinv)] # non-connections correlations matrix
    nonfunc.upper.tri <- (hi.cor.fixed + rnorm(length(nonfunc.adj[upper.tri(nonfunc.adj)]), sd = 0.1)) * nonfunc.adj[upper.tri(nonfunc.adj)]                  # non-functional connections correlations matrix

    mat1 <- matrix(0, nrow=dim(Adj)[1], ncol=dim(Adj)[1]) # initialize functional correlation matrix
    mat2 <- matrix(0, nrow=dim(Adj)[1], ncol=dim(Adj)[1]) # initialize non-connection correlation matrix
    mat3 <- matrix(0, nrow=dim(Adj)[1], ncol=dim(Adj)[1]) # initialize non-functional connection correlation matrix

    mat1[upper.tri(mat1)] <- func.upper.tri    # load upper triangle of functional matrix
    mat2[upper.tri(mat2)] <- noncon.upper.tri  # load upper triangle of non-connection matrix
    mat3[upper.tri(mat3)] <- nonfunc.upper.tri # load upper triangle of non-functional matrix

    full.mat <- mat1 + mat2 + mat3     # combine all correlations in single matrix
    full.mat <- full.mat + t(full.mat) # make symmetric
    diag(full.mat) <- 1                # 1's on diagonal

  }
  R <- full.mat

  # correct for negative eigenvalues to make matrix positive definite
  #
  if(use.Rcpp){ # compute eigenvalues and make diag matrix
    R.d <- diag(sort(c(getEigenValues(R)),decreasing=T))
  }else{
    R.d <- diag(eigen(R)$values)
  }

  eig.vals <- diag(R.d)            # vector of eigenvalues

  if (any(eig.vals<0)){            # if any eigenvalues are negative

    R.V <- eigen(R)$vectors        # compute eigenvectors,
    eig.vals[eig.vals < 0] <- 1e-7 # make negative into small positive,
    diag(R.d) <- eig.vals          # replace in R.d,
    R.fix <- R.V%*%R.d%*%t(R.V)    # compute new correlation matrix, and
    R <- R.fix                     # store in R

  }
  R <- as.matrix(R)

  # make 1's on diagonal of R
  #
  inv.diag <- 1/diag(R)            # multiplicative inverse of diag(R)
  mydiag <- diag(length(inv.diag)) # initialize diagonal matrix for inv.diag
  diag(mydiag) <- inv.diag         # swap 1's for inv.diag
  mydiag <- sqrt(mydiag)           # take sqrt of diag matrix
  R <- mydiag %*% R %*% mydiag     # compute corrected correlation matrix with 1's on diagonal (Still Pos. Def.)

  # return correlation matrix, degree vector, adjacency matrix, and functional variables
  list(corrmat = R, deg.vec = kvec, A.mat = Adj, sig.vars = diff.cor.vars, interacting.partners=int.partner.list)
}

###plot
ast_fn <- function(pval){

  if(is.na(pval)){
    ast <- ""
  }else if(pval < 1e-6){
    ast <- "******"
  }else if(pval < 1e-5){
    ast <- "*****"
  }else if(pval < 1e-4){
    ast <- "****"
  }else if(pval < 1e-3){
    ast <- "***"
  }else if(pval < 1e-2){
    ast <- "**"
  }else if(pval < .05){
    ast <- "*"
  }else{
    ast <- ""
  }

}

###plot?
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

###plot?g
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

#' make_sim_plot
#' @export
make_sim_plot <- function(case.sim, ctrl.sim, n1=30, n2=30, save.plot=F){

  case.sub.mat <- case.sim$corrmat
  colnames(case.sub.mat) <- paste("ROI",1:ncol(case.sub.mat))
  row.names(case.sub.mat) <- paste("ROI",1:ncol(case.sub.mat))

  ctrl.sub.mat <- ctrl.sim$corrmat
  colnames(ctrl.sub.mat) <- paste("ROI",1:ncol(ctrl.sub.mat))
  row.names(ctrl.sub.mat) <- paste("ROI",1:ncol(ctrl.sub.mat))

  p1 <- ggcorrplot(case.sub.mat, hc.order = FALSE, type = "lower", colors=c("darkblue","white","darkred"),
                   lab = TRUE, show.diag=T, title="Case Correlation Matrix")+theme(plot.title=element_text(size=14, hjust=0.5))
  p2 <- ggcorrplot(ctrl.sub.mat, hc.order = FALSE, type = "lower", colors=c("darkblue","white","darkred"),
                   lab = TRUE, show.diag=T, title="Control Correlation Matrix")+theme(plot.title=element_text(size=14, hjust=0.5))

  #z_ij_1 <- 0.5 * log((abs((1 + r_ij_1) / (1 - r_ij_1))))
  #z_ij_2 <- 0.5 * log((abs((1 + r_ij_2) / (1 - r_ij_2))))
  #Z_ij <- abs(z_ij_1 - z_ij_2) / sqrt((1/(n1 - 3) + 1 / (n2 - 3)))
  #pval <- 2 * pnorm(-abs(Z_ij))

  n1 <- n1
  n2 <- n2
  Z.case <- 0.5 * log((abs(1 + case.sub.mat) / (1 - case.sub.mat)))
  Z.ctrl <- 0.5 * log((abs(1 + ctrl.sub.mat) / (1 - ctrl.sub.mat)))
  Z <- abs(Z.case - Z.ctrl) / sqrt((1/(n1-3) + 1/(n2 - 3)))
  pvals <- 2 * pnorm(-abs(Z))
  diag(pvals) <- 1

  diag(Z) <- 0

  melted_zmat <- melt(Z, na.rm=T)
  melted_pvals <- melt(pvals, na.rm=T)

  my.ast <- character()
  for(i in 1:length(melted_pvals[,"value"])){

    my.ast[i] <- ast_fn(melted_pvals[i, "value"])

  }
  melted_zmat <- data.frame(melted_zmat, stars=my.ast)
  melted_zmat$Var1 <- factor(melted_zmat$Var1, levels=colnames(Z))

  # Reorder the correlation matrix
  upper_tri_zmat <- get_upper_tri(Z)

  # Melt the correlation matrix
  melted_zmat <- melt(upper_tri_zmat, na.rm = TRUE)

  upper_tri_pmat <- get_upper_tri(pvals)
  melted_pvals <- melt(upper_tri_pmat, na.rm = TRUE)

  amat <- case.sim$A.mat
  colnames(amat) <- colnames(Z)
  row.names(amat) <- colnames(Z)

  upper_tri_amat <- get_upper_tri(amat)
  melted_amat <- melt(upper_tri_amat, na.rm = TRUE)

  my.ast <- character()
  for(i in 1:length(melted_pvals[,"value"])){

    my.ast[i] <- ast_fn(melted_pvals[i, "value"])

  }
  melted_zmat <- data.frame(melted_zmat, stars=my.ast)
  melted_zmat$Var1 <- factor(melted_zmat$Var1, levels=colnames(Z))

  p3 <- ggplot(data=melted_zmat, aes(Var2, Var1, fill=value))+
    geom_tile(color="white")+
    geom_text(aes(Var2, Var1, label=stars), color="#ffeda0", size=6, nudge_y=-0.1)+
    scale_fill_gradient2(low="#fee0d2", high="darkred", mid="#fee0d2", midpoint=0, limit=c(0, max(Z)), space="Lab",
                         name="z-score")+
    theme_minimal()+
    ggtitle("Differential Correlation Test")+
    coord_equal() +
    theme(axis.title=element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.text.y=element_text(size=12),
          plot.title=element_text(size=14, hjust=0.5))

  my.adj <- character()
  sigvars <- paste("ROI", case.sim$sig.vars)
  for(i in 1:length(melted_amat[,"value"])){

    var1 <- as.character(melted_amat[i,1])
    var2 <- as.character(melted_amat[i,2])
    check1 <- sum(var1 %in% sigvars)
    check2 <- sum(var2 %in% sigvars)

    if(check1 == 1){
      idx.var1.sig <- which(sigvars == var1)
      my.pairs.tmp <- paste("ROI", case.sim$interacting.partners[[idx.var1.sig]])

      check3 <- sum(var2 %in% my.pairs.tmp)
      if(check3 == 1){
        my.adj[i] <- "+"
      }else{
        my.adj[i] <- "-"
      }
    }else if(check2 == 1){
      idx.var2.sig <- which(sigvars == var2)
      my.pairs.tmp <- paste("ROI", case.sim$interacting.partners[[idx.var2.sig]])

      check3 <- sum(var2 %in% my.pairs.tmp)
      if(check3 == 1){
        my.adj[i] <- "+"
      }else{
        my.adj[i] <- "-"
      }
    }else{
      my.adj[i] <- "-"
    }

  }

  melted_amat <- data.frame(melted_amat, CorrType=my.adj)
  melted_amat$Var1 <- factor(melted_amat$Var1, levels=colnames(amat))
  colnames(melted_amat)[3] <- "A_ij"
  melted_amat[,3] <- as.factor(as.character(melted_amat[,3]))

  p4 <- ggplot(data=melted_amat, aes(Var2, Var1, fill=A_ij))+
    geom_tile(color="white")+
    #geom_text(aes(Var2, Var1, label=Functional), color="#ffeda0", size=6, nudge_y=0)+
    geom_point(aes(Var2, Var1, shape=CorrType), color="#ffeda0", size=4)+
    scale_fill_manual(values=c("#525252","#a50f15"), breaks=c(0,1), name="A_ij",
                      guide = guide_legend(override.aes = list()))+
    #scale_color_discrete(name = "meh", guide="none") +
    #scale_shape_discrete(name  ="Corr Type",
    #                     breaks=c("+", "-"),
    #                     labels=c("Functional", "Irrelevant"))+
    scale_shape_manual(values=c(45,43), labels=c("Irrelevant","Functional"), breaks=c("-","+"), name="CorrType")+
    #scale_fill_discrete(name = "c", guide = "none")+
    #scale_size_manual(values=c(5,5))+
    theme_minimal()+
    coord_equal() +
    ggtitle("Adjacency Matrix from Graph")+
    #guides(colour = guide_legend(override.aes = list(color = NA)))+
    guides(fill = guide_legend(override.aes = list(shape = c(NA,NA) )))+
    theme(axis.title=element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.text.y=element_text(size=12),
          plot.title=element_text(size=14, hjust=0.5),
          legend.key=element_rect(fill='#525252'))

  p5 <- ggarrange(p1, p2, p3, p4, nrow=2, ncol=2, common.legend=F, align="hv")
  print(p5)

  if(save.plot==T){
    ggsave(plot=p5, "test_simcor_fn.pdf", height=13, width=13)
  }

}

#' make_case_ctrl_corrs
#' @export
make_case_ctrl_corrs <- function(user.adjacency=NULL,      # user-supplied binary adjacency matrix (Ex: could be based on real data)
                                 graph.type="Erdos-Renyi", # type of graph
                                 hi.cor=0.8, lo.cor=-0.8,   # controls the difference in correlations between cases and controls
                                 multiway=2,               # controls the number of pairwise differential correlations between functional variables
                                 num.variables=10,         # number of variables
                                 prob.connected=NULL,      # for Erdos-Renyi (prob of connection)
                                 out.degree=NULL,          # for scale-free (degree of each node)
                                 num.case.mats=1,          # how many case matrices to make
                                 num.ctrl.mats=1,          # how many control matrices to make
                                 pct.signals=0.5){         # fraction of total variables that are functionally related to case-control outcome
#case-control driver function starts here
  nbias <- pct.signals * num.variables

  if(graph.type == "Erdos-Renyi"){

    if(is.null(prob.connected)){

      prob <- 1 / (num.variables + 0.5)

    }else{

      prob <- prob.connected

    }

    g <- igraph::erdos.renyi.game(num.variables, prob)

  }else if(graph.type == "Scalefree"){

    if(is.null(out.degree)){

      g <- igraph::barabasi.game(num.variables, directed=F)

    }else{

      g <- igraph::barabasi.game(num.variables, m=out.degree, directed=F)

    }

  }else if(!is.null(user.graph)){

    g <- graph_from_adjacency_matrix(as.matrix(user.adjacency),
                                     mode = "undirected",
                                     weighted = NULL,
                                     diag = TRUE,
                                     add.colnames = NULL,
                                     add.rownames = NA)

  }

  all.stuff.for.corrs <- get_int_pairs(g, nbias, multiway)
  diff.cor.vars <- all.stuff.for.corrs$sig.vars
  int.partner.list <- all.stuff.for.corrs$interacting.partners

  case.mat.list <- list()
  for(i in 1:num.case.mats){

    # case correlation matrix
    case.sub <- generate_structured_corrmat(g=g,
                                            num.variables=num.variables,
                                            hi.cor.tmp=lo.cor,
                                            lo.cor.tmp=0.2,
                                            hi.cor.fixed=0.8,
                                            #diff.cor.vars=diff.cor.vars,
                                            #int.partner.list=int.partner.list,
                                            corr.structure=all.stuff.for.corrs,
                                            plot.graph=F,
                                            make.diff.cors=T,
                                            nbias=nbias, use.Rcpp=F,
                                            multiway=multiway)
    case.mat.list[[i]] <- case.sub

  }

  ctrl.mat.list <- list()
  for(i in 1:num.ctrl.mats){

    # control correlation matrix
    ctrl.sub <- generate_structured_corrmat(g=g,
                                            num.variables=num.variables,
                                            hi.cor.tmp=hi.cor,
                                            lo.cor.tmp=0.2,
                                            hi.cor.fixed=0.8,
                                            #diff.cor.vars=diff.cor.vars,
                                            #int.partner.list=int.partner.list,
                                            corr.structure=all.stuff.for.corrs,
                                            plot.graph=F,
                                            make.diff.cors=T,
                                            nbias=nbias, use.Rcpp=F,
                                            multiway=multiway)

    ctrl.mat.list[[i]] <- ctrl.sub

  }

  list(cases=case.mat.list, ctrls=ctrl.mat.list)
}

##############################
# functions for pre/post (paired) corr-mat-variable data
##############################

make_paired_corrs <- function(base.sub, cor.val, use.Rcpp=F){

  base.mat.tmp <- base.sub$corrmat
  followup.mat <- base.mat.tmp
  for(j in 1:length(base.sub$sig.vars)){

    n.pairs <- length(base.sub$interacting.partners[[j]])
    if(n.pairs > 0){

      for(k in 1:n.pairs){

        followup.mat[base.sub$sig.vars[j], base.sub$interacting.partners[[j]][k]] <- cor.val
        followup.mat[base.sub$interacting.partners[[j]][k], base.sub$sig.vars[j]] <- cor.val

      }

    }

  }

  noise.mat <- matrix(0, nrow=nrow(base.mat.tmp), ncol=ncol(base.mat.tmp))
  noise.mat[upper.tri(noise.mat)] <- rnorm((nrow(base.mat.tmp) * (nrow(base.mat.tmp) - 1))/2, sd=0.1)
  noise.mat <- noise.mat + t(noise.mat)
  diag(noise.mat) <- 0
  followup.mat <- followup.mat + noise.mat

  R <- followup.mat
  # correct for negative eigenvalues to make matrix positive definite
  #
  if(use.Rcpp){ # compute eigenvalues and make diag matrix
    R.d <- diag(sort(c(getEigenValues(R)),decreasing=T))
  }else{
    R.d <- diag(eigen(R)$values)
  }

  eig.vals <- diag(R.d)            # vector of eigenvalues

  if (any(eig.vals<0)){            # if any eigenvalues are negative

    R.V <- eigen(R)$vectors        # compute eigenvectors,
    eig.vals[eig.vals < 0] <- 1e-7 # make negative into small positive,
    diag(R.d) <- eig.vals          # replace in R.d,
    R.fix <- R.V%*%R.d%*%t(R.V)    # compute new correlation matrix, and
    R <- R.fix                     # store in R

  }
  R <- as.matrix(R)

  # make 1's on diagonal of R
  #
  inv.diag <- 1/diag(R)            # multiplicative inverse of diag(R)
  mydiag <- diag(length(inv.diag)) # initialize diagonal matrix for inv.diag
  diag(mydiag) <- inv.diag         # swap 1's for inv.diag
  mydiag <- sqrt(mydiag)           # take sqrt of diag matrix
  R <- mydiag %*% R %*% mydiag     # compute corrected correlation matrix with 1's on diagonal (Still Pos. Def.)

  return(R)
}

make_paired_corrmats <- function(user.adjacency=NULL,      # user-supplied binary adjacency matrix (Ex: could be based on real data)
                                 graph.type="Erdos-Renyi", # type of graph
                                 hi.cor=0.8, lo.cor=-0.8,   # controls the difference in correlations between cases and controls
                                 multiway=2,               # controls the number of pairwise differential correlations between functional variables
                                 num.variables=10,         # number of variables
                                 prob.connected=NULL,      # for Erdos-Renyi (prob of connection)
                                 out.degree=NULL,          # for scale-free (degree of each node)
                                 num.paired.mats=1,        # how many paired matrices to simulate (i.e. number of individuals)
                                 pct.signals=0.5){         # fraction of total variables that are functionally related to case-control outcome

  nbias <- pct.signals * num.variables

  if(graph.type == "Erdos-Renyi"){

    if(is.null(prob.connected)){

      prob <- 1 / (num.variables + 0.5)

    }else{

      prob <- prob.connected

    }

    g <- igraph::erdos.renyi.game(num.variables, prob)

  }else if(graph.type == "Scalefree"){

    if(is.null(out.degree)){

      g <- igraph::barabasi.game(num.variables, directed=F)

    }else{

      g <- igraph::barabasi.game(num.variables, m=out.degree, directed=F)

    }

  }else if(!is.null(user.graph)){

    g <- graph_from_adjacency_matrix(as.matrix(user.adjacency),
                                     mode = "undirected",
                                     weighted = NULL,
                                     diag = TRUE,
                                     add.colnames = NULL,
                                     add.rownames = NA)

  }

  all.stuff.for.corrs <- get_int_pairs(g, nbias, multiway)
  diff.cor.vars <- all.stuff.for.corrs$sig.vars
  int.partner.list <- all.stuff.for.corrs$interacting.partners

  baseline.mat.list <- list()
  followup.mat.list <- list()
  for(i in 1:num.paired.mats){

    # baseline correlation matrix
    base.sub <- generate_structured_corrmat(g=g,
                                            num.variables=num.variables,
                                            hi.cor.tmp=lo.cor,
                                            lo.cor.tmp=0.2,
                                            hi.cor.fixed=0.8,
                                            #diff.cor.vars=diff.cor.vars,
                                            #int.partner.list=int.partner.list,
                                            corr.structure=all.stuff.for.corrs,
                                            plot.graph=F,
                                            make.diff.cors=T,
                                            nbias=nbias, use.Rcpp=F,
                                            multiway=multiway)
    baseline.mat.list[[i]] <- base.sub

    follow.sub <- make_paired_corrs(base.sub, cor.val=hi.cor)

    followup.mat.list[[i]] <- follow.sub

  }

  list(baseline.mats=baseline.mat.list, followup.mats=followup.mat.list)
}




# the following two functions can be sourced at the end of the case-control and paired sim scripts, respectively
# there is an example after each function is defined that is commented out

################################################################################################################################################################
# for case-control sims
################################################################################################################################################################

#' make_big_case_ctrl_mat
#' Makes a big matrix containing all case and control correlation matrices from make_case_ctrl_corrs() function
#'
#' parameters:
#'
#' case.list - $cases list element of make_case_ctrl_corrs() function (Example: my.corrmats$cases)
#' ctrl.list - $ctrls list element of make_case_ctrl_corrs() function (Example: my.corrmats$ctrls)
#' save.file - (logical) will save .rds file if true
#' file.name - (character) name of file if you would like to save it (do not include .rds extension)
#'
#' output:
#'
#' big.case.ctrl.mat - (data.frame) each row is an ROI pair and each column is a subject with case/control identifier
#' @export
make_big_case_ctrl_mat <- function(case.list, ctrl.list,
                                   save.file=FALSE, file.name=NULL,
                                   rm.lower.tri = F){

  class.vec <- rep(0,length(case.list)+length(ctrl.list))
  big.case.mat <- NULL
  case.names <- character()
  sub.count <- 1
  for(i in 1:length(case.list)){

    case.mat.tmp <- case.list[[i]]$corrmat
    diag(case.mat.tmp) <- NA
    if (rm.lower.tri){
      case.mat.tmp[lower.tri(case.mat.tmp)] <- NA
    }
    melted.case.mat <- melt(case.mat.tmp, na.rm=T)

    case.names[i] <- paste("sub",sub.count, "_case",sep="")
    class.vec[i] <- 1  # controls are already 0's
    sub.count <- sub.count + 1
    if(i == 1){
      corr.pair.names <- apply(melted.case.mat[,c("Var1","Var2")], 1, function(x)paste("ROI", x[1], ".", "ROI", x[2], sep=""))
    }

    big.case.mat <- cbind(big.case.mat, melted.case.mat[,"value"])

  }
  big.case.mat <- as.data.frame(big.case.mat)
  colnames(big.case.mat) <- case.names
  row.names(big.case.mat) <- corr.pair.names

  big.ctrl.mat <- NULL
  ctrl.names <- character()
  for(i in 1:length(ctrl.list)){

    ctrl.mat.tmp <- ctrl.list[[i]]$corrmat
    diag(ctrl.mat.tmp) <- NA
    if (rm.lower.tri){
      ctrl.mat.tmp[lower.tri(ctrl.mat.tmp)] <- NA
    }
    melted.ctrl.mat <- melt(ctrl.mat.tmp, na.rm=T)

    ctrl.names[i] <- paste("sub",sub.count, "_ctrl",sep="")
    sub.count <- sub.count + 1

    big.ctrl.mat <- cbind(big.ctrl.mat, melted.ctrl.mat[,"value"])

  }
  big.ctrl.mat <- as.data.frame(big.ctrl.mat)
  colnames(big.ctrl.mat) <- ctrl.names
  row.names(big.ctrl.mat) <- corr.pair.names

  big.case.ctrl.mat <- data.frame(big.case.mat, big.ctrl.mat)
  wide.case.ctrl.dat <- data.frame(t(big.case.ctrl.mat),class=class.vec)
  output.list <- list(wide.case.ctrl.dat = wide.case.ctrl.dat,
                      case.ctrl.corrmats=big.case.ctrl.mat,
                      degree.vec=case.list[[1]]$deg.vec,
                      adjacency.mat=case.list[[1]]$A.mat,
                      functional.ROIs=case.list[[1]]$sig.vars,
                      diffcor.ROI.partners=case.list[[1]]$interacting.partners)

  if(save.file==T){
    if(!is.null(file.name)){
      saveRDS(output.list, paste(file.name, ".rds"))
    }else{
      saveRDS(output.list, "case_ctrl_corrmats.rds")
    }
  }

  return(output.list)
}

# example of creating large matrix (actually a data.frame), containing all case and control correlation matrices
# you can save the file if you want by setting save.file=TRUE

#case.ctrl.corrmats <- make_big_case_ctrl_mat(case.list=my.corrmats$cases, ctrl.list=my.corrmats$ctrls, save.file=F, file.name="example_case_ctrl_fmri_sim")
#str(case.ctrl.corrmats)

################################################################################################################################################################
# for paired sims
################################################################################################################################################################

# Makes a big matrix containing all pre and post correlation matrices from make_paired_corrmats() function
#
# parameters:
#
# pre.list - $baseline.mats list element of make_paired_corrmats() function (Example: my.corrmats$baseline.mats)
# post.list - $followup.mats list element of make_paired_corrmats() function (Example: my.corrmats$followup.mats)
# save.file - (logical) will save .rds file if true
# file.name - (character) name of file if you would like to save it (do not include .rds extension)
#
# output:
#
# big.paired.mat - (data.frame) each row is an ROI pair and each column is a subject with pre/post identifier
make_big_pre_post_mat <- function(pre.list, post.list,
                                  save.file=FALSE, file.name=NULL,
                                  rm.lower.tri=F){

  big.pre.mat <- NULL
  for(i in 1:length(pre.list)){

    pre.mat.tmp <- pre.list[[i]]$corrmat
    diag(pre.mat.tmp) <- NA
    if (rm.lower.tri){
      pre.mat.tmp[lower.tri(pre.mat.tmp)] <- NA
    }
    melted.pre.mat <- melt(pre.mat.tmp, na.rm=T)

    if(i == 1){
      corr.pair.names <- apply(melted.pre.mat[,c("Var1","Var2")], 1, function(x)paste("ROI", x[1], ".", "ROI", x[2], sep=""))
    }

    big.pre.mat <- cbind(big.pre.mat, melted.pre.mat[,"value"])

  }
  big.pre.mat <- as.data.frame(big.pre.mat)
  row.names(big.pre.mat) <- corr.pair.names

  big.post.mat <- NULL
  for(i in 1:length(post.list)){

    post.mat.tmp <- post.list[[i]]
    diag(post.mat.tmp) <- NA
    if (rm.lower.tri){
      post.mat.tmp[lower.tri(post.mat.tmp)] <- NA
    }
    melted.post.mat <- melt(post.mat.tmp, na.rm=T)

    big.post.mat <- cbind(big.post.mat, melted.post.mat[,"value"])

  }
  big.post.mat <- as.data.frame(big.post.mat)
  row.names(big.post.mat) <- corr.pair.names

  big.paired.mat <- data.frame(big.pre.mat, big.post.mat)

  # add columns in mixed effect format
  # subject column (repeated because of treatment pre/post)
  # treatment columns pre/post
  nsubs <- ncol(big.pre.mat)  # input rows are correlations, cols are pre/post subj
  subj_col <- c(1:nsubs, 1:nsubs)
  treat_col <- c(rep("pre", nsubs), rep("post",nsubs))
  paired.mixformat.dat <- data.frame(subjects=subj_col,
                                     treatment=treat_col,
                                     t(big.paired.mat))

  nsubs <- length(pre.list)
  colnames(big.paired.mat) <- c(paste("sub", 1:nsubs, "_pre", sep=""), paste("sub", 1:nsubs, "_post", sep=""))

  output.list <- list(paired.mixformat.dat=paired.mixformat.dat,
                      paired.corrmats=big.paired.mat,
                      degree.vec=pre.list[[1]]$deg.vec,
                      adjacency.mat=pre.list[[1]]$A.mat,
                      functional.ROIs=pre.list[[1]]$sig.vars,
                      diffcor.ROI.partners=pre.list[[1]]$interacting.partners)

  if(save.file==T){
    if(!is.null(file.name)){
      saveRDS(output.list, paste(file.name, ".rds"))
    }else{
      saveRDS(output.list, "paired_corrmats.rds")
    }
  }

  return(output.list)
}

# example of creating large matrix (actually a data.frame), containing all pre and post correlation matrices
# you can save the file if you want by setting save.file=TRUE

#paired.corrmats <- make_big_pre_post_mat(pre.list=my.corrmats$baseline.mats, post.list=my.corrmats$followup.mats, save.file=F, file.name="example_paired_fmri_sim")
#str(paired.corrmats)





