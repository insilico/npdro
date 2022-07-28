# =========================================================================#
#' consensus_nestedCV
#
# nested cross-validation for feature-selection/parameter-tuning
# Dawkins/McKinney based on cncv Parvandeh/McKinney, April 2018
#
#############################################
##### Consensus nested cross validation #####
#--------------------------------------------

#' Consensus nested cross validation for feature selection and parameter tuning
#' @param train.ds A training data frame with last column as outcome
#' @param validation.ds A validation data frame with last column as outcome
#' @param label A character vector of the outcome variable column name.
#' @param method.model Column name of outcome variable (string), classification or regression. If the analysis goal is classification make the column a factor type.
#' For regression, make outcome column numeric type.
#' @param is.simulated A TRUE or FALSE character for data type
#' @param ncv_folds A numeric vector to indicate nested cv folds: c(k_outer, k_inner)
#' @param param.tune A TRUE or FALSE character for tuning parameters
#' @param learning_method Name of the method: glmnet/xgbTree/rf
#' @param importance.algorithm A character vestor containing a specific importance algorithm subtype
#' @param wrapper feature selection algorithm including: rf, glmnet, t.test, centrality methods (PageRank, Katz,
#' EpistasisRank, and EpistasisKatz from Rinbix packages), ReliefF family, and etc.
#' @param inner_selection_percent = Percentage of features to be selected in each inner fold.
#' @param inner_selection_positivescores A TRUE or FALSE character to select positive scores (if the value is False, use the percentage method).
#' @param tune.inner_selection_percent A sequence vector of possible percentages for tuning
#' @param tune.k A sequence vector to tune k nearest neighbors in relief method, if TRUE the default grid is seq(1,kmax), where kmax=floor((m-1)/2)
#' and m is the number of samples. However, this kmax is for balanced data. If data are imbalance, where m_minority + m_majority = m, then kmax = floor(m_minority-1).
#' Default is FALSE.
#' @param tuneGrid A data frame with possible tuning values. The columns are named the same as the tuning parameters.
#' This caret library parameter, for more information refer to  http://topepo.github.io/caret/available-models.html.
#' @param relief.k.method A character of numeric to indicate number of nearest neighbors for relief algorithm.
#' Possible characters are: k_half_sigma (floor((num.samp-1)*0.154)), m6 (floor(num.samp/6)),
#' myopic (floor((num.samp-1)/2)), and m4 (floor(num.samp/4))
#' @param num_tree Number of trees in random forest and xgboost methods
#' @param verbose A flag indicating whether verbose output be sent to stdout
#' @return A list with:
#' \describe{
#'   \item{cv.acc}{Training data accuracy}
#'   \item{Validation}{Validation data accuracy}
#'   \item{Features}{number of variables detected correctly in nested cross validation}
#'   \item{Train_model}{Traing model to use for validation}
#'   \item{Elapsed}{total elapsed time}
#' }
#' num.samples <- 100
#' num.variables <- 100
#' pct.signals <- 0.1
#' label <- "class"
#' sim.data <- createSimulation(num.samples = num.samples,
#'                              num.variables = num.variables,
#'                              pct.signals = pct.signals,
#'                              sim.type = "mainEffect",
#'                              label = label,
#'                              verbose = FALSE)
#' cnCV.results <- consensus_nestedCV(train.ds = sim.data$train,
#'                                    validation.ds = sim.data$holdout,
#'                                    label = label,
#'                                    is.simulated = TRUE,
#'                                    ncv_folds = c(10, 10),
#'                                    param.tune = FALSE,
#'                                    learning_method = "rf",
#'                                    importance.algorithm = "ReliefFbestK",
#'                                    num_tree = 500,
#'                                    verbose = FALSE)
#' @family nestedCV
#' @export
consensus_nestedCV <- function(train.ds = NULL,
                               validation.ds = NULL,
                               label = "class",
                               method.model = "classification",
                               is.simulated = TRUE,
                               ncv_folds = c(10, 10),
                               param.tune = FALSE,
                               learning_method = "rf",
                               xgb.obj = "binary:logistic",
                               importance.algorithm = "ReliefFequalK",
                               wrapper = "relief",
                               inner_selection_percent = NULL,
                               inner_selection_positivescores = TRUE,
                               tune.inner_selection_percent = NULL,
                               tune.k = FALSE,
                               tuneGrid = NULL,
                               relief.k.method = "k_half_sigma",
                               num_tree = 500,
                               covars_vec = NULL,
                               covars.pval.adj = 0.05,
                               verbose = FALSE){
  # Define variables
  tune_params <- NULL; accu_vec <- NULL
  Train_accu <- NULL; Test_accu <- NULL
  relief_atts <- list()
  ptm <- proc.time()
  expr_data <- train.ds[, -ncol(train.ds)]
  num_vars <- ncol(expr_data)
  var_names <- colnames(expr_data)
  if(verbose){cat("\nRunning consensus nested cross-validation...\n")}
  if(verbose){cat("\n Create [",ncv_folds[1],"] outer folds\n")}
  outer_folds <- caret::createFolds(train.ds[, label], ncv_folds[1], list = FALSE)
  for (i in 1:ncv_folds[1]){
    atts <- list()
    if(verbose){cat("\n Create [",ncv_folds[2],"] inner folds of outer fold[",i,"]\n")}
    inner_folds <- caret::createFolds(train.ds[, label][outer_folds!=i], ncv_folds[2], list = TRUE)
    if(verbose){cat("\n Feature Selection...\n")}
    for (j in 1:length(inner_folds)){
      inner_idx <- which(outer_folds!=i)[-inner_folds[[j]]]
      if (tune.k){
        m <- table(train.ds[inner_idx, label])
        if (m[1]==m[2]){
          kmax <- floor((sum(m)-1)/2)
        } else {
          m_minority <- m[which.min(m)]
          kmax <- floor(m_minority-1)
        }
        tuneK <- seq(1,kmax)
      } else if (!tune.k){
        tuneK <- NULL
      } else {
        tuneK <- tune.k
      }

      # if someone specifies a numeric value (integer hopefully), use this value for k.
      # However, make sure it is not larger than floor((num.samp.min-1)/2), where
      # num.samp.min  is the min of the train, holdout.. sample sizes.
      # Or you could test the inequality when you encounter each data split.
      # If someone does exceed the threshold, set k to floor((num.samp.min-1)/2)
      # and writing warning that says
      # "ReliefF k too large. Using maximum."
      if (is.numeric(relief.k.method)) {
        if (relief.k.method > floor((dim(train.ds[inner_idx, ])[1]-1)/2)){
          warning("ReliefF k too large. Using maximum.")
          k <- floor((dim(train.ds[inner_idx, ])[1]-1)/2)
        } else {
          k <- relief.k.method
        }
      } else if (relief.k.method ==  "myopic"){
        k <- floor((dim(train.ds[inner_idx, ])[1]-1)/2)
        # where k may change based on num.samp for train, holdout...
      } else if (relief.k.method ==  "m6") { # default "m6" method
        k <- floor(dim(train.ds[inner_idx, ])[1]/6)
        # where k may change based on num.samp for train, holdout...
      } else if (relief.k.method == "m4") {
        k <- floor(dim(train.ds[inner_idx, ])[1]/4)
      } else {
        k <-  floor((dim(train.ds[inner_idx, ])[1]-1)*0.154)
      }
      if (sum(ncv_folds)>(dim(train.ds[inner_idx, ])[1])/3){
        stop("There are less than three observations in each fold")
      }

      k_inner_accs <- NULL
      if (wrapper == "relief" && !is.null(tuneK)){
        for (tk in tuneK) {
          ranked_vars <- CORElearn::attrEval(label, train.ds[inner_idx, ],
                                             estimator = importance.algorithm,
                                             costMatrix = NULL,
                                             outputNumericSplits=FALSE,
                                             kNearestEqual = tk)
          top_vars <- names(which(sort(ranked_vars, decreasing = TRUE)>0))

          # random forest
          inner_trn.data <- as.matrix(train.ds[inner_idx, top_vars])
          inner_tst.data <- as.matrix(train.ds[-inner_idx, top_vars])
          if(method.model == "classification"){
            inner_trn.pheno <- as.factor(train.ds[, label][inner_idx])
            inner_tst.pheno <- as.factor(train.ds[, label][-inner_idx])
          } else {
            inner_trn.pheno <- train.ds[, label][inner_idx]
            inner_tst.pheno <- train.ds[, label][-inner_idx]
          }
          tune_rf.model <- randomForest::randomForest(inner_trn.data,
                                                      y = if(method.model == "classification"){as.factor(inner_trn.pheno)}else{inner_trn.pheno},
                                                      mtry = if (method.model == "classification"){
                                                        max(floor(ncol(inner_trn.pheno)/3), 1)
                                                      } else {
                                                        floor(sqrt(ncol(inner_trn.pheno)))
                                                      },
                                                      ntree = num_tree)
          inner_Train_accu <- ifelse(method.model == "classification", 1 - mean(tune_rf.model$confusion[, "class.error"]),
                                     stats::cor(as.numeric(as.vector(tune_rf.model$predicted)), inner_trn.pheno)^2)
          # test
          inner_test.pred <- stats::predict(tune_rf.model, newdata = inner_tst.data)
          inner_Test_accu <- ifelse(method.model == "classification", confusionMatrix(as.factor(inner_test.pred), as.factor(inner_tst.pheno))$byClass["Balanced Accuracy"],
                                    stats::cor(as.numeric(as.vector(inner_test.pred)), inner_tst.pheno)^2)
          k_inner_accs <- rbind(k_inner_accs, data.frame(train=inner_Train_accu, test=inner_Test_accu))

        }
        # optimal k
        k_inner_accs$trade_off <- abs(k_inner_accs$train-k_inner_accs$test)
        k_inner_accs <- k_inner_accs[order(k_inner_accs$trade_off), ]
        ktuned <- as.integer(rownames(k_inner_accs)[1])
      }
      inner_accs <- NULL
      if (wrapper == "relief" && !is.null(tune.inner_selection_percent)){
        for (step in tune.inner_selection_percent){
          ranked_vars <- CORElearn::attrEval(label, train.ds[inner_idx, ],
                                             estimator = importance.algorithm,
                                             costMatrix = NULL,
                                             outputNumericSplits=FALSE,
                                             kNearestEqual = ifelse(!is.null(tuneK), ktuned, k))

          # evaluations
          wrapper.topN <- step*length(ranked_vars)
          top_vars <- names(sort(ranked_vars, decreasing = TRUE)[1:wrapper.topN])
          # random forest
          inner_trn.data <- as.matrix(train.ds[inner_idx, top_vars])
          inner_tst.data <- as.matrix(train.ds[-inner_idx, top_vars])
          if(method.model == "classification"){
            inner_trn.pheno <- as.factor(train.ds[, label][inner_idx])
            inner_tst.pheno <- as.factor(train.ds[, label][-inner_idx])
          } else {
            inner_trn.pheno <- train.ds[, label][inner_idx]
            inner_tst.pheno <- train.ds[, label][-inner_idx]
          }
          tune_rf.model <- randomForest::randomForest(inner_trn.data,
                                                      y = if(method.model == "classification"){as.factor(inner_trn.pheno)}else{inner_trn.pheno},
                                                      mtry = if (method.model == "classification"){
                                                        max(floor(ncol(inner_trn.pheno)/3), 1)
                                                      } else {
                                                        floor(sqrt(ncol(inner_trn.pheno)))
                                                      },
                                                      ntree = num_tree)
          inner_Train_accu <- ifelse(method.model == "classification", 1 - mean(tune_rf.model$confusion[, "class.error"]),
                                     stats::cor(as.numeric(as.vector(tune_rf.model$predicted)), inner_trn.pheno)^2)
          # test
          inner_test.pred <- stats::predict(tune_rf.model, newdata = inner_tst.data)
          inner_Test_accu <- ifelse(method.model == "classification", confusionMatrix(as.factor(inner_test.pred), as.factor(inner_tst.pheno))$byClass["Balanced Accuracy"],
                                    stats::cor(as.numeric(as.vector(inner_test.pred)), inner_tst.pheno)^2)
          inner_accs <- rbind(inner_accs, data.frame(train=inner_Train_accu, test=inner_Test_accu))
        }
      } else if (wrapper == "relief" && is.null(tune.inner_selection_percent)) {
        ranked_vars <- CORElearn::attrEval(label, train.ds[inner_idx, ],
                                           estimator = importance.algorithm,
                                           costMatrix = NULL,
                                           outputNumericSplits=FALSE,
                                           kNearestEqual = ifelse(!is.null(tuneK), ktuned, k))
      } else if (wrapper == "npdr" && importance.algorithm == "predk"){

        min.group.size <- min(as.numeric(table(train.ds[inner_idx, label])))
        #print(min.group.size)

        my.k <- npdr::knnSURF(2*min.group.size - 1, 0.5)

        #print(my.k)

        # npdr - predk
        npdr.predk <- npdr(label, train.ds[inner_idx, ], regression.type="binomial",
                           attr.diff.type="numeric-abs",
                           nbd.method="relieff",
                           #nbd.method="multisurf",
                           nbd.metric = "manhattan", msurf.sd.frac=.5, knn=my.k,
                           neighbor.sampling="none", separate.hitmiss.nbds=F,
                           dopar.nn = F, dopar.reg=F, padj.method="bonferroni", verbose=F)

        if(inner_selection_positivescores==T){
          ranked_vars <- npdr.predk %>% filter(pval.adj<.8) %>% pull(att)
          ranked_vars <- npdr.predk %>% filter(beta.Z.att>0) %>% pull(att)
        }else{
          ranked_vars <- as.character(npdr.predk$att)
        }

        #print(ranked_vars)

      } else if (wrapper == "npdr" && importance.algorithm == "vwok"){

        min.group.size <- min(as.numeric(table(train.ds[inner_idx, label])))

        vwok.grid <- unique(floor(seq(knnSURF(2*min.group.size-1, 1), knnSURF(2*min.group.size - 1, 0.5), length.out=6)))
        npdr.vwok <- vwok(dats=train.ds[inner_idx, ],
                          k.grid=vwok.grid, # default
                          verbose=F,
                          #separate.hitmiss.nbds=T,
                          #signal.names=signal.names,
                          label=label)

        if(inner_selection_positivescores==T){
          npdr.vwok.positives <- npdr.vwok$vwok.out %>% filter(pval.att<.8) %>% pull(att)
          ranked_vars <- npdr.vwok$vwok.out %>% filter(betas>0) %>% pull(att)
        }else{
          ranked_vars <- as.character(npdr.vwok$vwok.out[,"att"])
        }

      } else if (wrapper == "npdr" && importance.algorithm == "kpca"){

        min.group.size <- min(as.numeric(table(train.ds[inner_idx, label])))

        vwok.grid <- unique(floor(seq(knnSURF(2*min.group.size-1, 1), knnSURF(2*min.group.size - 1, 0.5), length.out=6)))
        npdr.kpca <- kPCA(dats=train.ds[inner_idx, ],
                          num.pcs=2,
                          plot=F,
                          #signal.names=signal.names,
                          scale.data=T,
                          k.grid=vwok.grid,
                          scoring.method="npdr")

        if(inner_selection_positivescores==T){
          npdr.kpca.positives <- npdr.kpca$out.mat %>% filter(pval<.8) %>% pull(att)
          npdr.kpca.positives <- npdr.kpca$out.mat %>% filter(beta>0) %>% pull(att)
          ranked_vars <- npdr.kpca.positives
        }else{
          ranked_vars <- as.character(npdr.kpca$out.mat[,"att"])
        }

      } else if (wrapper == "rf") {
        #rf_model <- CORElearn::CoreModel(label, train.ds[inner_idx, ], model = "rf")
        #ranked_vars <- CORElearn::rfAttrEval(rf_model)
        ranfor.fit <- randomForest(as.factor(class) ~ ., data = train.ds[inner_idx, ], ntree = num_tree)
        rf.importance <- importance(ranfor.fit)
        rf.sorted<-sort(rf.importance, decreasing=T, index.return=T)
        rf.df <-data.frame(att=rownames(rf.importance)[rf.sorted$ix],rf.scores=rf.sorted$x)
        ranked_vars <- as.character(rf.df[,"att"])
      } else if (wrapper == "glmnet") {
        if(method.model == "classification"){
          Class <- as.factor(train.ds[inner_idx, ncol(train.ds)])
        }
        glm_data <- data.frame(train.ds[inner_idx, -ncol(train.ds)], Class)
        ranked_glmnet <- Rinbix::rankGlmnet(glm_data)
        ranked_vars <- ranked_glmnet$score
        names(ranked_vars) <- ranked_glmnet$variable
      } else if (wrapper == "ttest") {
        t_test_pvals <- vector(mode = "numeric", length = num_vars)
        names(t_test_pvals) <- var_names
        for (var_idx in 1:num_vars) {
          t_test_pvals[var_idx] <-  t.test(train.ds[inner_idx, var_idx] ~ train.ds[inner_idx, ncol(train.ds)])$p.value
        }
        ranked_vars <- t_test_pvals
      } else if (wrapper == "PageRank") {
        Adj_mat <- ifelse(cor(train.ds[inner_idx, -ncol(train.ds)]) > 0, 1, 0)
        diag(Adj_mat) <- 0
        ranked_vars <- Rinbix::PageRank(Adj_mat)[, 1]
      } else if (wrapper == "Katz") {
        Adj_mat <- ifelse(cor(train.ds[inner_idx, -ncol(train.ds)]) > 0, 1, 0)
        a <- eigen(Adj_mat)
        beta <- rep(1, nrow(Adj_mat))/nrow(Adj_mat)
        alpha <- 1/max(a$values) - (1/max(a$values))/100
        ranked_vars <- Rinbix::EpistasisKatz(Adj_mat, alpha, beta)
        names(ranked_vars) <- colnames(Adj_mat)
      } else if (wrapper == "EpistasisKatz") {
        if(method.model == "classification"){
          Class <- as.factor(train.ds[inner_idx, ncol(train.ds)])
        }
        regain_data <- data.frame(train.ds[inner_idx, -ncol(train.ds)], Class)
        regain_mat <- Rinbix::regainParallel(regain_data, regressionFamily = ifelse(method.model == "classification",
                                                                                    "binomial", "gaussian"))
        alpha <- 1/mean(colSums(regain_mat))
        beta  <- diag(regain_mat)
        diag(regain_mat) <- 0
        ranked_vars <- Rinbix::EpistasisKatz(regain_mat, alpha, beta)
        names(ranked_vars) <- colnames(regain_mat)
      } else if (wrapper == "EpistasisRank") {
        if(method.model == "classification"){
          Class <- as.factor(train.ds[inner_idx, ncol(train.ds)])
        }
        regain_data <- data.frame(train.ds[inner_idx, -ncol(train.ds)], Class)
        regain_mat <- Rinbix::regainParallel(regain_data, regressionFamily = ifelse(method.model == "classification",
                                                                                    "binomial", "gaussian"))
        er_rank <- Rinbix::EpistasisRank(regain_mat, Gamma_vec = .85)
        ranked_vars <- er_rank$ER
        names(ranked_vars) <- er_rank$gene
      }

      if (wrapper == "relief" && inner_selection_positivescores){
        top_vars <- names(which(sort(ranked_vars, decreasing = TRUE)>0))
      } else if (wrapper == "relief" && !is.null(inner_selection_percent)){
        wrapper.topN <- inner_selection_percent*length(ranked_vars)
        top_vars <- names(sort(ranked_vars, decreasing = TRUE)[1:wrapper.topN])
      } else if (wrapper == "relief" && !is.null(tune.inner_selection_percent)){
        # find optimal percentage
        inner_accs$trade_off <- abs(inner_accs$train-inner_accs$test)
        inner_accs <- inner_accs[order(inner_accs$trade_off), ]
        inner_opt_idx <- as.integer(rownames(inner_accs)[1])
        # selected top optimal percentage genes
        wrapper.topN <- tune.inner_selection_percent[inner_opt_idx]*length(ranked_vars)
        top_vars <- names(sort(ranked_vars, decreasing = TRUE)[1:wrapper.topN])
      } else if (wrapper == "rf"){
        if(inner_selection_positivescores==F){
          wrapper.topN <- round(inner_selection_percent*length(ranked_vars))
          top_vars <- ranked_vars[1:wrapper.topN]
        }else{
          top_vars <- names(ranked_vars.num[which(ranked_vars.num>0)])
        }
      } else if (wrapper == "npdr"){
        if(inner_selection_positivescores==F){
          wrapper.topN <- round(inner_selection_percent*length(ranked_vars))
          top_vars <- ranked_vars[1:wrapper.topN]
        }else{
          top_vars <- ranked_vars
        }
      } else {
        num_ranked_vars <- length(ranked_vars)
        if (num_ranked_vars < wrapper.topN) {
          cat("WARNING glmnet selected less than specified top N:", wrapper.topN)
          cat(" setting top N to length glnmnet selection:", num_ranked_vars, "\n")
          wrapper.topN <- num_ranked_vars
        }
        if (num_ranked_vars < 1) {
          cat("No variable is selected:", num_ranked_vars, "\n")
        } else {
          if (wrapper == "ttest") {
            top_vars <- names(sort(ranked_vars, decreasing = FALSE)[1:wrapper.topN])
          } else {
            top_vars <- names(sort(ranked_vars, decreasing = TRUE)[1:wrapper.topN])
          }
        }
      }
      if(!is.null(covars_vec)){
        cat("Adjusting for covariates...\n")
        train.data_fltr <- train.ds[inner_idx, c(top_vars, label)]
        inner_covars_vec <- covars_vec[inner_idx]
        npdr.results <- npdr('class', train.data_fltr,
                             regression.type='glm', attr.diff.type='numeric-abs',
                             nbd.method='multisurf', nbd.metric = 'manhattan',
                             covars=inner_covars_vec,  # works with sex.covar.mat as well
                             covar.diff.type='match-mismatch', # for categorical covar like sex
                             msurf.sd.frac=0.5, padj.method='bonferroni')

        top_vars <- npdr.results[npdr.results$pval.adj<covars.pval.adj, "att"]
      }
      atts[[j]] <- top_vars
    }
    relief_atts[[i]] <- Reduce(intersect, atts)

    if(param.tune){
      outer_idx <- which(outer_folds!=i)
      trn.data <- as.matrix(train.ds[outer_idx, relief_atts[[i]]])
      tst.data <- as.matrix(train.ds[-outer_idx, relief_atts[[i]]])
      if(method.model == "classification"){
        trn.pheno <- as.factor(train.ds[, label][outer_idx])
        tst.pheno <- as.factor(train.ds[, label][-outer_idx])
      } else {
        trn.pheno <- train.ds[, label][outer_idx]
        tst.pheno <- train.ds[, label][-outer_idx]
      }
      if(verbose){cat("\n Parameter tuning...\n")}
      train_model <- caret::train(x = trn.data,
                                  y = trn.pheno,
                                  method = learning_method,
                                  metric = ifelse(is.factor(trn.pheno), "Accuracy", "RMSE"),
                                  trControl = caret::trainControl(method = "cv"),
                                  tuneGrid = tuneGrid)
      train_pred <- stats::predict(train_model, trn.data)
      train_acc <- ifelse(method.model == "classification",
                          confusionMatrix(train_pred, trn.pheno)$byClass["Balanced Accuracy"],
                          stats::cor(trn.pheno, train_pred)^2)

      test_pred <- stats::predict(train_model, tst.data)
      test_acc <- ifelse(method.model == "classification",
                         confusionMatrix(test_pred, tst.pheno)$byClass["Balanced Accuracy"],
                         stats::cor(tst.pheno, test_pred)^2)

      accu <- abs(train_acc-test_acc)
      tune_params <- rbind(tune_params, data.frame(train_model$bestTune, accu))
    }
  }
  nCV_atts <- Reduce(intersect, relief_atts)
  #print(nCV_atts)
  if(identical(nCV_atts, character(0))){stop("There is no consensus feature!")}
  if(param.tune) {
    tuneParam <- tune_params[which.min(tune_params$accu), ]
  }
  if(verbose){cat("\n Validating...\n")}
  train.data <- as.matrix(train.ds[, nCV_atts])
  if(!is.null(validation.ds)) {
    test.data  <- as.matrix(validation.ds[, nCV_atts])
  }
  if(method.model == "classification"){
    train.pheno <- as.integer(train.ds[, label])
    if(!is.null(validation.ds)) {
      test.pheno <- as.integer(validation.ds[, label])
    }
  } else if (method.model == "regression"){
    train.pheno <- train.ds[, label]
    if(!is.null(validation.ds)) {
      test.pheno <- validation.ds[, label]
    }
  }

  if(verbose){cat("\n Perform [",learning_method,"]\n")}
  if(learning_method == "glmnet"){
    # glmnet - train
    if(param.tune){alpha = tuneParam$alpha; lambda = tuneParam$lambda
    train.model <- glmnet::glmnet(train.data, train.pheno, family = ifelse(method.model == "classification", "binomial", "gaussian"),
                                  alpha = alpha, lambda = lambda)
    train.pred <- stats::predict(train.model, train.data, s = lambda, type = ifelse(method.model == "classification","class", "response"))
    }else{
      train.model <- glmnet::glmnet(train.data, train.pheno, family = ifelse(method.model == "classification", "binomial", "gaussian"))
      train.pred <- stats::predict(train.model, train.data, s = min(train.model$lambda), type = ifelse(method.model == "classification","class", "response"))
    }
    Train_accu <- ifelse(method.model == "classification" , confusionMatrix(as.factor(train.pred),
                                                                            as.factor(train.pheno))$byClass["Balanced Accuracy"],
                         stats::cor(as.numeric(train.pred), train.pheno)^2)
    # glmnet - test
    if(is.null(validation.ds)){
      Test_accu <- NA
    } else {
      if(param.tune){
        test.pred <- stats::predict(train.model, test.data, s = lambda, type = "class")
      } else {
        test.pred <- stats::predict(train.model, test.data, s = min(train.model$lambda), type = "class")
      }
      Test_accu <- ifelse(method.model == "classification", confusionMatrix(as.factor(test.pred),
                                                                            as.factor(test.pheno))$byClass["Balanced Accuracy"],
                          stats::cor(test.pred, test.pheno)^2)
    }
  } else if(learning_method == "xgbTree"){
    # xgboost - train
    dtrain <- xgboost::xgb.DMatrix(data = train.data, label = ifelse(train.pheno == 1, 0, 1))
    dtest <- xgboost::xgb.DMatrix(data = test.data, label = ifelse(test.pheno == 1, 0, 1))
    if(param.tune){
      shrinkage = tuneParam$eta; max_depth = tuneParam$max_depth; gamma = tuneParam$gamma
      subsample = tuneParam$subsample; colsample_bytree = tuneParam$colsample_bytree
      min_child_weight = tuneParam$min_child_weight
    } else {
      shrinkage = 0.3; max_depth = 2; gamma = 0; subsample = 1; colsample_bytree = 1;  min_child_weight = 1
    }
    train.model <- xgboost::xgboost(data = dtrain,
                                    eta = shrinkage,
                                    nrounds = 2,
                                    max_depth = max_depth,
                                    gamma = gamma,
                                    subsample = subsample,
                                    colsample_bytree = colsample_bytree,
                                    min_child_weight = min_child_weight,
                                    objective = ifelse(method.model == "classification", "binary:logistic", "reg:linear"))
    train.pred.prob <- stats::predict(train.model, dtrain)
    train.pred.bin <- as.numeric(train.pred.prob > 0.5)
    Train_accu <- ifelse(method.model == "classification", confusionMatrix(as.factor(ifelse(train.pred.bin == 0, 1, 2)),
                                                                           as.factor(train.pheno))$byClass["Balanced Accuracy"],
                         stats::cor(train.pred.prob, train.pheno)^2)
    # xgboost - test
    if(is.null(validation.ds)){
      Test_accu <- NA
    } else {
      test.pred.prob <- stats::predict(train.model, dtest)
      test.pred.bin <- as.numeric(test.pred.prob > 0.5)
      Test_accu <- ifelse(method.model == "classification", confusionMatrix(as.factor(ifelse(test.pred.bin == 0, 1, 2)),
                                                                            as.factor(test.pheno))$byClass["Balanced Accuracy"],
                          stats::cor(test.pred.prob, test.pheno)^2)
    }
  } else if(learning_method == "rf") {
    # random forest - train
    train.model <- randomForest::randomForest(train.data,
                                              y = if(method.model == "classification"){as.factor(train.pheno)}else{train.pheno},
                                              mtry = if(param.tune && tuneParam$mtry > 1 && tuneParam$mtry < ncol(train.data))
                                              {tuneParam$mtry} else if (method.model == "classification"){
                                                max(floor(ncol(train.data)/3), 1)
                                              } else {
                                                floor(sqrt(ncol(train.data)))
                                              },
                                              ntree = num_tree)
    Train_accu <- ifelse(method.model == "classification", 1 - mean(train.model$confusion[, "class.error"]),
                         stats::cor(as.numeric(as.vector(train.model$predicted)), train.pheno)^2)

    # random forest - test
    if(is.null(validation.ds)){
      Test_accu <- NA
    } else {
      test.pred <- stats::predict(train.model, newdata = test.data)
      Test_accu <- ifelse(method.model == "classification", confusionMatrix(as.factor(test.pred), as.factor(test.pheno))$byClass["Balanced Accuracy"],
                          stats::cor(as.numeric(as.vector(test.pred)), test.pheno)^2)
    }
  }
  elapsed.time <- (proc.time() - ptm)[3]
  if(verbose){cat("nestedCV elapsed time", elapsed.time, "\n")}
  list(cv.acc = Train_accu, Validation = Test_accu, Features = nCV_atts, Train_model = train.model, Elapsed = elapsed.time)
}
# End consensus_nestedCV =====================#
