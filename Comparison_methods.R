
pred_outs = list()
for (iter_num in 1:length(test_start)) {
  #options(warn = -1)
  # j=3
  # divide in train/test
  training = dff %>% filter(week < test_start[iter_num])
  testing = dff %>% filter(week >= test_start[iter_num] &
                             week < test_start[iter_num] + 70) # up to 10 weeks
  mmat = model.matrix(formula(reg_f), data = training)
  if (max(mmat[, 2]) == 0) {
    mmat = mmat[, -2]
    reg_f <-
      "covid_status ~ I(log(death+1)) + lagged_death + I(log(population)) + linear_trend + quadratic_trend + sine_trend + cosine_trend"
  }
  mmat = model.matrix(formula(reg_f),data = training)
  tmat = model.matrix(formula(reg_f), data = testing)
  
  # run GLM probit and make predictions
  begintime = Sys.time()
  glm_train <-
    glm(
      formula = formula(reg_f),
      family = binomial(link = "probit"),
      data = training
    )
  glm_pred <-
    predict.glm(glm_train, newdata = testing, type = "response")
  endtime <- Sys.time()
  glm_time <- as.numeric(difftime(endtime, begintime, units = "mins"))
  
  # run SAR probit and make predictions
  begintime = Sys.time()
  spdf <- training
  coordinates(spdf) <- ~ long + lat
  neib <- knn2nb(knearneigh(coordinates(spdf), longlat = TRUE))
  lw <- nb2listw(neib, style = "W")
  weightmat <- as(as_dgRMatrix_listw(lw), "CsparseMatrix")
  probit_sar_train <-
    ProbitSpatialFit(
      formula = formula(reg_f),
      data = training,
      W = weightmat,
      DGP = "SAR"
    )
  
  fulldf = rbind(training, testing)
  preddf <- fulldf
  coordinates(preddf) <- ~ long + lat
  predneib <- knn2nb(knearneigh(coordinates(preddf), longlat = TRUE))
  predlw <- nb2listw(predneib, style = "W")
  predweightmat <- as(as_dgRMatrix_listw(predlw), "CsparseMatrix")
  fullmat = model.matrix(formula(reg_f), data = fulldf)
  sar_pred <-
    predict(
      probit_sar_train,
      X = fullmat,
      type = "response",
      oos = TRUE,
      WSO = predweightmat
    )[-c(1:nrow(training))]
  
  endtime <- Sys.time()
  sar_time <- as.numeric(difftime(endtime, begintime, units = "mins"))
  
  # run SEM probit and make predictions
  begintime = Sys.time()
  spdf <- training
  coordinates(spdf) <- ~ long + lat
  neib <- knn2nb(knearneigh(coordinates(spdf), longlat = TRUE))
  lw <- nb2listw(neib, style = "W")
  weightmat <- as(as_dgRMatrix_listw(lw), "CsparseMatrix")
  probit_sem_train <-
    ProbitSpatialFit(
      formula = formula(reg_f),
      data = training,
      W = weightmat,
      DGP = "SEM"
    )
  
  fulldf = rbind(training, testing)
  preddf <- fulldf
  coordinates(preddf) <- ~ long + lat
  predneib <- knn2nb(knearneigh(coordinates(preddf), longlat = TRUE))
  predlw <- nb2listw(predneib, style = "W")
  predweightmat <- as(as_dgRMatrix_listw(predlw), "CsparseMatrix")
  fullmat = model.matrix(formula(reg_f), data = fulldf)
  sem_pred <-
    predict(
      probit_sem_train,
      X = fullmat,
      type = "response",
      oos = TRUE,
      WSO = predweightmat
    )[-c(1:nrow(training))]
  
  endtime <- Sys.time()
  sem_time <- as.numeric(difftime(endtime, begintime, units = "mins"))
  
  #run SEM probit and make predictions
  begintime = Sys.time()
  
  classifier = svm(formula = formula(reg_f),
                   data = training,
                   type = 'C-classification',
                   kernel = 'linear')
  svm_pred = predict(classifier, newdata = testing)
  endtime <- Sys.time()
  svm_time <- as.numeric(difftime(endtime,begintime,units = "mins"))
