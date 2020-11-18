  library(nlme)
  library(ggplot2)
  library(stringr)
  library(caret)
  library(RColorBrewer)
  library(tidyverse)
  library(ggfortify)
  library(ggpubr)
  library(factoextra)
  library(dendextend)
  
  ## set wd, set.seed, set digits
  setwd("~/Dokumente/fh/thesis/camda/")
  set.seed(42)
  options(digits=3)
  
  # load and format data and create bacteria-table ----------------------------------------
  # DOH, FAI, nyc, NYC -> 4 cities
  # files contain abundance data of all species found in samples
  
  ilr16_iev <- read.table("dat/abundance_all_S_ILR_16_IEV.tsv", header = T, sep = "\t", quote = "")
  # extract ilr16 
  ilr <- as.data.frame(ilr16_iev[, -48:-96])
  # # (1) delete starting "CSD16" in colnames, (2) rename 1st col to "species"
  names(ilr) <- substring(names(ilr), 6)
  names(ilr)[names(ilr) == ''] <- 'species'
  
  doh_fai_nyc <- read.table("dat/abundance_all_S_DOH_FAI_NYC.tsv", header = T, sep = "\t", quote = "")
  # # (1) delete starting "CSD16" in colnames, (2) rename 1st col to "species"
  names(doh_fai_nyc) <- substring(names(doh_fai_nyc), 6)
  names(doh_fai_nyc)[names(doh_fai_nyc) == ''] <- 'species'
  
  # merge data
  data_all <- merge(ilr, doh_fai_nyc, by = "species", all = T)
  # replace NAs' by zeros
  data_all[is.na(data_all)] <- 0
  
  # include BacDive List and filter for bacteria
  BacDive_Bacteria <- read.csv("dat/BacDive_Bacteria_All.csv")
  BacDive_Bacteria <- as.data.frame(unique(BacDive_Bacteria[, 2]))
  colnames(BacDive_Bacteria) <- "species"
  prokarya_all <- merge(BacDive_Bacteria, data_all, by = "species", all.x = T)
  prokarya_all <- na.omit(prokarya_all)
  
  # format data top abundant bacteria species of each city ----------------------------------------
  # 1st: add a column for the mean of each city to the prokarya-table
  prokarya_extended <- prokarya_all
  prokarya_extended$DOH_median <- apply(select(prokarya_all,starts_with("DOH")), 1, median)
  prokarya_extended$FAI_median <- apply(select(prokarya_all,starts_with("FAI")), 1, median)
  prokarya_extended$ILR_median <- apply(select(prokarya_all,starts_with("ILO")), 1, median)
  prokarya_extended$NYC_median <- apply(select(prokarya_all,starts_with("NYC")), 1, median)
  
  # 2nd step order by city median
  doh_ordered <- prokarya_extended[rev(order(prokarya_extended$DOH_median)),]
  fai_ordered <- prokarya_extended[rev(order(prokarya_extended$FAI_median)),]
  ilr_ordered <- prokarya_extended[rev(order(prokarya_extended$ILR_median)),]
  nyc_ordered <- prokarya_extended[rev(order(prokarya_extended$NYC_median)),]
  
  #sample count and "cities"-vector
  doh_sample_count <- grepl("DOH", names(prokarya_all))
  fai_sample_count <- grepl("FAI", names(prokarya_all))
  ilr_sample_count <- grepl("ILO", names(prokarya_all))
  nyc_sample_count <- grepl("NYC", names(prokarya_all))
  
  sample_count <- c(table(doh_sample_count)["TRUE"], 
                    table(fai_sample_count)["TRUE"], 
                    table(ilr_sample_count)["TRUE"], 
                    table(nyc_sample_count)["TRUE"])
  
  city_list <- c("Doha", "Fairbanks", "Ilorin", "New York")
  sample_count_table <- data.frame(city_list, sample_count)
  sample_count_table
  
  cities <- c(rep("Ilorin", sample_count_table[3, 2]),
              rep("Doha", sample_count_table[1, 2]),
              rep("Fairbanks", sample_count_table[2, 2]),
              rep("New York", sample_count_table[4, 2]))
  
  # choose samples for training and testset
  
  
  # create testset and trainingset -> 1/3 = testset;
  # testset size depending on sample size:48 -> 16, 46 -> 15, 44 -> 15 
  # define how many samples between which colums are taken, according to their positions in the "prokarya_all"-table
  ind_test_DOH <- sample(48:91, size = 15) ###
  ind_test_FAI <- sample(92:139, size = 16) ###
  ind_test_ILR <- sample(2:47, size = 15) ###
  ind_test_NYC <- sample(140:187, size = 16) ###
  ind_test <- c(ind_test_DOH, ind_test_FAI, ind_test_ILR, ind_test_NYC)
  
  # additional variables ----------------------------------------
  rf_colors <- c("darkblue", "darkred", "darkgreen", "darkviolet", "goldenrod4")
  
  # top 10 ----------------------------------------
  # extract top species
  doh_top10 <- as.data.frame(doh_ordered[1:10, 1])
  names(doh_top10)[names(doh_top10) == 'doh_ordered[1:10, 1]'] <- 'species'
  fai_top10 <- as.data.frame(fai_ordered[1:10, 1])
  names(fai_top10)[names(fai_top10) == 'fai_ordered[1:10, 1]'] <- 'species'
  ilr_top10 <- as.data.frame(ilr_ordered[1:10, 1])
  names(ilr_top10)[names(ilr_top10) == 'ilr_ordered[1:10, 1]'] <- 'species'
  nyc_top10 <- as.data.frame(nyc_ordered[1:10, 1])
  names(nyc_top10)[names(nyc_top10) == 'nyc_ordered[1:10, 1]'] <- 'species'
  
  # rbind the top species of all cities and than delete duplicates with "unique"
  top10 <- rbind(doh_top10, fai_top10, ilr_top10, nyc_top10)
  top10 <- as.data.frame(unique(top10))
  top10 <- merge(top10, prokarya_all, by = "species", all.x = T)
  
  # top 10 only species that occur in only one city -> species that appear in 2 cities are dismissed
  # if a species occures in the top list of more than one city it is deleted and not considered for further analysis
  top10_1 <- rbind(doh_top10, fai_top10, ilr_top10, nyc_top10)
  top10_1_species <- table(top10_1$species)
  top10_1 <- subset(top10_1, species %in% names(top10_1_species[top10_1_species < 2]))
  top10_1 <- merge(top10_1, prokarya_all, by = "species", all.x = T)
  
  # only species are kept, that occur at least twice
  top10_2 <- rbind(doh_top10, fai_top10, ilr_top10, nyc_top10)
  top10_2 <- top10_2[duplicated(top10_2),]
  top10_2 <- as.data.frame(unique(top10_2))
  names(top10_2)[names(top10_2) == 'unique(top10_2)'] <- 'species'
  top10_2 <- merge(top10_2, prokarya_all, by = "species", all.x = T)
  
  # top 20 ----------------------------------------
  # extract top species
  doh_top20 <- as.data.frame(doh_ordered[1:20, 1])
  names(doh_top20)[names(doh_top20) == 'doh_ordered[1:20, 1]'] <- 'species'
  fai_top20 <- as.data.frame(fai_ordered[1:20, 1])
  names(fai_top20)[names(fai_top20) == 'fai_ordered[1:20, 1]'] <- 'species'
  ilr_top20 <- as.data.frame(ilr_ordered[1:20, 1])
  names(ilr_top20)[names(ilr_top20) == 'ilr_ordered[1:20, 1]'] <- 'species'
  nyc_top20 <- as.data.frame(nyc_ordered[1:20, 1])
  names(nyc_top20)[names(nyc_top20) == 'nyc_ordered[1:20, 1]'] <- 'species'
  
  # rbind the top species of all cities and than delete duplicates with "unique"
  top20 <- rbind(doh_top20, fai_top20, ilr_top20, nyc_top20)
  top20 <- as.data.frame(unique(top20))
  top20 <- merge(top20, prokarya_all, by = "species", all.x = T)
  
  # top 10 only species that occur in only one city -> species that appear in 2 cities are dismissed
  # if a species occures in the top list of more than one city it is deleted and not considered for further analysis
  top20_1 <- rbind(doh_top20, fai_top20, ilr_top20, nyc_top20)
  top20_1_species <- table(top20_1$species)
  top20_1 <- subset(top20_1, species %in% names(top20_1_species[top20_1_species < 2]))
  top20_1 <- merge(top20_1, prokarya_all, by = "species", all.x = T)
  
  # only species are kept, that occur at least twice
  top20_2 <- rbind(doh_top20, fai_top20, ilr_top20, nyc_top20)
  top20_2 <- top20_2[duplicated(top20_2),]
  top20_2 <- as.data.frame(unique(top20_2))
  names(top20_2)[names(top20_2) == 'unique(top20_2)'] <- 'species'
  top20_2 <- merge(top20_2, prokarya_all, by = "species", all.x = T)
  
  # top 30 ----------------------------------------
  # extract top species
  doh_top30 <- as.data.frame(doh_ordered[1:30, 1])
  names(doh_top30)[names(doh_top30) == 'doh_ordered[1:30, 1]'] <- 'species'
  fai_top30 <- as.data.frame(fai_ordered[1:30, 1])
  names(fai_top30)[names(fai_top30) == 'fai_ordered[1:30, 1]'] <- 'species'
  ilr_top30 <- as.data.frame(ilr_ordered[1:30, 1])
  names(ilr_top30)[names(ilr_top30) == 'ilr_ordered[1:30, 1]'] <- 'species'
  nyc_top30 <- as.data.frame(nyc_ordered[1:30, 1])
  names(nyc_top30)[names(nyc_top30) == 'nyc_ordered[1:30, 1]'] <- 'species'
  
  # rbind the top species of all cities and than delete duplicates with "unique"
  top30 <- rbind(doh_top30, fai_top30, ilr_top30, nyc_top30)
  top30 <- as.data.frame(unique(top30))
  top30 <- merge(top30, prokarya_all, by = "species", all.x = T)
  
  # top 10 only species that occur in only one city -> species that appear in 2 cities are dismissed
  # if a species occures in the top list of more than one city it is deleted and not considered for further analysis
  top30_1 <- rbind(doh_top30, fai_top30, ilr_top30, nyc_top30)
  top30_1_species <- table(top30_1$species)
  top30_1 <- subset(top30_1, species %in% names(top30_1_species[top30_1_species < 2]))
  top30_1 <- merge(top30_1, prokarya_all, by = "species", all.x = T)
  
  # only species are kept, that occur at least twice
  top30_2 <- rbind(doh_top30, fai_top30, ilr_top30, nyc_top30)
  top30_2 <- top30_2[duplicated(top30_2),]
  top30_2 <- as.data.frame(unique(top30_2))
  names(top30_2)[names(top30_2) == 'unique(top30_2)'] <- 'species'
  top30_2 <- merge(top30_2, prokarya_all, by = "species", all.x = T)
  
  # top 100 ----------------------------------------
  # extract top species
  doh_top100 <- as.data.frame(doh_ordered[1:100, 1])
  names(doh_top100)[names(doh_top100) == 'doh_ordered[1:100, 1]'] <- 'species'
  fai_top100 <- as.data.frame(fai_ordered[1:100, 1])
  names(fai_top100)[names(fai_top100) == 'fai_ordered[1:100, 1]'] <- 'species'
  ilr_top100 <- as.data.frame(ilr_ordered[1:100, 1])
  names(ilr_top100)[names(ilr_top100) == 'ilr_ordered[1:100, 1]'] <- 'species'
  nyc_top100 <- as.data.frame(nyc_ordered[1:100, 1])
  names(nyc_top100)[names(nyc_top100) == 'nyc_ordered[1:100, 1]'] <- 'species'
  
  # rbind the top species of all cities and than delete duplicates with "unique"
  top100 <- rbind(doh_top100, fai_top100, ilr_top100, nyc_top100)
  top100 <- as.data.frame(unique(top100))
  top100 <- merge(top100, prokarya_all, by = "species", all.x = T)
  
  # top 10 only species that occur in only one city -> species that appear in 2 cities are dismissed
  # if a species occures in the top list of more than one city it is deleted and not considered for further analysis
  top100_1 <- rbind(doh_top100, fai_top100, ilr_top100, nyc_top100)
  top100_1_species <- table(top100_1$species)
  top100_1 <- subset(top100_1, species %in% names(top100_1_species[top100_1_species < 2]))
  top100_1 <- merge(top100_1, prokarya_all, by = "species", all.x = T)
  
  # only species are kept, that occur at least twice
  top100_2 <- rbind(doh_top100, fai_top100, ilr_top100, nyc_top100)
  top100_2 <- top100_2[duplicated(top100_2),]
  top100_2 <- as.data.frame(unique(top100_2))
  names(top100_2)[names(top100_2) == 'unique(top100_2)'] <- 'species'
  top100_2 <- merge(top100_2, prokarya_all, by = "species", all.x = T)
  # Top 10 analysis ----------------------------------------
  # Top 10 PCA and hierarchical clustering ----------------------------------------
  top10_TS <- as.matrix(top10[,-1])
  rownames(top10_TS) <- top10$species
  pca_top10 <- prcomp(t(top10_TS), scale=TRUE)
  biplot(prcomp(t(top10_TS)), main = "Top 10 species of each city")
  top10_TS2 <- t(top10_TS)
  
  top10_TS3 <- data.frame(top10_TS2)
  top10_TS3$cities <- as.factor(cities)
  
  #library(ggfortify)
  
  autoplot(pca_top10, main = "Top 10 species of each city")
  autoplot(pca_top10, data = top10_TS3, colour = 'cities', main = "Top 10 species of each city")
  autoplot(pca_top10, data = top10_TS3, colour = 'cities', main = "Top 10 species of each city",
           loadings = TRUE, loadings.colour = 'blue',
           loadings.label = TRUE, loadings.label.size = 3)
  
  # loadings of species
  pca_top10$rotation
  
  ### --- hierarchical clustering
  dd_top10 <- dist(scale(top10_TS2), method = "euclidean")
  
  hc_top10 <- hclust(dd_top10, method = "complete")
  plot(hc_top10, main = "Top 10 species of each city")
  
  # top 10 RF ----------------------------------------
  rf_train_top10 <- top10_TS3[-ind_test, ]
  rf_test_top10 <- top10_TS3[ind_test, ]
  rf_test_top10 <- na.omit(rf_test_top10)
  
  # replace spaces and "-" by "_"
  #names(train_set_top10) <- str_replace_all(names(train_set_top10), c(" " = "_", "-" = "_"))
  #names(test_set_top10) <- str_replace_all(names(test_set_top10), c(" " = "_", "-" = "_"))
  
  # use caret with random forest as my model with 5 fold cross validation
  rf_model_top10<-train(cities~.,data=rf_train_top10,method="rf",
                        trControl=trainControl(method="cv",number=5),
                        prox=TRUE,allowParallel=TRUE)
  
  # print rf model
  print(rf_model_top10)
  print(rf_model_top10$finalModel)
  
  # plot rf model
  plot(rf_model_top10)
  plot(rf_model_top10$finalModel, main = "Top 10", col =rf_colors, lwd = 1)
  #legend("topright", ,fill =rf_colors,cex=0.8)
  
  # variable importance plot
  plot(varImp(rf_model_top10, scale = FALSE))
  plot(varImp(rf_model_top10, scale = TRUE))
  
  rf_pred_top10 <- predict(rf_model_top10, newdata = rf_test_top10)
  rf_pred_top10
  
  confusionMatrix(rf_pred_top10, rf_test_top10$cities)
  
  # top10 kNN ----------------------------------------
  knn_train_top10 <- top10_TS3[-ind_test, ]
  knn_test_top10 <- top10_TS3[ind_test, ]
  knn_test_top10 <- na.omit(knn_test_top10)
  
  # start knn
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  knn_fit_model_top10 <- train(cities ~., data = knn_train_top10, method = "knn",
                               trControl=trctrl,
                               preProcess = c("center", "scale"),
                               tuneLength = 10)
  
  # watch output
  knn_fit_model_top10
  plot(knn_fit_model_top10)
  
  # test set prediction
  knn_pred_top10 <- predict(knn_fit_model_top10, newdata = knn_test_top10)
  knn_pred_top10
  
  # confusion matrix
  confusionMatrix(knn_pred_top10, knn_test_top10$cities)
  
  # top10 support vector machines ----------------------------------------
  
  svm_fit_model_top10 <- train(cities ~., data = knn_train_top10, method = "svmRadial", trControl = trctrl, preProcess = c("center","scale"), tuneLength = 10)
  
  # watch output
  svm_fit_model_top10
  plot(svm_fit_model_top10)
  
  # test set prediction
  svm_pred_top10 <- predict(svm_fit_model_top10, newdata = knn_test_top10)
  svm_pred_top10
  
  # confusion matrix
  confusionMatrix(svm_pred_top10, knn_test_top10$cities)
  
  # Top 10_1 analysis ----------------------------------------
  # Top 10_1 PCA and hierarchical clustering ----------------------------------------
  top10_1_TS <- as.matrix(top10_1[,-1])
  rownames(top10_1_TS) <- top10_1$species
  pca_top10_1 <- prcomp(t(top10_1_TS), scale=TRUE)
  biplot(prcomp(t(top10_1_TS)), main = "Top 10 species of each city")
  top10_1_TS2 <- t(top10_1_TS)
  
  top10_1_TS3 <- data.frame(top10_1_TS2)
  top10_1_TS3$cities <- as.factor(cities)
  
  #library(ggfortify)
  
  autoplot(pca_top10_1, main = "Top 10 species of each city")
  autoplot(pca_top10_1, data = top10_1_TS3, colour = 'cities', main = "Top 10 species of each city")
  autoplot(pca_top10_1, data = top10_1_TS3, colour = 'cities', main = "Top 10 species of each city",
           loadings = TRUE, loadings.colour = 'blue',
           loadings.label = TRUE, loadings.label.size = 3)
  
  # loadings of species
  pca_top10_1$rotation
  
  ### --- hierarchical clustering
  dd_top10_1 <- dist(scale(top10_1_TS2), method = "euclidean")
  
  hc_top10_1 <- hclust(dd_top10_1, method = "complete")
  plot(hc_top10_1, main = "Top 10 species of each city")
  
  # top 10_1 RF ----------------------------------------
  rf_train_top10_1 <- top10_1_TS3[-ind_test, ]
  rf_test_top10_1 <- top10_1_TS3[ind_test, ] # why do I get one row with only NAs?
  rf_test_top10_1 <- na.omit(rf_test_top10_1)
  
  # replace spaces and "-" by "_"
  #names(train_set_top10_1) <- str_replace_all(names(train_set_top10_1), c(" " = "_", "-" = "_"))
  #names(test_set_top10_1) <- str_replace_all(names(test_set_top10_1), c(" " = "_", "-" = "_"))
  
  # use caret with random forest as my model with 5 fold cross validation
  rf_model_top10_1<-train(cities~.,data=rf_train_top10_1,method="rf",
                          trControl=trainControl(method="cv",number=5),
                          prox=TRUE,allowParallel=TRUE)
  
  print(rf_model_top10_1)
  print(rf_model_top10_1$finalModel)
  plot(rf_model_top10_1)
  plot(rf_model_top10_1$finalModel)
  
  rf_pred_top10_1 <- predict(rf_model_top10_1, newdata = rf_test_top10_1)
  rf_pred_top10_1
  
  confusionMatrix(rf_pred_top10_1, rf_test_top10_1$cities)
  # top10_1 kNN ----------------------------------------
  knn_train_top10_1 <- top10_1_TS3[-ind_test, ]
  knn_test_top10_1 <- top10_1_TS3[ind_test, ]
  knn_test_top10_1 <- na.omit(knn_test_top10_1)
  
  # start knn
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  knn_fit_model_top10_1 <- train(cities ~., data = knn_train_top10_1, method = "knn",
                                 trControl=trctrl,
                                 preProcess = c("center", "scale"),
                                 tuneLength = 10)
  
  # watch output
  knn_fit_model_top10_1
  plot(knn_fit_model_top10_1)
  
  # test set prediction
  knn_pred_top10_1 <- predict(knn_fit_model_top10_1, newdata = knn_test_top10_1)
  knn_pred_top10_1
  
  # confusion matrix
  confusionMatrix(knn_pred_top10_1, knn_test_top10_1$cities)
  
  # top10_1 support vector machines ----------------------------------------
  
  svm_fit_model_top10_1 <- train(cities ~., data = knn_train_top10_1, method = "svmRadial", trControl = trctrl, preProcess = c("center","scale"), tuneLength = 10)
  
  # watch output
  svm_fit_model_top10_1
  plot(svm_fit_model_top10_1)
  
  # test set prediction
  svm_pred_top10_1 <- predict(svm_fit_model_top10_1, newdata = knn_test_top10_1)
  svm_pred_top10_1
  
  # confusion matrix
  confusionMatrix(svm_pred_top10_1, knn_test_top10_1$cities)
  
  # Top 10_2 analysis ----------------------------------------
  # Top 10_2 PCA and hierarchical clustering ----------------------------------------
  top10_2_TS <- as.matrix(top10_2[,-1])
  rownames(top10_2_TS) <- top10_2$species
  pca_top10_2 <- prcomp(t(top10_2_TS), scale=TRUE)
  biplot(prcomp(t(top10_2_TS)), main = "Top 10 species of each city")
  top10_2_TS2 <- t(top10_2_TS)
  
  top10_2_TS3 <- data.frame(top10_2_TS2)
  top10_2_TS3$cities <- as.factor(cities)
  
  #library(ggfortify)
  
  autoplot(pca_top10_2, main = "Top 10 species of each city")
  autoplot(pca_top10_2, data = top10_2_TS3, colour = 'cities', main = "Top 10 species of each city")
  autoplot(pca_top10_2, data = top10_2_TS3, colour = 'cities', main = "Top 10 species of each city",
           loadings = TRUE, loadings.colour = 'blue',
           loadings.label = TRUE, loadings.label.size = 3)
  
  # loadings of species
  pca_top10_2$rotation
  
  ### --- hierarchical clustering
  dd_top10_2 <- dist(scale(top10_2_TS2), method = "euclidean")
  
  hc_top10_2 <- hclust(dd_top10_2, method = "complete")
  plot(hc_top10_2, main = "Top 10 species of each city")
  
  # top 10_2 RF ----------------------------------------
  rf_train_top10_2 <- top10_2_TS3[-ind_test, ]
  rf_test_top10_2 <- top10_2_TS3[ind_test, ] # why do I get one row with only NAs?
  rf_test_top10_2 <- na.omit(rf_test_top10_2)
  
  # replace spaces and "-" by "_"
  #names(train_set_top10_2) <- str_replace_all(names(train_set_top10_2), c(" " = "_", "-" = "_"))
  #names(test_set_top10_2) <- str_replace_all(names(test_set_top10_2), c(" " = "_", "-" = "_"))
  
  # use caret with random forest as my model with 5 fold cross validation
  rf_model_top10_2<-train(cities~.,data=rf_train_top10_2,method="rf",
                          trControl=trainControl(method="cv",number=5),
                          prox=TRUE,allowParallel=TRUE)
  
  print(rf_model_top10_2)
  print(rf_model_top10_2$finalModel)
  plot(rf_model_top10_2)
  plot(rf_model_top10_2$finalModel)
  
  rf_pred_top10_2 <- predict(rf_model_top10_2, newdata = rf_test_top10_2)
  rf_pred_top10_2
  
  confusionMatrix(rf_pred_top10_2, rf_test_top10_2$cities)
  # top10_2 kNN ----------------------------------------
  knn_train_top10_2 <- top10_2_TS3[-ind_test, ]
  knn_test_top10_2 <- top10_2_TS3[ind_test, ]
  knn_test_top10_2 <- na.omit(knn_test_top10_2)
  
  # start knn
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  knn_fit_model_top10_2 <- train(cities ~., data = knn_train_top10_2, method = "knn",
                                 trControl=trctrl,
                                 preProcess = c("center", "scale"),
                                 tuneLength = 10)
  
  # watch output
  knn_fit_model_top10_2
  plot(knn_fit_model_top10_2)
  
  # test set prediction
  knn_pred_top10_2 <- predict(knn_fit_model_top10_2, newdata = knn_test_top10_2)
  knn_pred_top10_2
  
  # confusion matrix
  confusionMatrix(knn_pred_top10_2, knn_test_top10_2$cities)
  
  # top10_2 support vector machines ----------------------------------------
  
  svm_fit_model_top10_2 <- train(cities ~., data = knn_train_top10_2, method = "svmRadial", trControl = trctrl, preProcess = c("center","scale"), tuneLength = 10)
  
  # watch output
  svm_fit_model_top10_2
  plot(svm_fit_model_top10_2)
  
  # test set prediction
  svm_pred_top10_2 <- predict(svm_fit_model_top10_2, newdata = knn_test_top10_2)
  svm_pred_top10_2
  
  # confusion matrix
  confusionMatrix(svm_pred_top10_2, knn_test_top10_2$cities)
  
  # Top 30 analysis --------------------------------
  # Top 20 PCA and hierarchical clustering ----------------------------------------
  top20_TS <- as.matrix(top20[,-1])
  rownames(top20_TS) <- top20$species
  pca_top20 <- prcomp(t(top20_TS), scale=TRUE)
  biplot(prcomp(t(top20_TS)), main = "Top 10 species of each city")
  top20_TS2 <- t(top20_TS)
  
  top20_TS3 <- data.frame(top20_TS2)
  top20_TS3$cities <- as.factor(cities)
  
  #library(ggfortify)
  
  autoplot(pca_top20, main = "Top 10 species of each city")
  autoplot(pca_top20, data = top20_TS3, colour = 'cities', main = "Top 10 species of each city")
  autoplot(pca_top20, data = top20_TS3, colour = 'cities', main = "Top 10 species of each city",
           loadings = TRUE, loadings.colour = 'blue',
           loadings.label = TRUE, loadings.label.size = 3)
  
  # loadings of species
  pca_top20$rotation
  
  ### --- hierarchical clustering
  dd_top20 <- dist(scale(top20_TS2), method = "euclidean")
  
  hc_top20 <- hclust(dd_top20, method = "complete")
  plot(hc_top20, main = "Top 10 species of each city")
  
  # top 20 RF ----------------------------------------
  rf_train_top20 <- top20_TS3[-ind_test, ]
  rf_test_top20 <- top20_TS3[ind_test, ] # why do I get one row with only NAs?
  rf_test_top20 <- na.omit(rf_test_top20)
  
  # replace spaces and "-" by "_"
  #names(train_set_top20) <- str_replace_all(names(train_set_top20), c(" " = "_", "-" = "_"))
  #names(test_set_top20) <- str_replace_all(names(test_set_top20), c(" " = "_", "-" = "_"))
  
  # use caret with random forest as my model with 5 fold cross validation
  rf_model_top20<-train(cities~.,data=rf_train_top20,method="rf",
                        trControl=trainControl(method="cv",number=5),
                        prox=TRUE,allowParallel=TRUE)
  
  print(rf_model_top20)
  print(rf_model_top20$finalModel)
  plot(rf_model_top20)
  plot(rf_model_top20$finalModel)
  
  rf_pred_top20 <- predict(rf_model_top20, newdata = rf_test_top20)
  rf_pred_top20
  
  confusionMatrix(rf_pred_top20, rf_test_top20$cities)
  
  # top20 kNN ----------------------------------------
  knn_train_top20 <- top20_TS3[-ind_test, ]
  knn_test_top20 <- top20_TS3[ind_test, ]
  knn_test_top20 <- na.omit(knn_test_top20)
  
  # start knn
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  knn_fit_model_top20 <- train(cities ~., data = knn_train_top20, method = "knn",
                               trControl=trctrl,
                               preProcess = c("center", "scale"),
                               tuneLength = 10)
  
  # watch output
  knn_fit_model_top20
  plot(knn_fit_model_top20)
  
  # test set prediction
  knn_pred_top20 <- predict(knn_fit_model_top20, newdata = knn_test_top20)
  knn_pred_top20
  
  # confusion matrix
  confusionMatrix(knn_pred_top20, knn_test_top20$cities)
  
  # top20 support vector machines ----------------------------------------
  
  svm_fit_model_top20 <- train(cities ~., data = knn_train_top20, method = "svmRadial", trControl = trctrl, preProcess = c("center","scale"), tuneLength = 10)
  
  # watch output
  svm_fit_model_top20
  plot(svm_fit_model_top20)
  
  # test set prediction
  svm_pred_top20 <- predict(svm_fit_model_top20, newdata = knn_test_top20)
  svm_pred_top20
  
  # confusion matrix
  confusionMatrix(svm_pred_top20, knn_test_top20$cities)
  
  # Top 20_1 analysis ----------------------------------------
  # Top 20_1 PCA and hierarchical clustering ----------------------------------------
  top20_1_TS <- as.matrix(top20_1[,-1])
  rownames(top20_1_TS) <- top20_1$species
  pca_top20_1 <- prcomp(t(top20_1_TS), scale=TRUE)
  biplot(prcomp(t(top20_1_TS)), main = "Top 10 species of each city")
  top20_1_TS2 <- t(top20_1_TS)
  
  top20_1_TS3 <- data.frame(top20_1_TS2)
  top20_1_TS3$cities <- as.factor(cities)
  
  #library(ggfortify)
  
  autoplot(pca_top20_1, main = "Top 10 species of each city")
  autoplot(pca_top20_1, data = top20_1_TS3, colour = 'cities', main = "Top 10 species of each city")
  autoplot(pca_top20_1, data = top20_1_TS3, colour = 'cities', main = "Top 10 species of each city",
           loadings = TRUE, loadings.colour = 'blue',
           loadings.label = TRUE, loadings.label.size = 3)
  
  # loadings of species
  pca_top20_1$rotation
  
  ### --- hierarchical clustering
  dd_top20_1 <- dist(scale(top20_1_TS2), method = "euclidean")
  
  hc_top20_1 <- hclust(dd_top20_1, method = "complete")
  plot(hc_top20_1, main = "Top 10 species of each city")
  
  # top 20_1 RF ----------------------------------------
  rf_train_top20_1 <- top20_1_TS3[-ind_test, ]
  rf_test_top20_1 <- top20_1_TS3[ind_test, ] # why do I get one row with only NAs?
  rf_test_top20_1 <- na.omit(rf_test_top20_1)
  
  # replace spaces and "-" by "_"
  #names(train_set_top20_1) <- str_replace_all(names(train_set_top20_1), c(" " = "_", "-" = "_"))
  #names(test_set_top20_1) <- str_replace_all(names(test_set_top20_1), c(" " = "_", "-" = "_"))
  
  # use caret with random forest as my model with 5 fold cross validation
  rf_model_top20_1<-train(cities~.,data=rf_train_top20_1,method="rf",
                          trControl=trainControl(method="cv",number=5),
                          prox=TRUE,allowParallel=TRUE)
  
  print(rf_model_top20_1)
  print(rf_model_top20_1$finalModel)
  plot(rf_model_top20_1)
  plot(rf_model_top20_1$finalModel)
  
  rf_pred_top20_1 <- predict(rf_model_top20_1, newdata = rf_test_top20_1)
  rf_pred_top20_1
  
  confusionMatrix(rf_pred_top20_1, rf_test_top20_1$cities)
  # top20_1 kNN ----------------------------------------
  knn_train_top20_1 <- top20_1_TS3[-ind_test, ]
  knn_test_top20_1 <- top20_1_TS3[ind_test, ]
  knn_test_top20_1 <- na.omit(knn_test_top20_1)
  
  # start knn
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  knn_fit_model_top20_1 <- train(cities ~., data = knn_train_top20_1, method = "knn",
                                 trControl=trctrl,
                                 preProcess = c("center", "scale"),
                                 tuneLength = 10)
  
  # watch output
  knn_fit_model_top20_1
  plot(knn_fit_model_top20_1)
  
  # test set prediction
  knn_pred_top20_1 <- predict(knn_fit_model_top20_1, newdata = knn_test_top20_1)
  knn_pred_top20_1
  
  # confusion matrix
  confusionMatrix(knn_pred_top20_1, knn_test_top20_1$cities)
  
  # top20_1 support vector machines ----------------------------------------
  
  svm_fit_model_top20_1 <- train(cities ~., data = knn_train_top20_1, method = "svmRadial", trControl = trctrl, preProcess = c("center","scale"), tuneLength = 10)
  
  # watch output
  svm_fit_model_top20_1
  plot(svm_fit_model_top20_1)
  
  # test set prediction
  svm_pred_top20_1 <- predict(svm_fit_model_top20_1, newdata = knn_test_top20_1)
  svm_pred_top20_1
  
  # confusion matrix
  confusionMatrix(svm_pred_top20_1, knn_test_top20_1$cities)
  
  # Top 20_2 analysis ----------------------------------------
  # Top 20_2 PCA and hierarchical clustering ----------------------------------------
  top20_2_TS <- as.matrix(top20_2[,-1])
  rownames(top20_2_TS) <- top20_2$species
  pca_top20_2 <- prcomp(t(top20_2_TS), scale=TRUE)
  biplot(prcomp(t(top20_2_TS)), main = "Top 10 species of each city")
  top20_2_TS2 <- t(top20_2_TS)
  
  top20_2_TS3 <- data.frame(top20_2_TS2)
  top20_2_TS3$cities <- as.factor(cities)
  
  #library(ggfortify)
  
  autoplot(pca_top20_2, main = "Top 10 species of each city")
  autoplot(pca_top20_2, data = top20_2_TS3, colour = 'cities', main = "Top 10 species of each city")
  autoplot(pca_top20_2, data = top20_2_TS3, colour = 'cities', main = "Top 10 species of each city",
           loadings = TRUE, loadings.colour = 'blue',
           loadings.label = TRUE, loadings.label.size = 3)
  
  # loadings of species
  pca_top20_2$rotation
  
  ### --- hierarchical clustering
  dd_top20_2 <- dist(scale(top20_2_TS2), method = "euclidean")
  
  hc_top20_2 <- hclust(dd_top20_2, method = "complete")
  plot(hc_top20_2, main = "Top 10 species of each city")
  
  # top 20_2 RF ----------------------------------------
  rf_train_top20_2 <- top20_2_TS3[-ind_test, ]
  rf_test_top20_2 <- top20_2_TS3[ind_test, ] # why do I get one row with only NAs?
  rf_test_top20_2 <- na.omit(rf_test_top20_2)
  
  # replace spaces and "-" by "_"
  #names(train_set_top20_2) <- str_replace_all(names(train_set_top20_2), c(" " = "_", "-" = "_"))
  #names(test_set_top20_2) <- str_replace_all(names(test_set_top20_2), c(" " = "_", "-" = "_"))
  
  # use caret with random forest as my model with 5 fold cross validation
  rf_model_top20_2<-train(cities~.,data=rf_train_top20_2,method="rf",
                          trControl=trainControl(method="cv",number=5),
                          prox=TRUE,allowParallel=TRUE)
  
  print(rf_model_top20_2)
  print(rf_model_top20_2$finalModel)
  plot(rf_model_top20_2)
  plot(rf_model_top20_2$finalModel)
  
  rf_pred_top20_2 <- predict(rf_model_top20_2, newdata = rf_test_top20_2)
  rf_pred_top20_2
  
  confusionMatrix(rf_pred_top20_2, rf_test_top20_2$cities)
  # top20_2 kNN ----------------------------------------
  knn_train_top20_2 <- top20_2_TS3[-ind_test, ]
  knn_test_top20_2 <- top20_2_TS3[ind_test, ]
  knn_test_top20_2 <- na.omit(knn_test_top20_2)
  
  # start knn
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  knn_fit_model_top20_2 <- train(cities ~., data = knn_train_top20_2, method = "knn",
                                 trControl=trctrl,
                                 preProcess = c("center", "scale"),
                                 tuneLength = 10)
  
  # watch output
  knn_fit_model_top20_2
  plot(knn_fit_model_top20_2)
  
  # test set prediction
  knn_pred_top20_2 <- predict(knn_fit_model_top20_2, newdata = knn_test_top20_2)
  knn_pred_top20_2
  
  # confusion matrix
  confusionMatrix(knn_pred_top20_2, knn_test_top20_2$cities)
  
  # top20_2 support vector machines ----------------------------------------
  
  svm_fit_model_top20_2 <- train(cities ~., data = knn_train_top20_2, method = "svmRadial", trControl = trctrl, preProcess = c("center","scale"), tuneLength = 10)
  
  # watch output
  svm_fit_model_top20_2
  plot(svm_fit_model_top20_2)
  
  # test set prediction
  svm_pred_top20_2 <- predict(svm_fit_model_top20_2, newdata = knn_test_top20_2)
  svm_pred_top20_2
  
  # confusion matrix
  confusionMatrix(svm_pred_top20_2, knn_test_top20_2$cities)
  # top 30 analysis ----------------------------------------
  # Top 30 PCA and hierarchical clustering ----------------------------------------
  top30_TS <- as.matrix(top30[,-1])
  rownames(top30_TS) <- top30$species
  pca_top30 <- prcomp(t(top30_TS), scale=TRUE)
  biplot(prcomp(t(top30_TS)), main = "Top 10 species of each city")
  top30_TS2 <- t(top30_TS)
  
  top30_TS3 <- data.frame(top30_TS2)
  top30_TS3$cities <- as.factor(cities)
  
  #library(ggfortify)
  
  autoplot(pca_top30, main = "Top 10 species of each city")
  autoplot(pca_top30, data = top30_TS3, colour = 'cities', main = "Top 10 species of each city")
  autoplot(pca_top30, data = top30_TS3, colour = 'cities', main = "Top 10 species of each city",
           loadings = TRUE, loadings.colour = 'blue',
           loadings.label = TRUE, loadings.label.size = 3)
  
  # loadings of species
  pca_top30$rotation
  
  ### --- hierarchical clustering
  dd_top30 <- dist(scale(top30_TS2), method = "euclidean")
  
  hc_top30 <- hclust(dd_top30, method = "complete")
  plot(hc_top30, main = "Top 10 species of each city")
  
  # top 30 RF ----------------------------------------
  rf_train_top30 <- top30_TS3[-ind_test, ]
  rf_test_top30 <- top30_TS3[ind_test, ] # why do I get one row with only NAs?
  rf_test_top30 <- na.omit(rf_test_top30)
  
  # replace spaces and "-" by "_"
  #names(train_set_top30) <- str_replace_all(names(train_set_top30), c(" " = "_", "-" = "_"))
  #names(test_set_top30) <- str_replace_all(names(test_set_top30), c(" " = "_", "-" = "_"))
  
  # use caret with random forest as my model with 5 fold cross validation
  rf_model_top30<-train(cities~.,data=rf_train_top30,method="rf",
                        trControl=trainControl(method="cv",number=5),
                        prox=TRUE,allowParallel=TRUE)
  
  print(rf_model_top30)
  print(rf_model_top30$finalModel)
  plot(rf_model_top30)
  plot(rf_model_top30$finalModel)
  
  rf_pred_top30 <- predict(rf_model_top30, newdata = rf_test_top30)
  rf_pred_top30
  
  confusionMatrix(rf_pred_top30, rf_test_top30$cities)
  
  # top30 kNN ----------------------------------------
  knn_train_top30 <- top30_TS3[-ind_test, ]
  knn_test_top30 <- top30_TS3[ind_test, ]
  knn_test_top30 <- na.omit(knn_test_top30)
  
  # start knn
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  knn_fit_model_top30 <- train(cities ~., data = knn_train_top30, method = "knn",
                               trControl=trctrl,
                               preProcess = c("center", "scale"),
                               tuneLength = 10)
  
  # watch output
  knn_fit_model_top30
  plot(knn_fit_model_top30)
  
  # test set prediction
  knn_pred_top30 <- predict(knn_fit_model_top30, newdata = knn_test_top30)
  knn_pred_top30
  
  # confusion matrix
  confusionMatrix(knn_pred_top30, knn_test_top30$cities)
  
  # top30 support vector machines ----------------------------------------
  
  svm_fit_model_top30 <- train(cities ~., data = knn_train_top30, method = "svmRadial", trControl = trctrl, preProcess = c("center","scale"), tuneLength = 10)
  
  # watch output
  svm_fit_model_top30
  plot(svm_fit_model_top30)
  
  # test set prediction
  svm_pred_top30 <- predict(svm_fit_model_top30, newdata = knn_test_top30)
  svm_pred_top30
  
  # confusion matrix
  confusionMatrix(svm_pred_top30, knn_test_top30$cities)
  
  # Top 30_1 analysis ----------------------------------------
  # Top 30_1 PCA and hierarchical clustering ----------------------------------------
  top30_1_TS <- as.matrix(top30_1[,-1])
  rownames(top30_1_TS) <- top30_1$species
  pca_top30_1 <- prcomp(t(top30_1_TS), scale=TRUE)
  biplot(prcomp(t(top30_1_TS)), main = "Top 10 species of each city")
  top30_1_TS2 <- t(top30_1_TS)
  
  top30_1_TS3 <- data.frame(top30_1_TS2)
  top30_1_TS3$cities <- as.factor(cities)
  
  #library(ggfortify)
  
  autoplot(pca_top30_1, main = "Top 10 species of each city")
  autoplot(pca_top30_1, data = top30_1_TS3, colour = 'cities', main = "Top 10 species of each city")
  autoplot(pca_top30_1, data = top30_1_TS3, colour = 'cities', main = "Top 10 species of each city",
           loadings = TRUE, loadings.colour = 'blue',
           loadings.label = TRUE, loadings.label.size = 3)
  
  # loadings of species
  pca_top30_1$rotation
  
  ### --- hierarchical clustering
  dd_top30_1 <- dist(scale(top30_1_TS2), method = "euclidean")
  
  hc_top30_1 <- hclust(dd_top30_1, method = "complete")
  plot(hc_top30_1, main = "Top 10 species of each city")
  
  # top 30_1 RF ----------------------------------------
  rf_train_top30_1 <- top30_1_TS3[-ind_test, ]
  rf_test_top30_1 <- top30_1_TS3[ind_test, ] # why do I get one row with only NAs?
  rf_test_top30_1 <- na.omit(rf_test_top30_1)
  
  # replace spaces and "-" by "_"
  #names(train_set_top30_1) <- str_replace_all(names(train_set_top30_1), c(" " = "_", "-" = "_"))
  #names(test_set_top30_1) <- str_replace_all(names(test_set_top30_1), c(" " = "_", "-" = "_"))
  
  # use caret with random forest as my model with 5 fold cross validation
  rf_model_top30_1<-train(cities~.,data=rf_train_top30_1,method="rf",
                          trControl=trainControl(method="cv",number=5),
                          prox=TRUE,allowParallel=TRUE)
  
  print(rf_model_top30_1)
  print(rf_model_top30_1$finalModel)
  plot(rf_model_top30_1)
  plot(rf_model_top30_1$finalModel)
  
  rf_pred_top30_1 <- predict(rf_model_top30_1, newdata = rf_test_top30_1)
  rf_pred_top30_1
  
  confusionMatrix(rf_pred_top30_1, rf_test_top30_1$cities)
  # top30_1 kNN ----------------------------------------
  knn_train_top30_1 <- top30_1_TS3[-ind_test, ]
  knn_test_top30_1 <- top30_1_TS3[ind_test, ]
  knn_test_top30_1 <- na.omit(knn_test_top30_1)
  
  # start knn
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  knn_fit_model_top30_1 <- train(cities ~., data = knn_train_top30_1, method = "knn",
                                 trControl=trctrl,
                                 preProcess = c("center", "scale"),
                                 tuneLength = 10)
  
  # watch output
  knn_fit_model_top30_1
  plot(knn_fit_model_top30_1)
  
  # test set prediction
  knn_pred_top30_1 <- predict(knn_fit_model_top30_1, newdata = knn_test_top30_1)
  knn_pred_top30_1
  
  # confusion matrix
  confusionMatrix(knn_pred_top30_1, knn_test_top30_1$cities)
  
  # top30_1 support vector machines ----------------------------------------
  
  svm_fit_model_top30_1 <- train(cities ~., data = knn_train_top30_1, method = "svmRadial", trControl = trctrl, preProcess = c("center","scale"), tuneLength = 10)
  
  # watch output
  svm_fit_model_top30_1
  plot(svm_fit_model_top30_1)
  
  # test set prediction
  svm_pred_top30_1 <- predict(svm_fit_model_top30_1, newdata = knn_test_top30_1)
  svm_pred_top30_1
  
  # confusion matrix
  confusionMatrix(svm_pred_top30_1, knn_test_top30_1$cities)
  
  # Top 30_2 analysis ----------------------------------------
  # Top 30_2 PCA and hierarchical clustering ----------------------------------------
  top30_2_TS <- as.matrix(top30_2[,-1])
  rownames(top30_2_TS) <- top30_2$species
  pca_top30_2 <- prcomp(t(top30_2_TS), scale=TRUE)
  biplot(prcomp(t(top30_2_TS)), main = "Top 10 species of each city")
  top30_2_TS2 <- t(top30_2_TS)
  
  top30_2_TS3 <- data.frame(top30_2_TS2)
  top30_2_TS3$cities <- as.factor(cities)
  
  #library(ggfortify)
  
  autoplot(pca_top30_2, main = "Top 10 species of each city")
  autoplot(pca_top30_2, data = top30_2_TS3, colour = 'cities', main = "Top 10 species of each city")
  autoplot(pca_top30_2, data = top30_2_TS3, colour = 'cities', main = "Top 10 species of each city",
           loadings = TRUE, loadings.colour = 'blue',
           loadings.label = TRUE, loadings.label.size = 3)
  
  # loadings of species
  pca_top30_2$rotation
  
  ### --- hierarchical clustering
  dd_top30_2 <- dist(scale(top30_2_TS2), method = "euclidean")
  
  hc_top30_2 <- hclust(dd_top30_2, method = "complete")
  plot(hc_top30_2, main = "Top 10 species of each city")
  
  # top 30_2 RF ----------------------------------------
  rf_train_top30_2 <- top30_2_TS3[-ind_test, ]
  rf_test_top30_2 <- top30_2_TS3[ind_test, ] # why do I get one row with only NAs?
  rf_test_top30_2 <- na.omit(rf_test_top30_2)
  
  # replace spaces and "-" by "_"
  #names(train_set_top30_2) <- str_replace_all(names(train_set_top30_2), c(" " = "_", "-" = "_"))
  #names(test_set_top30_2) <- str_replace_all(names(test_set_top30_2), c(" " = "_", "-" = "_"))
  
  # use caret with random forest as my model with 5 fold cross validation
  rf_model_top30_2<-train(cities~.,data=rf_train_top30_2,method="rf",
                          trControl=trainControl(method="cv",number=5),
                          prox=TRUE,allowParallel=TRUE)
  
  print(rf_model_top30_2)
  print(rf_model_top30_2$finalModel)
  plot(rf_model_top30_2)
  plot(rf_model_top30_2$finalModel)
  
  rf_pred_top30_2 <- predict(rf_model_top30_2, newdata = rf_test_top30_2)
  rf_pred_top30_2
  
  confusionMatrix(rf_pred_top30_2, rf_test_top30_2$cities)
  # top30_2 kNN ----------------------------------------
  knn_train_top30_2 <- top30_2_TS3[-ind_test, ]
  knn_test_top30_2 <- top30_2_TS3[ind_test, ]
  knn_test_top30_2 <- na.omit(knn_test_top30_2)
  
  # start knn
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  knn_fit_model_top30_2 <- train(cities ~., data = knn_train_top30_2, method = "knn",
                                 trControl=trctrl,
                                 preProcess = c("center", "scale"),
                                 tuneLength = 10)
  
  # watch output
  knn_fit_model_top30_2
  plot(knn_fit_model_top30_2)
  
  # test set prediction
  knn_pred_top30_2 <- predict(knn_fit_model_top30_2, newdata = knn_test_top30_2)
  knn_pred_top30_2
  
  # confusion matrix
  confusionMatrix(knn_pred_top30_2, knn_test_top30_2$cities)
  
  # top30_2 support vector machines ----------------------------------------
  
  svm_fit_model_top30_2 <- train(cities ~., data = knn_train_top30_2, method = "svmRadial", trControl = trctrl, preProcess = c("center","scale"), tuneLength = 10)
  
  # watch output
  svm_fit_model_top30_2
  plot(svm_fit_model_top30_2)
  
  # test set prediction
  svm_pred_top30_2 <- predict(svm_fit_model_top30_2, newdata = knn_test_top30_2)
  svm_pred_top30_2
  
  # confusion matrix
  confusionMatrix(svm_pred_top30_2, knn_test_top30_2$cities)
  
  # Top 100 analysis ----------------------------------------
  # Top 100 PCA and hierarchical clustering ----------------------------------------
  top100_TS <- as.matrix(top100[,-1])
  rownames(top100_TS) <- top100$species
  pca_top100 <- prcomp(t(top100_TS), scale=TRUE)
  biplot(prcomp(t(top100_TS)), main = "Top 10 species of each city")
  top100_TS2 <- t(top100_TS)
  
  top100_TS3 <- data.frame(top100_TS2)
  top100_TS3$cities <- as.factor(cities)
  
  #library(ggfortify)
  
  autoplot(pca_top100, main = "Top 10 species of each city")
  autoplot(pca_top100, data = top100_TS3, colour = 'cities', main = "Top 10 species of each city")
  autoplot(pca_top100, data = top100_TS3, colour = 'cities', main = "Top 10 species of each city",
           loadings = TRUE, loadings.colour = 'blue',
           loadings.label = TRUE, loadings.label.size = 3)
  
  # loadings of species
  pca_top100$rotation
  
  ### --- hierarchical clustering
  dd_top100 <- dist(scale(top100_TS2), method = "euclidean")
  
  hc_top100 <- hclust(dd_top100, method = "complete")
  plot(hc_top100, main = "Top 10 species of each city")
  
  # top 100 RF ----------------------------------------
  rf_train_top100 <- top100_TS3[-ind_test, ]
  rf_test_top100 <- top100_TS3[ind_test, ] # why do I get one row with only NAs?
  rf_test_top100 <- na.omit(rf_test_top100)
  
  # replace spaces and "-" by "_"
  #names(train_set_top100) <- str_replace_all(names(train_set_top100), c(" " = "_", "-" = "_"))
  #names(test_set_top100) <- str_replace_all(names(test_set_top100), c(" " = "_", "-" = "_"))
  
  # use caret with random forest as my model with 5 fold cross validation
  rf_model_top100<-train(cities~.,data=rf_train_top100,method="rf",
                         trControl=trainControl(method="cv",number=5),
                         prox=TRUE,allowParallel=TRUE)
  
  print(rf_model_top100)
  print(rf_model_top100$finalModel)
  plot(rf_model_top100)
  plot(rf_model_top100$finalModel)
  
  rf_pred_top100 <- predict(rf_model_top100, newdata = rf_test_top100)
  rf_pred_top100
  
  confusionMatrix(rf_pred_top100, rf_test_top100$cities)
  
  # top100 kNN ----------------------------------------
  knn_train_top100 <- top100_TS3[-ind_test, ]
  knn_test_top100 <- top100_TS3[ind_test, ]
  knn_test_top100 <- na.omit(knn_test_top100)
  
  # start knn
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  knn_fit_model_top100 <- train(cities ~., data = knn_train_top100, method = "knn",
                                trControl=trctrl,
                                preProcess = c("center", "scale"),
                                tuneLength = 10)
  
  # watch output
  knn_fit_model_top100
  plot(knn_fit_model_top100)
  
  # test set prediction
  knn_pred_top100 <- predict(knn_fit_model_top100, newdata = knn_test_top100)
  knn_pred_top100
  
  # confusion matrix
  confusionMatrix(knn_pred_top100, knn_test_top100$cities)
  
  # top100 support vector machines ----------------------------------------
  
  svm_fit_model_top100 <- train(cities ~., data = knn_train_top100, method = "svmRadial", trControl = trctrl, preProcess = c("center","scale"), tuneLength = 10)
  
  # watch output
  svm_fit_model_top100
  plot(svm_fit_model_top100)
  
  # test set prediction
  svm_pred_top100 <- predict(svm_fit_model_top100, newdata = knn_test_top100)
  svm_pred_top100
  
  # confusion matrix
  confusionMatrix(svm_pred_top100, knn_test_top100$cities)
  
  # Top 100_1 analysis ----------------------------------------
  # Top 100_1 PCA and hierarchical clustering ----------------------------------------
  top100_1_TS <- as.matrix(top100_1[,-1])
  rownames(top100_1_TS) <- top100_1$species
  pca_top100_1 <- prcomp(t(top100_1_TS), scale=TRUE)
  biplot(prcomp(t(top100_1_TS)), main = "Top 10 species of each city")
  top100_1_TS2 <- t(top100_1_TS)
  
  top100_1_TS3 <- data.frame(top100_1_TS2)
  top100_1_TS3$cities <- as.factor(cities)
  
  #library(ggfortify)
  
  autoplot(pca_top100_1, main = "Top 10 species of each city")
  autoplot(pca_top100_1, data = top100_1_TS3, colour = 'cities', main = "Top 10 species of each city")
  autoplot(pca_top100_1, data = top100_1_TS3, colour = 'cities', main = "Top 10 species of each city",
           loadings = TRUE, loadings.colour = 'blue',
           loadings.label = TRUE, loadings.label.size = 3)
  
  # loadings of species
  pca_top100_1$rotation
  
  ### --- hierarchical clustering
  dd_top100_1 <- dist(scale(top100_1_TS2), method = "euclidean")
  
  hc_top100_1 <- hclust(dd_top100_1, method = "complete")
  plot(hc_top100_1, main = "Top 10 species of each city")
  
  # top 100_1 RF ----------------------------------------
  rf_train_top100_1 <- top100_1_TS3[-ind_test, ]
  rf_test_top100_1 <- top100_1_TS3[ind_test, ] # why do I get one row with only NAs?
  rf_test_top100_1 <- na.omit(rf_test_top100_1)
  
  # replace spaces and "-" by "_"
  #names(train_set_top100_1) <- str_replace_all(names(train_set_top100_1), c(" " = "_", "-" = "_"))
  #names(test_set_top100_1) <- str_replace_all(names(test_set_top100_1), c(" " = "_", "-" = "_"))
  
  # use caret with random forest as my model with 5 fold cross validation
  rf_model_top100_1<-train(cities~.,data=rf_train_top100_1,method="rf",
                           trControl=trainControl(method="cv",number=5),
                           prox=TRUE,allowParallel=TRUE)
  
  print(rf_model_top100_1)
  print(rf_model_top100_1$finalModel)
  plot(rf_model_top100_1)
  plot(rf_model_top100_1$finalModel)
  
  rf_pred_top100_1 <- predict(rf_model_top100_1, newdata = rf_test_top100_1)
  rf_pred_top100_1
  
  confusionMatrix(rf_pred_top100_1, rf_test_top100_1$cities)
  # top100_1 kNN ----------------------------------------
  knn_train_top100_1 <- top100_1_TS3[-ind_test, ]
  knn_test_top100_1 <- top100_1_TS3[ind_test, ]
  knn_test_top100_1 <- na.omit(knn_test_top100_1)
  
  # start knn
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  knn_fit_model_top100_1 <- train(cities ~., data = knn_train_top100_1, method = "knn",
                                  trControl=trctrl,
                                  preProcess = c("center", "scale"),
                                  tuneLength = 10)
  
  # watch output
  knn_fit_model_top100_1
  plot(knn_fit_model_top100_1)
  
  # test set prediction
  knn_pred_top100_1 <- predict(knn_fit_model_top100_1, newdata = knn_test_top100_1)
  knn_pred_top100_1
  
  # confusion matrix
  confusionMatrix(knn_pred_top100_1, knn_test_top100_1$cities)
  
  # top100_1 support vector machines ----------------------------------------
  
  svm_fit_model_top100_1 <- train(cities ~., data = knn_train_top100_1, method = "svmRadial", trControl = trctrl, preProcess = c("center","scale"), tuneLength = 10)
  
  # watch output
  svm_fit_model_top100_1
  plot(svm_fit_model_top100_1)
  
  # test set prediction
  svm_pred_top100_1 <- predict(svm_fit_model_top100_1, newdata = knn_test_top100_1)
  svm_pred_top100_1
  
  # confusion matrix
  confusionMatrix(svm_pred_top100_1, knn_test_top100_1$cities)
  
  # Top 100_2 analysis ----------------------------------------
  # Top 100_2 PCA and hierarchical clustering ----------------------------------------
  top100_2_TS <- as.matrix(top100_2[,-1])
  rownames(top100_2_TS) <- top100_2$species
  pca_top100_2 <- prcomp(t(top100_2_TS), scale=TRUE)
  biplot(prcomp(t(top100_2_TS)), main = "Top 10 species of each city")
  top100_2_TS2 <- t(top100_2_TS)
  
  top100_2_TS3 <- data.frame(top100_2_TS2)
  top100_2_TS3$cities <- as.factor(cities)
  
  #library(ggfortify)
  
  autoplot(pca_top100_2, main = "Top 10 species of each city")
  autoplot(pca_top100_2, data = top100_2_TS3, colour = 'cities', main = "Top 10 species of each city")
  autoplot(pca_top100_2, data = top100_2_TS3, colour = 'cities', main = "Top 10 species of each city",
           loadings = TRUE, loadings.colour = 'blue',
           loadings.label = TRUE, loadings.label.size = 3)
  
  # loadings of species
  pca_top100_2$rotation
  
  ### --- hierarchical clustering
  dd_top100_2 <- dist(scale(top100_2_TS2), method = "euclidean")
  
  hc_top100_2 <- hclust(dd_top100_2, method = "complete")
  plot(hc_top100_2, main = "Top 10 species of each city")
  
  # top 100_2 RF ----------------------------------------
  rf_train_top100_2 <- top100_2_TS3[-ind_test, ]
  rf_test_top100_2 <- top100_2_TS3[ind_test, ] # why do I get one row with only NAs?
  rf_test_top100_2 <- na.omit(rf_test_top100_2)
  
  # replace spaces and "-" by "_"
  #names(train_set_top100_2) <- str_replace_all(names(train_set_top100_2), c(" " = "_", "-" = "_"))
  #names(test_set_top100_2) <- str_replace_all(names(test_set_top100_2), c(" " = "_", "-" = "_"))
  
  # use caret with random forest as my model with 5 fold cross validation
  rf_model_top100_2<-train(cities~.,data=rf_train_top100_2,method="rf",
                           trControl=trainControl(method="cv",number=5),
                           prox=TRUE,allowParallel=TRUE)
  
  print(rf_model_top100_2)
  print(rf_model_top100_2$finalModel)
  plot(rf_model_top100_2)
  plot(rf_model_top100_2$finalModel)
  
  rf_pred_top100_2 <- predict(rf_model_top100_2, newdata = rf_test_top100_2)
  rf_pred_top100_2
  
  confusionMatrix(rf_pred_top100_2, rf_test_top100_2$cities)
  # top100_2 kNN ----------------------------------------
  knn_train_top100_2 <- top100_2_TS3[-ind_test, ]
  knn_test_top100_2 <- top100_2_TS3[ind_test, ]
  knn_test_top100_2 <- na.omit(knn_test_top100_2)
  
  # start knn
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  knn_fit_model_top100_2 <- train(cities ~., data = knn_train_top100_2, method = "knn",
                                  trControl=trctrl,
                                  preProcess = c("center", "scale"),
                                  tuneLength = 10)
  
  # watch output
  knn_fit_model_top100_2
  plot(knn_fit_model_top100_2)
  
  # test set prediction
  knn_pred_top100_2 <- predict(knn_fit_model_top100_2, newdata = knn_test_top100_2)
  knn_pred_top100_2
  
  # confusion matrix
  confusionMatrix(knn_pred_top100_2, knn_test_top100_2$cities)
  
  # top100_2 support v  ector machines ----------------------------------------
  
  svm_fit_model_top100_2 <- train(cities ~., data = knn_train_top100_2, method = "svmRadial", trControl = trctrl, preProcess = c("center","scale"), tuneLength = 10)
  
  # watch output
  svm_fit_model_top100_2
  plot(svm_fit_model_top100_2)
  
  # test set prediction
  svm_pred_top100_2 <- predict(svm_fit_model_top100_2, newdata = knn_test_top100_2)
  svm_pred_top100_2
  
  # confusion matrix
  confusionMatrix(svm_pred_top100_2, knn_test_top100_2$cities)
  
  # Results Section Thesis ----------------------------------------
  # RF Error rates table ----------------------------------------
  rf_error_table <- as.data.frame(rf_model_top10$finalModel$confusion[, 5])
  rf_error_table["Top 10_1"] <- as.data.frame(rf_model_top10_1$finalModel$confusion[, 5])
  rf_error_table["Top 10_2"] <- as.data.frame(rf_model_top10_2$finalModel$confusion[, 5])
  rf_error_table["Top 20"] <- as.data.frame(rf_model_top20$finalModel$confusion[, 5])
  rf_error_table["Top 20_1"] <- as.data.frame(rf_model_top20_1$finalModel$confusion[, 5])
  rf_error_table["Top 20_2"] <- as.data.frame(rf_model_top20_2$finalModel$confusion[, 5])
  rf_error_table["Top 30"] <- as.data.frame(rf_model_top30$finalModel$confusion[, 5])
  rf_error_table["Top 30_1"] <- as.data.frame(rf_model_top30_1$finalModel$confusion[, 5])
  rf_error_table["Top 30_2"] <- as.data.frame(rf_model_top30_2$finalModel$confusion[, 5])
  rf_error_table["Top 100"] <- as.data.frame(rf_model_top100$finalModel$confusion[, 5])
  rf_error_table["Top 100_1"] <- as.data.frame(rf_model_top100_1$finalModel$confusion[, 5])
  rf_error_table["Top 100_2"] <- as.data.framerre(rf_model_top100_2$finalModel$confusion[, 5])
  names(rf_error_table)[names(rf_error_table) == 'rf_model_top10$finalModel$confusion[, 5]'] <- 'Top 10'
  
  rf_mean_error <- c(colMeans(rf_error_table))
  rf_oob_error <- c(10.5, 12.9, 11.3, 8.87, 7.26, 11.3, 9.68, 8.87, 11.3, 12.1, 12.1, 12.1)
  rf_error_table <- rbind(rf_error_table, rf_mean_error, rf_oob_error)
  rownames(rf_error_table)[rownames(rf_error_table) == "5"] <- "mean error"
  rownames(rf_error_table)[rownames(rf_error_table) == "6"] <- "oob [%]"
  
  #write.csv(rf_error_table, file = "/home/rene/Dokumente/fh/thesis/writing/csd16_rf_error.csv")
  
  # bar plots to show how many species are in each group ----------------------------------------
  
  no_species <- c(nrow(top10), nrow(top10_1), nrow(top10_2),
                  nrow(top20), nrow(top20_1), nrow(top20_2),
                  nrow(top30), nrow(top30_1), nrow(top30_2),
                  nrow(top100), nrow(top100_1), nrow(top100_2)
  )
  
  bp <- barplot(no_species,
                names.arg = c("Top 10", "Top 10_1", "Top 10_2",
                              "Top 20", "Top 20_1", "Top 20_2",
                              "Top 30", "Top 30_1", "Top 30_2",
                              "Top 100", "Top 100_1", "Top 100_2"
                ),
                col = brewer.pal(n = 3, name = "Set2"),
                las = 2,
                ylab = "Number of Species",
                main = "Species per data set",
                ylim = c(0,250),
                cex.lab = 1,
                cex.main = 2
  )
  text(bp, no_species+8, paste(no_species), cex =1.5)
  # 1000 x 550
  
  PCA1 <- autoplot(pca_top10, data = top10_TS3, colour = 'cities', main = "Top 10 species of each city")
  PCA2 <- autoplot(pca_top10_1, data = top10_1_TS3, colour = 'cities', main = "Top 10_1 species of each city")
  PCA3 <- autoplot(pca_top10_2, data = top10_2_TS3, colour = 'cities', main = "Top 10_2 species of each city")
  PCA4 <- autoplot(pca_top20, data = top20_TS3, colour = 'cities', main = "Top 20 species of each city")
  PCA5 <- autoplot(pca_top20_1, data = top20_1_TS3, colour = 'cities', main = "Top 20_1 species of each city")
  PCA6 <- autoplot(pca_top20_2, data = top20_2_TS3, colour = 'cities', main = "Top 20_2 species of each city")
  PCA7 <- autoplot(pca_top30, data = top30_TS3, colour = 'cities', main = "Top 30 species of each city")
  PCA8 <- autoplot(pca_top30_1, data = top30_1_TS3, colour = 'cities', main = "Top 30_1 species of each city")
  PCA9 <- autoplot(pca_top30_2, data = top30_2_TS3, colour = 'cities', main = "Top 30_2 species of each city")
  PCA10 <- autoplot(pca_top100, data = top100_TS3, colour = 'cities', main = "Top 100 species of each city")
  PCA11 <- autoplot(pca_top100_1, data = top100_1_TS3, colour = 'cities', main = "Top 100_1 species of each city")
  PCA12 <- autoplot(pca_top100_2, data = top10_TS3, colour = 'cities', main = "Top 100_2 species of each city")
  
  ggarrange(PCA1, PCA2, PCA3,
            ncol = 3,
            nrow = 1)
  
  ggarrange(PCA4, PCA5, PCA6,
            ncol = 3,
            nrow = 1)
  
  ggarrange(PCA7, PCA8, PCA9,
            ncol = 3,
            nrow = 1)
  
  ggarrange(PCA10, PCA11, PCA12,
            ncol = 3,
            nrow = 1)
  
  # actual PCA plots in thesis
  ggarrange(PCA3, PCA8,
            ncol = 2,
            nrow = 1)
  # 1200 x 700
  # PCA plots ----------------------------------------
  
  PCA1 <- autoplot(pca_top10, data = top10_TS3, colour = 'cities', main = "Top 10 species of each city")
  PCA2 <- autoplot(pca_top10_1, data = top10_1_TS3, colour = 'cities', main = "Top 10_1 species of each city")
  PCA3 <- autoplot(pca_top10_2, data = top10_2_TS3, colour = 'cities', main = "Top 10_2 species of each city")
  PCA4 <- autoplot(pca_top20, data = top20_TS3, colour = 'cities', main = "Top 20 species of each city")
  PCA5 <- autoplot(pca_top20_1, data = top20_1_TS3, colour = 'cities', main = "Top 20_1 species of each city")
  PCA6 <- autoplot(pca_top20_2, data = top20_2_TS3, colour = 'cities', main = "Top 20_2 species of each city")
  PCA7 <- autoplot(pca_top30, data = top30_TS3, colour = 'cities', main = "Top 30 species of each city")
  PCA8 <- autoplot(pca_top30_1, data = top30_1_TS3, colour = 'cities', main = "Top 30_1 species of each city")
  PCA9 <- autoplot(pca_top30_2, data = top30_2_TS3, colour = 'cities', main = "Top 30_2 species of each city")
  PCA10 <- autoplot(pca_top100, data = top100_TS3, colour = 'cities', main = "Top 100 species of each city")
  PCA11 <- autoplot(pca_top100_1, data = top100_1_TS3, colour = 'cities', main = "Top 100_1 species of each city")
  PCA12 <- autoplot(pca_top100_2, data = top10_TS3, colour = 'cities', main = "Top 100_2 species of each city")
  
  ggarrange(PCA1, PCA2, PCA3,
            ncol = 3,
            nrow = 1)
  
  ggarrange(PCA4, PCA5, PCA6,
            ncol = 3,
            nrow = 1)
  
  ggarrange(PCA10, PCA11, PCA12,
            ncol = 3,
            nrow = 1)
  
  # actual PCA plots in thesis
  ggarrange(PCA3, PCA8,
            ncol = 2,
            nrow = 1)
  # 1200 x 550
  # 1200 x 700 (_2)
  # RF  plots ----------------------------------------
  # top10_1 and top20_1 -> least and highest OOB
  
  # RF model
  # 1200 x 700
  plot(rf_model_top10_1$finalModel, main = "Random forest model for the \"Top 10_1\" dataset", cex.main = 2, col =rf_colors, lwd = 1.5)
  legend("topright", colnames(rf_model_top10$finalModel$err.rate),col=rf_colors,cex=0.8, fill = rf_colors)
  
  plot(rf_model_top20_1$finalModel, main = "Random forest model for the \"Top 20_1\" dataset", cex.main = 2, col =rf_colors, lwd = 1.5)
  legend("topright", colnames(rf_model_top10$finalModel$err.rate),col=rf_colors,cex=0.8, fill = rf_colors)
  
  # VarImp
  # 1200 x 500
  varimp_top10_1 <- plot(varImp(rf_model_top10_1, scale = FALSE), top = 10, main = "Variable importance of \"Top 10_1\"", cex.main = 2, xlim = c(0, 13))
  varimp_top20_1 <- plot(varImp(rf_model_top20_1, scale = FALSE), top = 10, main = "Variable importance of \"Top 20_1\"", cex.main = 2, xlim = c(0,8))
  
  ggarrange(varimp_top10_1, varimp_top20_1,
            ncol = 2,
            nrow = 1)
  
  # KNN plots ----------------------------------------
  
  par(mfrow = c(2,2))
  # 1300 x 800
  plot(knn_fit_model_top10$results$k, knn_fit_model_top10$results$Accuracy, type = "b", lty =1, lwd = 1.5, col = "red", pch = 1,
       main = "KNN model accuracy of top 10 species",
       xlab = "No. of Neighbours",
       ylab = "Accuracy (Repeated Cross-Validation)",
       xlim = c(5, 25),
       ylim = c(0.45, 0.8))
  lines(knn_fit_model_top10_1$results$k, knn_fit_model_top10_1$results$Accuracy, type = "b", lty =2, lwd = 1.5, col = "green", pch = 2)
  lines(knn_fit_model_top10_2$results$k, knn_fit_model_top10_2$results$Accuracy, type = "b", lty =3, lwd = 1.5, col = "blue", pch = 3)
  legend("topright", legend = c("Top 10", "Top 10_1", "Top 10_2"),
         col = c("red","green", "blue"), lty = 1:3, cex = 0.8)
  
  plot(knn_fit_model_top20$results$k, knn_fit_model_top20$results$Accuracy, type = "b", lty =1, lwd = 1.5, col = "red", pch = 1,
       main = "KNN model accuracy of top 20 species",
       xlab = "No. of Neighbours",
       ylab = "Accuracy (Repeated Cross-Validation)",
       xlim = c(5, 25),
       ylim = c(0.45, 0.8))
  lines(knn_fit_model_top20_1$results$k, knn_fit_model_top20_1$results$Accuracy, type = "b", lty =2, lwd = 1.5, col = "green", pch = 2)
  lines(knn_fit_model_top20_2$results$k, knn_fit_model_top20_2$results$Accuracy, type = "b", lty =3, lwd = 1.5, col = "blue", pch = 3)
  legend("topright", legend = c("Top 20", "Top 20_1", "Top 20_2"),
         col = c("red","green", "blue"), lty = 1:3, cex = 0.6)
  
  plot(knn_fit_model_top30$results$k, knn_fit_model_top30$results$Accuracy, type = "b", lty =1, lwd = 1.5, col = "red", pch = 1,
       main = "KNN model accuracy of top 30 species",
       xlab = "No. of Neighbours",
       ylab = "Accuracy (Repeated Cross-Validation)",
       xlim = c(5, 25),
       ylim = c(0.45, 0.8))
  lines(knn_fit_model_top30_1$results$k, knn_fit_model_top30_1$results$Accuracy, type = "b", lty =2, lwd = 1.5, col = "green", pch = 2)
  lines(knn_fit_model_top30_2$results$k, knn_fit_model_top30_2$results$Accuracy, type = "b", lty =3, lwd = 1.5, col = "blue", pch = 3)
  legend("topright", legend = c("Top 30", "Top 30_1", "Top 30_2"),
         col = c("red","green", "blue"), lty = 1:3, cex = 0.6)
  
  plot(knn_fit_model_top100$results$k, knn_fit_model_top100$results$Accuracy, type = "b", lty =1, lwd = 1.5, col = "red", pch = 1,
       main = "KNN model accuracy of top 100 species",
       xlab = "No. of Neighbours",
       ylab = "Accuracy (Repeated Cross-Validation)",
       xlim = c(5, 25),
       ylim = c(0.45, 0.8))
  lines(knn_fit_model_top100_1$results$k, knn_fit_model_top100_1$results$Accuracy, type = "b", lty =2, lwd = 1.5, col = "green", pch = 2)
  lines(knn_fit_model_top100_2$results$k, knn_fit_model_top100_2$results$Accuracy, type = "b", lty =3, lwd = 1.5, col = "blue", pch = 3)
  legend("topright", legend = c("Top 100", "Top 100_1", "Top 100_2"),
         col = c("red","green", "blue"), lty = 1:3, cex = 0.8)
  
  par(mfrow = c(1,1))
  
  # SVM plots ----------------------------------------
  
  
  par(mfrow = c(2,2))
  # 1300 x 800
  plot(svm_fit_model_top10$results$C, svm_fit_model_top10$results$Accuracy, type = "b", lty =1, lwd = 1.5, col = "red", pch = 1,
       main = "SVM model accuracy of top 10 species",
       xlab = "Cost",
         ylab = "Accuracy (Repeated Cross-Validation)",
       xlim = c(0, 140),
       ylim = c(0.45, 0.8))
  lines(svm_fit_model_top10_1$results$C, svm_fit_model_top10_1$results$Accuracy, type = "b", lty =2, lwd = 1.5, col = "green", pch = 2)
  lines(svm_fit_model_top10_2$results$C, svm_fit_model_top10_2$results$Accuracy, type = "b", lty =3, lwd = 1.5, col = "blue", pch = 3)
  legend("bottomright", legend = c("Top 10", "Top 10_1", "Top 10_2"),
         col = c("red","green", "blue"), lty = 1:3, cex = 0.5)
  
  plot(svm_fit_model_top20$results$C, svm_fit_model_top20$results$Accuracy, type = "b", lty =1, lwd = 1.5, col = "red", pch = 1,
       main = "SVM model accuracy of top 20 species",
       xlab = "Cost",
       ylab = "Accuracy (Repeated Cross-Validation)",
       xlim = c(0, 140),
       ylim = c(0.45, 0.8))
  lines(svm_fit_model_top20_1$results$C, svm_fit_model_top20_1$results$Accuracy, type = "b", lty =2, lwd = 1.5, col = "green", pch = 2)
  lines(svm_fit_model_top20_2$results$C, svm_fit_model_top20_2$results$Accuracy, type = "b", lty =3, lwd = 1.5, col = "blue", pch = 3)
  legend("bottomright", legend = c("Top 20", "Top 20_1", "Top 20_2"),
         col = c("red","green", "blue"), lty = 1:3, cex = 0.8)
  
  plot(svm_fit_model_top30$results$C, svm_fit_model_top30$results$Accuracy, type = "b", lty =1, lwd = 1.5, col = "red", pch = 1,
       main = "SVM model accuracy of top 30 species",
       xlab = "Cost",
       ylab = "Accuracy (Repeated Cross-Validation)",
       xlim = c(0, 140),
       ylim = c(0.45, 0.8))
  lines(svm_fit_model_top30_1$results$C, svm_fit_model_top30_1$results$Accuracy, type = "b", lty =2, lwd = 1.5, col = "green", pch = 2)
  lines(svm_fit_model_top30_2$results$C, svm_fit_model_top30_2$results$Accuracy, type = "b", lty =3, lwd = 1.5, col = "blue", pch = 3)
  legend("bottomright", legend = c("Top 30", "Top 30_1", "Top 30_2"),
         col = c("red","green", "blue"), lty = 1:3, cex = 0.8)
  
  plot(svm_fit_model_top100$results$C, svm_fit_model_top100$results$Accuracy, type = "b", lty =1, lwd = 1.5, col = "red", pch = 1,
       main = "SVM model accuracy of top 100 species",
       xlab = "Cost",
       ylab = "Accuracy (Repeated Cross-Validation)",
       xlim = c(0, 140),
       ylim = c(0.45, 0.8))
  lines(svm_fit_model_top100_1$results$C, svm_fit_model_top100_1$results$Accuracy, type = "b", lty =2, lwd = 1.5, col = "green", pch = 2)
  lines(svm_fit_model_top100_2$results$C, svm_fit_model_top100_2$results$Accuracy, type = "b", lty =3, lwd = 1.5, col = "blue", pch = 3)
  legend("bottomright", legend = c("Top 100", "Top 100_1", "Top 100_2"),
         col = c("red","green", "blue"), lty = 1:3, cex = 0.8)
  
  par(mfrow = c(1,1))
  
  # table extract ----------------------------------------
  
  write.csv(top10, "extract/processed_data/csd16/top10.csv")
  write.csv(top10_1, "extract/processed_data/csd16/top10_1.csv")
  write.csv(top10_2, "extract/processed_data/csd16/top10_2.csv")
  
  write.csv(top20, "extract/processed_data/csd16/top20.csv")
  write.csv(top20_1, "extract/processed_data/csd16/top20_1.csv")
  write.csv(top20_2, "extract/processed_data/csd16/top20_2.csv")
  
  write.csv(top30, "extract/processed_data/csd16/top30.csv")
  write.csv(top30_1, "extract/processed_data/csd16/top30_1.csv")
  write.csv(top30_2, "extract/processed_data/csd16/top30_2.csv")
  
  write.csv(top100, "extract/processed_data/csd16/top100.csv")
  write.csv(top100_1, "extract/processed_data/csd16/top100_1.csv")
  write.csv(top100_2, "extract/processed_data/csd16/top100_2.csv")
  
  write.csv(prokarya_all, "extract/processed_data/csd16/bacteria")
  
  # hc ----------------------------------------
  #1200x700
  tdf <- as.data.frame(top10_TS2)
  top10_dend <- as.dendrogram(hclust(dist(tdf)))
  c_type <- rep("Other", length(rownames(tdf)))
  is_x <- grepl("ILO", rownames(tdf))
  c_type[is_x] <- "Ilorin"
  is_x <- grepl("DOH", rownames(tdf))
  c_type[is_x] <- "Doha"
  is_x <- grepl("NYC", rownames(tdf))
  c_type[is_x] <- "New York"
  is_x <- grepl("FAI", rownames(tdf))
  c_type[is_x] <- "Fairbanks"
  c_type <- factor(c_type)
  n_c_types <- length(unique(c_type))
  #col_c_type <- colorspace::rainbow_hcl(n_c_types, c = 90, l = 60)[c_type]
  col_c_type <- c("blue", "red", "darkgreen", "yellow")[c_type]
  # color labels by car company:
  labels_colors(top10_dend) <- col_c_type[order.dendrogram(top10_dend)]
  # color branches based on cutting the tree into 4 clusters:
  #top10_dend <- color_branches(top10_dend)
  
  ### plots
  par(mar = c(12, 4, 1, 1))
  plot(top10_dend, main = "Hierarchical clustering Top 10")
  colored_bars(cbind(col_c_type), top10_dend,
               rowLabels = "City"
  )
  legend("topright", legend = c("Ilorin", "Doha", "Fairbanks", "New York"), fill = c("darkgreen", "blue", "red", "yellow"))
  
  
  col_c_type <- c("red", "blue", "gray", "darkgreen")