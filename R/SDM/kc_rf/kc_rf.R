library(ranger)
library(rsample)
library(spatialsample)
library(dplyr)
library(pROC)
library(sf)
library(ggplot2)
library(pscl)
library(corrplot)
library(xts)
library(ISOweek)
library(glmnet)
library(raster)
library(PerformanceAnalytics)
set.seed(123)



# Spatial CV

sf <- read_sf("Data/3_clusterExport/clusterSHP/clusterPolys.shp")

clusters <- spatial_block_cv(sf, method = "random", n = 10, relevant_only = T, v = 5)

splits_df <- data.frame()
for (i in 1:5) {
  new_df <- assessment(clusters$splits[[i]])
  new_df$fold <- i
  splits_df <- rbind(splits_df, new_df)  # Bind all points and fold IDs together
}

splits_df <- st_drop_geometry(splits_df)
splits_df <- distinct(splits_df, clust, .keep_all = TRUE)


MIRPIR <- read.csv("Data/9_joinGEE/wnvMIRPIR1500LagWeeks0.csv")

MIRPIR_meanTemp <- MIRPIR %>%
  mutate(Temp.Mean.K = rowMeans(dplyr::select(., Temp.Max.K, Temp.Min.K), na.rm = TRUE)) %>%
  filter(!is.na(Temp.Max.K) & !is.na(Temp.Min.K)) %>%
  dplyr::select(-Temp.Max.K, -Temp.Min.K) %>%
  merge(splits_df, by = "clust") %>%
  mutate(MIRAll = ifelse(MIRAll > 0, 1, 0))

form <- MIRAll ~ chirpsMean + Humidity.Max + Humidity.Min + Specific.Humidity +
  Temp.Mean.K + VPD.kPa + Precip.Mean.mm...day + EVI + NDVI + jrcStandingWater +
  irrWater + eddi14d + eddi180d + eddi1y + eddi270d + eddi2y + eddi30d + eddi5y +
  eddi90d + pdsi + spei14d + spei180d + spei1y + spei270d + spei2y + spei30d + spei5y +
  spei90d + spi14d + spi180d + spi1y + spi270d + spi2y + spi30d + spi90d + z


mirpir_split <- initial_split(MIRPIR_meanTemp, prop = 0.8, strata = "clust")
mirpir_train <- training(mirpir_split)
mirpir_test <- testing(mirpir_split)

mirpir <- mirpir_train[c(c("MIRAll", "fold", "chirpsMean" , "Humidity.Max" ,
                           "Humidity.Min" , "Specific.Humidity" , "Temp.Mean.K" ,
                           "VPD.kPa" , "Precip.Mean.mm...day" ,"EVI", "NDVI" ,
                           "jrcStandingWater" , "irrWater" , "eddi14d" , "eddi180d" ,
                           "eddi1y" , "eddi270d" , "eddi2y" , "eddi30d" , "eddi5y" ,
                           "eddi90d" , "pdsi" , "spei14d" , "spei180d" , "spei1y" ,
                           "spei270d" , "spei2y" , "spei30d" , "spei5y" , "spei90d" ,
                           "spi14d" , "spi180d" , "spi1y" , "spi270d" , "spi2y" ,
                           "spi30d", "spi90d" , "z"))]


mirpir_final <- mirpir_test[c(c("MIRAll", "fold", "chirpsMean" , "Humidity.Max" ,
                                "Humidity.Min" , "Specific.Humidity" , "Temp.Mean.K" ,
                                "VPD.kPa" , "Precip.Mean.mm...day" ,"EVI", "NDVI" ,
                                "jrcStandingWater" , "irrWater" , "eddi14d" , "eddi180d" ,
                                "eddi1y" , "eddi270d" , "eddi2y" , "eddi30d" , "eddi5y" ,
                                "eddi90d" , "pdsi" , "spei14d" , "spei180d" , "spei1y" ,
                                "spei270d" , "spei2y" , "spei30d" , "spei5y" , "spei90d" ,
                                "spi14d" , "spi180d" , "spi1y" , "spi270d" , "spi2y" ,
                                "spi30d", "spi90d" , "z"))]

final_mirpir <- MIRPIR_meanTemp[c(c("MIRAll", "fold", "chirpsMean" , "Humidity.Max" ,
                                    "Humidity.Min" , "Specific.Humidity" , "Temp.Mean.K" ,
                                    "VPD.kPa" , "Precip.Mean.mm...day" ,"EVI", "NDVI" ,
                                    "jrcStandingWater" , "irrWater" , "eddi14d" , "eddi180d" ,
                                    "eddi1y" , "eddi270d" , "eddi2y" , "eddi30d" , "eddi5y" ,
                                    "eddi90d" , "pdsi" , "spei14d" , "spei180d" , "spei1y" ,
                                    "spei270d" , "spei2y" , "spei30d" , "spei5y" , "spei90d" ,
                                    "spi14d" , "spi180d" , "spi1y" , "spi270d" , "spi2y" ,
                                    "spi30d", "spi90d" , "z"))]

# Tuning grid

rf_performance <- data.frame(model  =  rep("RF", 5),
                             fold_id  =  1:5,
                             auc  =  rep(NA, 5),
                             presence  =  rep(NA, 5),
                             background  =  rep(NA, 5))

hypergrid_final <- data.frame(mtry  =  rep(NA,5),
                              node_size    =  rep(NA, 5),
                              sample_size  =  rep(NA, 5))





for(i  in  1:5){

  train  <-  mirpir[mirpir$fold !=  1, ];  train  <-  train[-2]
  test  <-  mirpir[mirpir$fold ==  1, ];  test  <-  test[-2]

  train_complete  <-  train[complete.cases(train), ]
  test_complete  <-  test[complete.cases(test), ]

  hyper_grid  <-  expand.grid(
    mtry =  seq(1, 10, by  =  1),
    node_size =  seq(1,4, by  =  1),
    sample_size  =  c(.6, .70, .80),
    OOB_RMSE =  0
  )

  for(j  in  1:nrow(hyper_grid)){

    model  <-  ranger(
      formula  =  MIRAll  ~  .,
      data  =  train_complete,
      num.trees  =  2000,
      mtry  =  hyper_grid$mtry[j],
      min.node.size  =  hyper_grid$node_size[j],
      sample.fraction  =  hyper_grid$sample_size[j],
      probability  =  TRUE,
      replace = TRUE,
      splitrule = "extratrees",
      seed  =  123
    )

    hyper_grid$OOB_RMSE[j]  <-  sqrt(model$prediction.error)
  }

  hyper_grid2  <-  hyper_grid  %>%
    dplyr::arrange(OOB_RMSE)

  #train  model
  train_model  <-  ranger(
    formula  =  MIRAll  ~  .,
    data  =  train_complete,
    num.trees  =  2000,
    mtry  =  hyper_grid2$mtry[1],
    min.node.size  =  hyper_grid2$node_size[1],
    sample.fraction  =  hyper_grid2$sample_size[1],
    probability  =  TRUE,
    replace = TRUE,
    splitrule = "hellinger",
    seed  =  123)

  #save  model  performance  results
  pred0  <-  predict(train_model, data=test_complete);  pred  <-  pred0$predictions[,1]
  auc  <-  pROC::roc(response=test_complete[,"MIRAll"], predictor=pred, levels=c(0, 1), auc  =  TRUE)
  rf_performance[i, "auc"]  <-  auc$auc
  rf_performance[i, "presence"] <- nrow(subset(test, MIRAll  ==  1))
  rf_performance[i, "background"] <-nrow(subset(test, MIRAll  ==  0))

  hypergrid_final[i, "mtry"]  <-  hyper_grid2$mtry[1]
  hypergrid_final[i, "node_size"]  <-  hyper_grid2$node_size[1]
  hypergrid_final[i, "sample_size"]  <-  hyper_grid2$sample_size[1]

}

saveRDS(hypergrid_final, file = "Data/SDM/kc_rf/hypergrid_final")
saveRDS(rf_performance, file = "Data/SDM/kc_rf/rf_performance")


final_model  <-  ranger(
  formula  =  MIRAll  ~  .,
  data  =  mirpir[complete.cases(mirpir), -2],
  num.trees  =  2000 ,
  mtry  =  1,
  min.node.size  =  1,
  sample.fraction  =  0.6,
  probability  =  TRUE,
  replace = TRUE,
  splitrule = "hellinger",
  importance  =  'permutation',
  seed  =  123)

pred0  <-  predict(final_model, data=mirpir_final[complete.cases(mirpir_final), -2]);  pred  <-  pred0$predictions[,1]
auc  <-  pROC::roc(response=mirpir_final[complete.cases(mirpir_final), -2][,"MIRAll"], predictor=pred, levels=c(0, 1), auc  =  TRUE)

plot(auc)

saveRDS(final_model, file = "Data/SDM/rf/final_model")


# Variable Importance
permutation_importance  <-  data.frame(variable  =  rownames(as.data.frame(final_model$variable.importance)),
                                       importance  = as.vector(final_model$variable.importance))

#plot  importance
ggplot(permutation_importance, aes(x  =  reorder(variable, -importance), y  =  importance))  +
  geom_bar(stat="identity")  +
  ggtitle("Permutation Importance")  +
  ylab("importance (change in model error)") +
  coord_flip()  +
  theme_classic()


# PDPs
var_names  <-  names(final_mirpir[complete.cases(final_mirpir),  -c(1, 2)])

pd_df  =  data.frame(matrix(vector(), 0, 3, dimnames=list(c(), c('variable', 'value', 'yhat'))),
                     row.names  =  NULL, stringsAsFactors=F)


for  (i  in  1:length(var_names))  {

  output  <-  as.data.frame(pdp::partial(final_model, pred.var  =  var_names[i], prob  =  TRUE, train  = final_mirpir[complete.cases(final_mirpir), -2]))

  loop_df  <-  data.frame(variable = rep(var_names[i], length(output[[1]])),
                          value  =  output[[1]],
                          yhat  =  output[[2]])

  pd_df  <-  rbind(pd_df, loop_df)
}

saveRDS(pd_df, file = "Data/SDM/rf/kc_rd_pd")


kc_pd_df <- readRDS("Data/SDM/rf/kc_pd_df_v2")
output_directory <- "Data/SDM/rf/kc_pdps"
pdp_plots <- list()

if (!dir.exists(output_directory)) {
  dir.create(output_directory)
}

for (variable_name in unique(kc_pd_df$variable)) {

  variable_data <- kc_pd_df[kc_pd_df$variable == variable_name, ]

  pdp_plot <- ggplot(variable_data, aes(x = value, y = yhat)) +
    geom_smooth() +
    ylab("probability") +
    ggtitle(paste("Partial Dependence Plot for", variable_name)) +
    theme_bw(base_size = 14)

  pdp_plots[[variable_name]] <- pdp_plot

  output_file <- file.path(output_directory, paste0("pdp_", variable_name, ".png"))
  ggsave(output_file, pdp_plot, width = 8, height = 6)
}

for (i in 1:length(pdp_plots)) {
  print(pdp_plots[[i]])
}