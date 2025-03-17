args = commandArgs(TRUE)

outdir=as.character(args[1])
value=as.character(args[2])
m=as.numeric(args[3])
model_name=as.character(args[4])

impute_outlier_data = function(outdir,value,m,model_name ) {
  library(mice)
  set.seed(123)
  filename = file.path(outdir, 
                       paste("data",
                             model_name,
                             "rda", 
                             sep = "."))
  load(filename)
  if(value == "delta") {
    db = df$df$delta$df
    idx = df$outliers$delta$idx
    rm.subs = df$outliers$delta$rmsubs.miss.data
  } else if (value == "mean") {
    db = df$df$mean$df
    idx = df$outliers$mean$idx
    rm.subs = df$outliers$mean$rmsubs.miss.data
  }
  names(db) = paste0("V", 1:length(names(db)))
  db[!idx] = NaN
  db = db[rm.subs,]
  imputed = mice(db, m=m)
  
  load(filename)
  if(value == "delta") {
    df$df$imputed.delta = imputed
  } else if (value == "mean") {
    df$df$imputed.mean = imputed
  }
  save(df, file = filename)
}


impute_outlier_data(outdir,value,m, model_name) 

