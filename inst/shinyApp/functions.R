

ui_tabOptions <- function(options, algo){
  
  docs <- paste0("icon ('book'), a('Documentation for ", algo, "', href = '", link_documentation(algo), "')")
  
  name_algo <- grep(algo, options@models, value = T)
  if(algo == "RF"){name_algo <- name_algo[1]}
  
  # names_options <- names(options@options[[name_algo]]@args.values$`_allData_allRun`)
  # 
  # ## Remove formula and other problems
  # to_remove <- c("formula", "data", "weights", "...", "subset", "family", "strata", 
  #                "resp.var", "expl.var", "new.env", "glm", "na.action")
  # names_options <- names_options[!names_options %in% to_remove]
  
  names_options <- OptionsNames(algo)
  
  vec <- c()
  for (i in names_options){
    value <- options@options[[name_algo]]@args.values$`_allData_allRun`[[i]]
    if (missing(value)){
      value <- ""
    }
    if (is.list(value)){
      l.names <- names(value)
      for (l in l.names){
        vec <- c(vec, paste0("textInput('", i, "_", l, "', '", i, "$", l, "', value = '", value[[l]], "')"))
      }
    } else {
      vec <- c(vec, paste0("textInput('", i, "', '", i, "', value = '", value, "')"))
    }
  }
  
  texte <- paste0( vec, collapse = ", \t")
  texte <- paste0("column(6, br(), ", docs, ",br() , br(), \t", texte, ")")
  return(texte)
}


hide_several_tabs <- function(inputId, vec){
  for (i in 1:length(vec)){
    hideTab(inputId = inputId, target = vec[i])
  }
}


link_documentation <- function(algo){
  switch(algo, 
         "ANN" = "https://cran.r-project.org/web/packages/nnet/nnet.pdf",
         "CTA" = "https://cran.r-project.org/web/packages/rpart/rpart.pdf",
         "FDA" = "https://cran.r-project.org/web/packages/mda/mda.pdf",
         "GAM" = "https://cran.r-project.org/web/packages/mgcv/mgcv.pdf",
         "GBM" = "https://cran.r-project.org/web/packages/gbm/gbm.pdf",
         "GLM" = "https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/glm",
         "MARS" = "https://cran.r-project.org/web//packages/earth/earth.pdf",
         "MAXENT" = "https://biomodhub.github.io/biomod2/reference/bm_ModelingOptions.html",
         "MAXNET" = "https://biomodhub.github.io/biomod2/reference/bm_ModelingOptions.html",
         "RF" = "https://cran.r-project.org/web/packages/randomForest/randomForest.pdf",
         "RFd" = "https://cran.r-project.org/web/packages/randomForest/randomForest.pdf",
         "SRE" = "https://biomodhub.github.io/biomod2/reference/bm_SRE.html",
         "XGBOOST" = "https://cran.r-project.org/web/packages/xgboost/xgboost.pdf"
         )
}

get_inputs <- function(base, models){
  texte <- "list_inputs <- list("
  for (m in models){
    names_options <- OptionsNames(m)
    list_algo <- c()
    for (i in names_options){
      list_algo <- c(list_algo, paste0(i , "= input$", i))
    }
    list_algo <- paste0(list_algo, collapse = ", ")
    list_algo <- paste0("list(",list_algo, ")")
    
    texte <- paste0(texte, m , " = ", list_algo, ",")
  }
  texte <- substr(texte, 1, nchar(texte)-1)
  texte <- paste0(texte, ")")
  return(texte)
}




create_options <- function(list_inputs, base, datasets, models, data.type, strategy, name_obj){
  
  ##Create user.val 
  user.val <- list()
  for (algo in models){
    name_algo <- grep(algo, base@models, value = T)
    if(algo == "RF"){name_algo <- name_algo[1]}
    
    opt <- list_inputs[[algo]]
    for (i in names(opt)){
      ### cas des list 
      if (grepl("_", i)){
        level1 <- unlist(strsplit(i, "_"))[1]
        level2 <- unlist(strsplit(i, "_"))[2]

        if (is.null(opt[[level1]])){opt[[level1]] <- list()}

        value <- base@options[[name_algo]]@args.values$`_allData_allRun`[[level1]][[level2]]
        if (missing(value)){
          value <- NULL
        }
        new <- opt[[i]]
        class(new) <- class(value)
        
        opt <- opt[names(opt) != i]
        opt[[level1]][[level2]] <- new
        
      } else {
        value <- base@options[[name_algo]]@args.values$`_allData_allRun`[[i]]
        if (missing(value)){
          value <- NULL
        }
        new <- opt[[i]]
        class(new) <- class(value)
        opt[[i]] <- new
      }
    }
    
    # user.val[[name_algo]] <- list()
    # user.val[[name_algo]][[datasets]] <- opt
    user.val[[name_algo]] <- list("for_all_datasets" = opt)
  }
  
  bm_options <- biomod2::bm_ModelingOptions(data.type = data.type,
                                            models = models,
                                            strategy = "user.defined",
                                            user.base = strategy,
                                            user.val = user.val)
  
  cat_code_used(user.val, data.type, strategy, models, name_obj)
  
  return(bm_options)
}


## Get the names of the options for each algorithm 
## 2 options : 
## une fonction 

OptionsNames_byAlgo <- function(base, algo){
  
  name_algo <- grep(algo, base@models, value = T)
  if(algo == "RF"){name_algo <- name_algo[1]}
  names_options <- names(base@options[[name_algo]]@args.values$`_allData_allRun`)
  
  ## Remove formula and other problems
  to_remove <- c("formula", "data", "weights", "...", "subset", "family", "strata", 
                 "resp.var", "expl.var", "new.env", "glm")
  names_options <- names_options[!names_options %in% to_remove]
  
  for (i in names_options){
    value <- base@options[[name_algo]]@args.values$`_allData_allRun`[[i]]
    if (missing(value)){
      value <- ""
    }
    if (is.list(value)){
      names_options <- names_options[!names_options == i]
      l.names <- names(value)
      for (l in l.names){
        names_options <- c(names_options, paste0(i,"_",l))
        
      }
    }
  }
  return(names_options)
}


## Ou un switch (Pour l'instant ça permet un meilleur contrôle de ce qu'on veut mettre ou non)

OptionsNames <- function(algo){
  switch(algo, 
         "ANN" = c("size", "decay", "rang", "maxit", "trace"),
         "CTA" = c("method", "parms", "cost",
                   "control_xval", "control_minbucket", "control_minsplit", "control_cp", "control_maxdepth",
                   "model", "x", "y"),
         "FDA" = c("theta", "eps", "method"),
         "GAM" = c("method", "optimizer", "scale", "select", "gamma",
                   "control_epsilon", "control_trace", "control_maxit"),
         "GBM" = c( "distribution", "n.trees", "interaction.depth", "n.minobsinnode", "shrinkage", "bag.fraction",
                    "train.fraction", "cv.folds", "keep.data", "verbose", "n.cores"),
         "GLM" = c( "etastart", "mustart", "offset", "method",         
                   "control_epsilon", "control_maxit", "control_trace",
                   "singular.ok", "model", "x", "y"),
         "MARS" = c("ncross", "penalty", "thresh", "pmethod"),
         "MAXENT" = c("path_to_maxent.jar", "memory_allocated", "background_data_dir", "visible", "linear", 
                      "quadratic", "product", "threshold", "hinge",  "lq2lqptthreshold",  "l2lqthreshold",
                      "hingethreshold", "beta_threshold", "beta_categorical", "beta_lqp", "beta_hinge", "betamultiplier",
                      "defaultprevalence", "togglelayerselected", "maximumbackground", "maximumiterations", "convergencethreshold",
                      "autofeature", "jackknife", "writeclampgrid", "writemess", "logfile", "verbose"),
         "MAXNET" = c("regmult", "regfun", "addsamplestobackground"),
         "RF" = c("type", "mtry", "ntree", "nodesize"),
         "RFd" = c("type", "mtry", "ntree", "nodesize"),
         "SRE" = c("quant", "do.extrem"),
         "XGBOOST" = c("objective", "nrounds", "nthread", "params_max_depth", "params_eta",
                       "callbacks", "missing", "save_name", "verbose", "print_every_n")
  )
}


## Je pense quand même que cette app pourrait être bien plus simple


## ¨Print used code
cat_code_used <- function(user.val, data.type, base, models, name_obj){
  texte_1 <- "To obtain these options, you could run this code :\n\n"
  texte_2 <- paste0("user.val <-", print_list(user.val), ")\n\n")
  texte_3 <- paste0(name_obj, " <- bm_ModelingOptions(", 
                    "\n\t data.type = '", data.type, "',",
                    "\n\t models = c('", paste0(models, collapse = "', '"), "'),",
                    "\n\t strategy = 'user.defined',",
                    "\n\t user.base = '", base, "',",
                    "\n\t user.val = user.val)"
  )
  cat(texte_1, texte_2, texte_3)
}

print_list <- function(l, n = 1){
  texte <- "list("
  for (i in 1:length(l)){
    if(is.list(l[[i]])){
      texte <- paste0(texte, "\n", paste(rep("\t",n), collapse = ""), names(l)[i],  " = ", print_list(l[[i]], n = n+1), ", ")
    } else {
      if(!is.null(l[[i]]) & l[[i]] != "" & !is.na(l[[i]])){
        texte <- paste0(texte, "\n", paste(rep("\t",n), collapse = ""), names(l)[i],  " = ", l[[i]], ", ")
      }
    }
  }
  texte <- substr(texte, 1, nchar(texte) - 2)
  texte <- paste0(texte, ")")
  texte
}



### Fonction for gestion biomod.formated.data and CVtable

namesObjectsByType <- function(type){
  names_obj <- objects(envir = globalenv())
  text  <- paste0(" inherits(", names_obj, ", '", type, "' )", sep = "")
  yes <- sapply(text , FUN = function(x){eval(parse(text = x), envir = globalenv())})
  names_type <- names_obj[yes]
  if(rlang::is_empty(names_type)) {names_type <- paste0("No ", type, " object available")}
  names_type
}

get_PA_datasets <- function(bm.format){
  if (inherits(bm.format, "BIOMOD.formated.data.PA")){
    choices <- c("All the same", paste0("PA ", 1:ncol(bm.format@PA.table)), "allData")  
  } else {
    choices <- c("allData")
  }
  choices 
}

get_CV_datasets <- function(bm.format, cvtable){
  ##Check if number of PA in bm.format and cvtable are coherent ? 
  if (is.null(cvtable)){
    names_run <- c()
  } else {
    names_run <- unlist(strsplit(colnames(cvtable), split = "_"))
    names_run <- names_run[grep(names_run, pattern = "RUN")]
  }
  return(c("All the same", names_run))
}


## Update user.val for one datasets 

update_userval <- function(user.val, whichPA, whichCV, list_inputs,
                           names_PA, names_CV, models, base){# list_inputs, base, datasets, models, data.type, strategy, name_obj
  ## Which dataset
  if (whichPA == "All the same"){
    if (whichCV == "All the same"){
      datasets <- "for_all_datasets"
    } else {
      datasets <- paste0("_", names_PA, "_", whichCV)
    }
  } else {
    if (whichCV == "All the same"){
      datasets <- paste0("_", whichPA, "_", names_CV[-1])
    } else {
      datasets <- paste0("_", whichPA, "_", whichCV)
    }
  }
  
  ## 
  for (algo in models){
    name_algo <- grep(algo, base@models, value = T)
    if(algo == "RF"){name_algo <- name_algo[1]}
    
    opt <- list_inputs[[algo]]
    for (i in names(opt)){
      ### cas des list 
      if (grepl("_", i)){
        level1 <- unlist(strsplit(i, "_"))[1]
        level2 <- unlist(strsplit(i, "_"))[2]
        
        if (is.null(opt[[level1]])){opt[[level1]] <- list()}
        
        value <- base@options[[name_algo]]@args.values$`_allData_allRun`[[level1]][[level2]]
        if (missing(value)){
          value <- NULL
        }
        new <- opt[[i]]
        class(new) <- class(value)
        
        opt <- opt[names(opt) != i]
        opt[[level1]][[level2]] <- new
        
      } else {
        value <- base@options[[name_algo]]@args.values$`_allData_allRun`[[i]]
        if (missing(value)){
          value <- NULL
        }
        new <- opt[[i]]
        class(new) <- class(value)
        opt[[i]] <- new
      }
    }
    
    if (is.null(user.val[[name_algo]])){
      user.val[[name_algo]] <- list()
    }
    for (dataset in datasets){
      user.val[[name_algo]][[dataset]] <- opt
    }
    
  }
  user.val
} 
  

### Progress functions
progress_plot <- function(progress_table){
  if (!is.null(progress_table)){
    p <- ggplot(progress_table, aes(x = RUN, y = PA, fill = progress))+
      geom_tile(color = "white",
                lwd = 1.5,
                linetype = 1) +
      coord_fixed()+
      scale_fill_manual(values = c("TRUE" = "#9BCD9B", "FALSE" = "#eaa704"))
  } else {
    p <- NULL
  }
  p
}
  

init_progress <- function(PA, CV){
  PA <- PA[PA != "All the same"]
  if (length(CV) == 1 && CV == "All the same"){
    CV <- "allRun"
  } else{
    CV <- CV[CV != "All the same"]
  }
  progress_table <- expand.grid(PA = PA, RUN = CV)
  progress_table$progress <- FALSE
  progress_table
}

update_progress <- function(progress_table, whichPA, whichCV){
  if (whichPA == "All the same"){
    if (whichCV == "All the same"){
      progress_table$progress <- TRUE
    } else {
      progress_table[progress_table$RUN == whichCV, "progress"] <- TRUE
    }
  } else {
    if (whichCV == "All the same"){
      progress_table[progress_table$PA == whichPA, "progress"] <- TRUE
    } else {
      progress_table[(progress_table$PA == whichPA & progress_table$RUN == whichCV), "progress"] <- TRUE
    }
  }
  progress_table
}