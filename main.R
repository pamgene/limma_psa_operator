library(tercen)
library(dplyr)
library(limma)
library(reshape2)

ctx <- tercenCtx()

getData = function(ctx){
  if(length(ctx$labels) < 1) stop("Define samples for each UKA using labels in Tercen")
  
  if(!ctx$hasXAxis) stop("Define grouping using the x-axis in Tercen")
  
  ctx %>% 
    select(.ri, .ci, .y, .x) %>%
    dplyr::mutate( grp  = .x  %>%  as.factor ) %>% 
    dplyr::mutate(obs = ctx$select(ctx$labels) %>%  interaction() %>%  droplevels())
  
}

trtContrasts = function(grp){
  lvgrp = levels(grp)
  CM = matrix(nrow = length(lvgrp), ncol = length(lvgrp)) %>% 
    upper.tri()
  dimnames(CM) = list(lvgrp, lvgrp)
  tc = CM %>% 
    reshape2::melt() %>% 
    filter(value) %>% 
    select(-value) %>% 
    mutate(contrast = paste(Var1, Var2, sep="-")) 
}

limmaFun = function(df){
  X = df %>% 
    acast(.ri ~obs, value.var = ".y")
  
  grp = df %>% 
    acast(.ri ~ obs, value.var = "grp")
  
  grp = grp[1,] %>% 
    make.names() %>% 
    as.factor() 
  
  mm = model.matrix(~ 0 +grp)
  colnames(mm) = levels(grp)
  
  fit = lmFit(X, mm)
  contr.df = trtContrasts(grp)
  contr.mat = makeContrasts(contrasts = contr.df %>%  pull(contrast), levels = levels(grp))
  
  fresult = lmFit(X, mm) %>% 
    contrasts.fit(contr.mat) %>% 
    eBayes()
  
  result = contr.df %>% 
    group_by(contrast) %>% 
    do({
      tt = topTable(fresult, coef = .$contrast, number = dim(X)[1])
      tt %>% 
        mutate(.ri = rownames(tt))
    }) %>% 
    ungroup()
}

result = ctx %>% 
  getData() %>% 
  limmaFun() %>% 
  select(.ri, contrast, logFC, AveExpr, t, P.Value) %>% 
  mutate(.ri = as.integer(.ri))

result %>%   
  ctx$addNamespace() %>% 
  ctx$save()
