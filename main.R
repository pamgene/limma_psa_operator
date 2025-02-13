library(tercen)
library(dplyr)
library(limma)
library(statmod)
library(reshape2)

ctx <- tercenCtx()

getData = function(ctx){
  if(length(ctx$labels) < 1) stop("Define samples using labels (e.g. Barcode, Row)")
  
  if(!ctx$hasXAxis) stop("Define grouping using the x-axis")
  
  if(length(ctx$colors)>1) stop("More than one pairing factor is not allowed")
  
  df = ctx %>% 
    select(.ri, .ci, .y, .x) %>%
    dplyr::mutate( grp  = .x  %>%  as.factor ) %>% 
    dplyr::mutate(obs = ctx$select(ctx$labels) %>%  interaction() %>%  droplevels())
  
  if(length(ctx$colors) > 0){
    df = df %>% 
      dplyr::mutate(.pf = ctx$select(ctx$colors) %>% pull(1))
  }
  return(df)
}

CONTR_REVERSE <- ctx$op.value('ReverseContrast', as.logical, TRUE)

trtContrasts = function(grp){
  lvgrp = levels(grp)
  CM = matrix(nrow = length(lvgrp), ncol = length(lvgrp))
  
  if(CONTR_REVERSE){
    CM = lower.tri(CM)
  } else {
    CM = upper.tri(CM)
  }
  dimnames(CM) = list(lvgrp, lvgrp)
  tc = CM %>% 
    reshape2::melt() %>% 
    filter(value) %>% 
    select(-value) %>% 
    mutate(contrast = paste(Var1, Var2, sep="-")) 
}

supergroups = function(ctx){
  ctx$cselect() %>% 
    dplyr::mutate(.ci = 0:(n()-1))
}

limmaFun = function(df){
  if (".pf" %in% colnames(df)){
    result = doLimmaWithPairing(df)
  } else {
    result = doLimmaStraight(df)
  }
  return(result)
}

doLimmaStraight = function(df){
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
      tt = topTable(fresult, coef = .$contrast, number = dim(X)[1], adjust.method = "BH")
      tt %>% 
        mutate(.ri = rownames(tt))
    }) %>% 
    ungroup()
}

doLimmaWithPairing = function(df){
  X = df %>% 
    acast(.ri ~obs, value.var = ".y")
  
  grp = df %>% 
    acast(.ri ~ obs, value.var = "grp")
  
  grp = grp[1,] %>% 
    make.names() %>% 
    as.factor() 
  
  pf = df %>% 
    acast(.ri ~obs, value.var = ".pf")
  
  pf =pf[1,] %>% 
    make.names() %>% 
    as.factor()
  
  mm = model.matrix(~ 0 + grp)
  colnames(mm) = levels(grp)
  contr.df = trtContrasts(grp)
  contr.mat = makeContrasts(contrasts = contr.df %>%  pull(contrast), levels = levels(grp))
  
  dc = duplicateCorrelation(X, design = mm, block = pf)
  
  fit = lmFit(X, mm, block = pf, correlation = dc$consensus.correlation) 
  
  fresult = fit %>% 
    contrasts.fit(contr.mat) %>% 
    eBayes()
  
  result = contr.df %>% 
    group_by(contrast) %>% 
    do({
      tt = topTable(fresult, coef = .$contrast, number = dim(X)[1], adjust.method = "BH")
      tt %>% 
        mutate(.ri = rownames(tt))
    }) %>% 
    ungroup()
}


result = ctx %>% 
  getData() %>% 
  group_by(.ci) %>% 
  do(limmaFun(.))%>% 
  ungroup() %>% 
  select(.ci, .ri, contrast, logFC, AveExpr, t, pvalue = P.Value, FDR = adj.P.Val) %>% 
  mutate(logp= -log10(pvalue)) %>% 
  group_by(.ci, contrast) %>% 
  mutate(rankp = rank(pvalue)) %>% 
  ungroup() %>% 
  mutate(.ri = as.integer(.ri),
         contrast = sub("-", " vs ", .$contrast))

if(max(result$.ci) >0) {
  result = result %>%  
    left_join(supergroups(ctx), by = ".ci")
}

result %>%
  ungroup() %>% 
  ctx$addNamespace() %>% 
  ctx$save()

