cplot <- function(matrix) {
  
  loadings <- F
  
  if (class(matrix)[1] == "loadings") {
    loadings <- T
    matrix <- unclass(matrix)
  }
  
  out <- matrix %>% 
    as.data.frame() %>% 
    magrittr::set_colnames(colnames(matrix)) %>% 
    tibble::rownames_to_column("var2") %>% 
    gather("var1", "value", -var2) %>% 
    mutate(lbl = as.character(format(round(value, 2), nsmall = 2)),
           lbl = case_when(lbl == " 1.00" ~ "1",
                           value < 0 ~ paste("-", substr(lbl, 3, 5), sep = ""),
                           value > 0 ~ substr(lbl, 3, 5))) %>% 
    ggplot(aes(x = var1, y = var2, fill = value)) + 
    geom_tile() +
    geom_text(aes(label = lbl)) +
    scale_fill_gradient2(low = scales::muted("blue"), mid = "white", high = scales::muted("red"), 
                         limits = c(-1.001, 1.001)) +
    scale_y_discrete(limits = rev) +
    scale_x_discrete(position = "top") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  if (loadings) {
    out + 
      labs(fill = "Ladung") +
      theme(axis.title = element_blank())
  } else {
    out +
      labs(fill = "Korrelation") +
      theme(axis.title = element_blank())
  }
  
}

modelcor <- function(loadings, residuals) {
  
  m.cor <- matrix(rep(0, nrow(residuals)^2), nrow = nrow(residuals), 
                  dimnames = list(rownames(loadings), rownames(loadings)))
  
  for (i in 1:nrow(m.cor)) {
    
    for (j in 1:nrow(m.cor)) {
      
      m.cor[i,j] <- residuals[i,j]
      
      for (k in 1:ncol(loadings)) {
        m.cor[i,j] <- m.cor[i,j] + loadings[i,k]*loadings[j,k]
      }
      
    }
    
  }
  
  return(m.cor)
  
}

complexity <- function(loadings, kappa) {
  
  # Zeilenkomplexität
  zeilenkomp <- 0
  
  for (i in 1:nrow(loadings)) {
    
    for (j in 1:ncol(loadings)) {
      
      for (k in 1:ncol(loadings)) {
        
        if (j == k) {
          
          next
          
        } else {
          
          zeilenkomp <- zeilenkomp + loadings[i,j]^2 * loadings[i,k]^2
          
        }
        
      }
      
    }
    
  }
  
  
  # Spaltenkomplexität
  spaltenkomp <- 0
  for (j in 1:ncol(loadings)) {
    
    for (i in 1:nrow(loadings)) {
      
      for (k in 1:nrow(loadings)) {
        
        if (i == k) {
          
          next
          
        } else {
          
          spaltenkomp <- spaltenkomp + loadings[i,j]^2 * loadings[k,j]^2
          
        }
        
      }
      
    }
    
  }
  
  komp <- (1-kappa) * zeilenkomp + k * spaltenkomp 
  return(komp)
  
}




run.efa.unrotated <- function(method, data, m.cor) {
  
  out <- list()
  
  if (method == "pfa") {
    
    model1 <- fa(r = m.cor, nfactors = 3, n.obs = nrow(data), rotate = "none")
    
  } else if (method == "pca") {
    
    model1 <- pca(r = m.cor, nfactors = 3, n.obs = nrow(data), rotate = "none")
    
  }
  
  modelcor1 <- modelcor(loadings = model1$loadings, residuals = model1$residual)
  
  out$load.unrotated <- cplot(model1$loadings)
  out$cor.unrotated <- cplot(modelcor1)
  
  return(out)
  
}


run.efa.rotated <- function(method, rotation, data, m.cor) {
  
  out <- list()
  
  if (method == "pfa") {
    
    model2 <- fa(r = m.cor, nfactors = 3, n.obs = nrow(data), rotate = rotation)
    
  } else if (method == "pca") {
    
    model2 <- pca(r = m.cor, nfactors = 3, n.obs = nrow(data), rotate = rotation)
    
  }
  
  modelcor2 <- modelcor(loadings = model2$loadings, residuals = model2$residual)
  
  out$load.rotated <- cplot(model2$loadings)
  out$cor.rotated <- cplot(modelcor2)
  
  return(out)
  
}


run.simulations <- function(method, rotation, desired.comp, data, m.cor) {
  
  out <- list()
  
  nrep <- 200
  
  if (rotation == "varimax") {
    kappa <- 1/nrow(m.cor)
  } else if (rotation == "equamax") {
    kappa <- 0
  }
  
  
  if (method == "pfa") {
    
    model1 <- fa(r = m.cor, nfactors = 3, n.obs = nrow(data), rotate = "none")
    model2 <- fa(r = m.cor, nfactors = 3, n.obs = nrow(data), rotate = rotation)
    
  } else if (method == "pca") {
    
    model1 <- pca(r = m.cor, nfactors = 3, n.obs = nrow(data), rotate = "none")
    model2 <- pca(r = m.cor, nfactors = 3, n.obs = nrow(data), rotate = rotation)
    
  }
  
  modelcor1 <- modelcor(loadings = model1$loadings, residuals = model1$residual)
  modelcor2 <- modelcor(loadings = model2$loadings, residuals = model2$residual)
  
  
  # orthogonale Matrix erstellen
  ni <- nrow(model1$loadings)
  nf <- ncol(model1$loadings)
  
  
  rotationmatrix <- list()
  loadings.rot <- list()
  model.cor <- list()
  complexities <- c()
  
  for (i in 1:nrep) {
    
    rotationmatrix[[i]] <- rstiefel::rustiefel(m = nf, R = nf)
    loadings.rot[[i]] <- model1$loadings %*% rotationmatrix[[i]]
    colnames(loadings.rot[[i]]) <- colnames(model1$loadings)
    model.cor[[i]] <- modelcor(loadings = loadings.rot[[i]], residuals = model1$residual)
    complexities[i] <- complexity(loadings = loadings.rot[[i]], kappa = kappa)
    
  }
  
  
  # Sort complexities in descending order and get indices
  reps <- 1:nrep
  idx <- reps[order(complexities, decreasing = T)]
  
  
  #  Reorder the rotation matrices, loading matrices, and model correlation matrices
  rotationmatrix <- rotationmatrix[idx]
  loadings.rot <- loadings.rot[idx]
  model.cor <- model.cor[idx]
  complexities <- complexities[idx]
  
  # Find closest complexity to desired complexity
  comp.pct <- cumsum(complexities) / sum(complexities) * 100
  comp.idx <- which.min(abs(comp.pct - desired.comp))
  
  # Output
  out$load.random <- cplot(loadings.rot[[comp.idx]])
  out$cor.random <- cplot(model.cor[[comp.idx]])
  
  return(out)
  
}
