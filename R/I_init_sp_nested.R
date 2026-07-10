# Function to initialise (log)-smoothing parameters for GAMs with nested effects.
# Code taken from mgcv::estimate.gam from mgcv version 1.9.4 (9th July 2026)
# Should initialised smoothing parameters as mgcv for standard effects, but for
# nested effects .my.initial.spg does something different.
.init.sp.nested <- function(G, start) {
  
  if (!is.null(G$family$preinitialize)) {
    stop("gam_nl does not know how to use family$preinitialise")
  }
  
  if( !is.null(G$L) ){
    stop("gam_nl does not know how to handle G$L")
  }

  if (length(G$sp) == 0) return(NULL)

  is_nested <- sapply(G$smooth, function(sm) inherits(sm, "nested"))
  if (!any(is_nested)) return(NULL)

  nested_range <- lapply(G$smooth[is_nested], function(sm) c(sm$first.para, sm$last.para))

  block_is_nested <- sapply(G$off, function(o) {
    any(vapply(nested_range, function(r) o >= r[1] && o <= r[2], logical(1)))
  })

  if (!any(block_is_nested)) return(NULL)
  
  #
  ### START: Code of the first part of mgcv::estimate.gam as of mgcv version 1.9.4, 9th of July 2026.
  # We are interested in the code up to 
  #
  # lsp2 <- if (length(G$sp)>0) log(initial.spg(G$X,G$y,G$w,G$family,G$S,G$rank,G$off,
  # offset=G$offset,L=G$L,lsp0=G$lsp0,E=G$Eb,...))  else rep(0,0)
  #
  # (excluded below). We want to initialise the smoothing parametes in a different way, but we need to prepare the G
  # object in the same way as mgcv does, otherwise our smoothing parameters will not be compatible with what 
  # estimate.gam does.
  
  # Same reparametrisation estimate.gam() applies to G$X before calling initial.spg() for
  # general.family models - see the header comment above.
  G$Sl <- Sl.setup(G) ## prepare penalty sequence
  
  G$X <- Sl.initial.repara(G$Sl,G$X,both.sides=FALSE) ## re-parameterize accordingly
  
  if (!is.null(start)){ 
    start <- Sl.initial.repara(G$Sl,start,inverse=FALSE,both.sides=FALSE)
  }
  
  G$Eb <- totalPenaltySpace(G$S,G$H,G$off,ncol(G$X))$E  ## balanced penalty square root for rank determination purposes 
  
  ### END
  
  lambda <- .my.initial.spg(G$X, G$y, G$w, G$family, G$S, G$rank, G$off, nested_idx = block_is_nested,
                           offset = G$offset, L = G$L, lsp0 = G$lsp0, E = G$Eb, start = start)
  
  return( lambda )
}
