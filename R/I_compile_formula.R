
.compile_formula <- function(form){
  if( is.list(form) ){ # Recursion on list elements
    form <- lapply(form, .compile_formula)
  } else{
    # Decompose formula into terms
    dec <- terms(form, keep.order = TRUE)
    ter <- attr(dec, "term.labels")
    
    # Is there an offset()? If so add it to terms
    fact <- rownames( attr(dec, "factors") )
    faci <- which(grepl("offset(", fact, fixed = TRUE))
    if( length(faci) ){
      ter <- c(fact[ faci ], ter)
    }
    
    # Has the intercept been removed? If so add "-1" to the terms
    if(attr(dec, "intercept") == 0){ ter <- c("-1", ter)  }
    
    # Look for nested effects and replace "s_nest(...)" with standard "s(...)" 
    ine <- which(startsWith(ter, "s_nest"))
    if(length(ine)){
      for(ii in ine){
        tmp <- ter[ii]
        ter[ii] <- eval(parse(text=tmp))$form_term
      }
      # Rebuild model formula with modified terms
      form <- as.formula(paste0(all.vars(form)[1], "~ ", paste0(ter, collapse = " + ")))
    }
  }
  return(form)
}