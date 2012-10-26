



make.dist <- function(DATA, METHOD) {
  as.matrix(dist(t(DATA),   #  to get the distance matrix
                 #  we have to transpose the expression matrix
                 #  to get genes in columns and then apply 
                 #  the distance function with selected method
                 method = METHOD))    #  46 x 46 matrix
}


#  the METHOD could also be set as"maximum", "manhattan", "canberra", "binary" 
#or "minkowski"
