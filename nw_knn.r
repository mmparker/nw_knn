


nw_knn <- function(x, true.class, k = 3, dist.method = "euclidean") {


    # Create distance matrix 
    Dist.matrix <- make.dist(x, dist.method)

    # Preallocate neighbor results vector
    n.pos.neighbors <- vector(mode = "numeric", length = nrow(Dist.matrix))

    # Preallocate the results vector
    KNN.class <- vector(mode = "character", length = nrow(Dist.matrix))
    
    # Could be an apply(Dist.matrix, INDEX = 1, FUN = )
    for (i in seq_len(nrow(Dist.matrix))) {   #  for each distance metric

        #  STEP 2 - identify the k profiles closest to sample i.  The smallest distance (zero)
        #  is always going to be the sample i.  This line identifies the the 2nd 
        #  through kth value in each distance vector
        nearest.indices <- as.integer(substring(labels(Dist.matrix[i, order(Dist.matrix[,i])[2:(k+1)]]), 2, 3))
    
        # so that we can use it to extract class info      
        nearest.classes <- true.class[nearest.indices]   # this gives the class labels for the k nearest neighbors, 
        # ie [1] "ER+" "ER+" "ER+"
        n.pos.neighbors[i] <- sum(nearest.classes %in% "ER+")     #  determines many of the nearest neighbors are ER+

        #  T[S/k == .5] <- "TIE" 
        KNN.class[i] <- predicted.class  # A vector of "predicted ER+/ER- class labels

    }

    KNN.class[n.pos.neighbors / k > 0.5] <- "ER+"
    KNN.class[n.pos.neighbors / k < 0.5] <- "ER-"


    # Return the classification
    KNN.class

}
