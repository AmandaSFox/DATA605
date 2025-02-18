"In the Vector unit, section SSNS discusses the Spanning Sets of Null Spaces. 
This chapter helped me understand more geometrically what the set of solutions to a matrix Ax = 0
looks like by not showing the variables as the scalars (because they can be anything when creating the set, they are implied).
"

library(tidyverse)
library(MASS)
library(pracma)


"First I created the matrix by using the c() function to create a vector and then use nrow & byrow 
to wrap the vector by rows into a matrix with the correct dimensions. R can also wrap column-wise (the default)."

# create matrix A
A <- matrix(c(1,3,1,9,2,1,-3,8,1,1,-1,5), 
            nrow = 3,
            byrow = TRUE)
A

"The MASS package uses the Null() function to find the null space of A in a single step. 
The output is the set of basis vectors."

A_null_space_vectors <- Null(A) 
A_null_space_vectors

"
Here it went off the rails. You can see that this single vector is not the correct solution 
as there are two free variables in the row reduced echelon format below, using the pracma package:"
A_rref <- rref(A)
A_rref

"Doing this manually on paper, I got this set of two basis vectors, which I verified with an online calculator: "
solution <- matrix(c(2,-1,1,0,-3,-2,0,1),ncol=2)
solution

"The pracma package did a little better:"
null_space <- nullspace(A)  # Robustly computes the null space
print(null_space)

"This gave me two vectors but still as decimals. I believe this has to do with scaling
but I was unable to figure out how to get it to integers.

For a quick fix, my friend ChatGPT 
provided some complicated code to scale it to integers, which 
unsurprisingly gave me a glaringly wrong result, so I tried a third approach."

integer_null_space <- apply(null_space, 2, function(col) {
  scale_factor <- 1 / min(abs(col[col != 0]))  # Scale to smallest non-zero element
  round(col * scale_factor)  # Convert to integers
})
integer_null_space

"Finally I found a solution on stackoverflow.com. Some kind soul had created the below
formula which worked perfectly. I'm still picking it apart to really understand it, but
I wanted to share here for thoughts."

# Function from kind soul
# "https://stackoverflow.com/questions/43223579/solve-homogenous-system-ax-0-for-any-m-n-matrix-a-in-r-find-null-space-basi)"

NullSpace <- function (A) {
  m <- dim(A)[1]; n <- dim(A)[2]
  ## QR factorization and rank detection
  QR <- base::qr.default(A)
  r <- QR$rank
  ## cases 2 to 4
  if ((r < min(m, n)) || (m < n)) {
    R <- QR$qr[1:r, , drop = FALSE]
    P <- QR$pivot
    F <- R[, (r + 1):n, drop = FALSE]
    I <- base::diag(1, n - r)
    B <- -1.0 * base::backsolve(R, F, r)
    Y <- base::rbind(B, I)
    X <- Y[base::order(P), , drop = FALSE]
    return(X)
  }
  ## case 1
  return(base::matrix(0, n, 1))
}

#Call function
X <- NullSpace(A)
round(X, 15)

X







