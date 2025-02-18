library(tidyverse)

#define column vectors as matrices 
a <- matrix(c(2,4,6), ncol = 1)
b <- matrix(c(1,3,5), ncol = 1)

a
b

#scale b by 3
scaled_b <- 3* b
scaled_b

#add a + scaled_b
sum_a_scaled_b <- a + scaled_b
sum_a_scaled_b

# dot product of a and b
product_a_b <- t(a) %*% b
product_a_b

t(a) #transpose a and make it horizontal
a*b #multiply vectors

# dot products two ways
sum(a*b) # multiply vectors and sum the resulting vector, same as this:
t(a) %*% b #transpose one vector and multiply, use matrix math notation??

#Matrix multiplication requires the number of columns in the first matrix to equal the number of rows in the second matrix. Here, both 

