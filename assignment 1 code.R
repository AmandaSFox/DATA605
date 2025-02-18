library(tidyverse)
library(ggplot2)
library(gganimate)
# library(gifski)

# --------------
# plot a triangle
# --------------

# create matrix defining a triangle:
my_shape <- matrix(c(0,0,
                     1,0,
                     0,1,
                     0,0), 
                   nrow = 2, byrow = TRUE)
my_shape

# prepare to plot: convert to df and label columns x & y
df_my_shape <- as.data.frame(t(my_shape))  
colnames(df_my_shape) <- c("x","y")

df_my_shape

# plot
plot_my_shape <- df_my_shape %>% 
  ggplot(aes(x = x, y = y)) +
         geom_path(color = 'blue', linewidth = 1.2) +
         geom_hline(yintercept = 0, color = "gray30") +
         geom_vline(xintercept = 0, color = "gray30") +
         coord_fixed(xlim = c(-2,2), ylim = c(-2,2)) +
         theme_minimal()

plot_my_shape + plot_annotation("Original Triangle")

# --------------
# scale by 2
# --------------

# create scaling transformation matrix multiplying x and y by 2 
matrix_scale <- matrix(c(2,0,
                         0,2),
                       nrow = 2, byrow = TRUE)
matrix_scale

# multiply: transpose my_shape matrix so # rows match # columns 
my_shape_scaled <- t(my_shape) %*% matrix_scale 
my_shape_scaled

# convert to df with column names x,y
df_my_shape_scaled <- as.data.frame(my_shape_scaled)  
colnames(df_my_shape_scaled) <- c("x","y")

df_my_shape_scaled

# plot
plot_my_shape_scaled <- df_my_shape_scaled %>% 
  ggplot(aes(x = x, y = y)) +
  geom_path(color = 'blue', linewidth = 1.2) +
  geom_hline(yintercept = 0, color = "gray30") +
  geom_vline(xintercept = 0, color = "gray30") +
  coord_fixed(xlim = c(-2,2), ylim = c(-2,2)) +
  theme_minimal()

plot_my_shape + plot_my_shape_scaled  + plot_annotation("Original and Scaled Triangle")

# --------------
# flip over x axis
# --------------

# create transformation matrix (x, y) → (x, -y)
matrix_flip <- matrix(c(1, 0,
                        0,-1),
                       nrow = 2, byrow = TRUE)
matrix_flip

# multiply: transpose my_shape matrix so # rows match # columns 
my_shape_flip <- t(my_shape) %*% matrix_flip 
my_shape_flip

# convert to df with column names x,y
df_my_shape_flip <- as.data.frame(my_shape_flip)  
colnames(df_my_shape_flip) <- c("x","y")

df_my_shape_flip

# plot
plot_my_shape_flip <- df_my_shape_flip %>% 
  ggplot(aes(x = x, y = y)) +
  geom_path(color = 'blue', linewidth = 1.2) +
  geom_hline(yintercept = 0, color = "gray30") +
  geom_vline(xintercept = 0, color = "gray30") +
  coord_fixed(xlim = c(-2,2), ylim = c(-2,2)) +
  theme_minimal()

plot_my_shape + plot_my_shape_flip + plot_annotation("Original and Flipped Triangle Over x Axis")

# ------------
# rotate and store fixed plot images
# ------------

# create transformation matrix (x, y) → (-y, x)
matrix_rotate_90 <- matrix(c(0,-1,
                             1, 0),
                           nrow = 2, byrow = TRUE)
matrix_rotate_90

# create list to hold plots for display
rotated_shapes <- list()
rotated_shapes[[1]] <- plot_my_shape

current_shape <- t(my_shape)
current_shape

# loop: multiply, convert to df, create and store plot

for (i in 2:5) {
  
  # multiply matrix
  current_shape <- current_shape %*% matrix_rotate_90 
  current_shape
  
  # convert transformed matrix to df with column names x,y
  df_rotate_90 <- as.data.frame(current_shape)  
  colnames(df_rotate_90) <- c("x","y")
  
  df_rotate_90
  
  # create and store plot
  rotated_shapes[[i]] <- df_rotate_90 %>% 
    ggplot(aes(x = x, y = y)) +
    geom_path(color = 'blue', linewidth = 1.2) +
    geom_hline(yintercept = 0, color = "gray30") +
    geom_vline(xintercept = 0, color = "gray30") +
    coord_fixed(xlim = c(-2,2), ylim = c(-2,2)) +
    theme_minimal() 
  }

plot_display <- rotated_shapes[[1]] +
                rotated_shapes[[2]]+
                rotated_shapes[[3]]+
                rotated_shapes[[4]]+
                rotated_shapes[[5]]

plot_display + 
  plot_layout(ncol = 5) +
  plot_annotation("Original and Rotated Triangle: Clockwise 360 Degrees")



















# --------------
# rotation: counterclockwise rotation animation sequence
# using loop and gganimation package
# --------------

# create transformation matrix (x, y) → (-y, x)
matrix_rotate_90 <- matrix(c(0,-1,
                             1, 0),
                        nrow = 2, byrow = TRUE)
matrix_rotate_90

# create empty df
df_animation <- data.frame()

# current shape
current_shape <- my_shape

# loop: append current shape matrix to df and rotate
for (i in 0:3) {
  temp <- as.data.frame(t(current_shape))
  colnames(temp) <- c("x","y")
  temp$frame <- i
  
  df_animation <- rbind(df_animation,temp)
  
  # check i and rotate
  if(i<3)  {current_shape <- matrix_rotate_90 %*% current_shape}
}

#df_animation

plot_animated <- df_animation %>% 
  ggplot(aes(x=x,y=y, group = frame))+
  geom_path(color = "blue", linewidth = 1.5) +
  coord_fixed(xlim = c(-2,2), ylim = c(-2,2)) +
  theme_minimal() +
  transition_states(frame, transition_length = 2, state_length = 1, wrap = TRUE) 

animate(plot_animated, 
        nframes = 5, 
        fps = 1, 
        renderer = file_renderer(dir = "C:/Users/amand/Documents/", 
                                 prefix = "triangle_rotation_", 
                                 overwrite = TRUE))


















# --------------
# NOT USED - rotate four times -90 degrees
# --------------

# create list of matrices with original in position 1
rotated_shapes <- list()
rotated_shapes[[1]] <- my_shape

# loop and multiply to create rotated matrices
for (i in 2:5) {rotated_shapes[[i]] <- matrix_rotate_90 %*% rotated_shapes[[i-1]]}

# convert new matrices to dfs for plotting
df_my_shape_rotate_90 <- as.data.frame(t(rotated_shapes[[2]]))
df_my_shape_rotate_180 <- as.data.frame(t(rotated_shapes[[3]]))
df_my_shape_rotate_270 <- as.data.frame(t(rotated_shapes[[4]]))
df_my_shape_rotate_360 <- as.data.frame(t(rotated_shapes[[5]]))
