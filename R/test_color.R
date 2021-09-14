
hcl_palettes("Sequential (single-hue)", n = 7, plot = TRUE)

viridis_hcl = colorspace::sequential_hcl(num_color_slices,
                          
                           power = c(0.8, 1.2)
                           
                           , palette = 'Blues 2'
                           )



library(scales)

show_col(viridis_hcl)


p <- ggplot(mtcars, aes(wt, mpg))
p + geom_point(size=4, aes(colour = factor(cyl))) +
  scale_color_viridis(discrete=T) +
  theme_bw()

p <- ggplot(mtcars, aes(wt, mpg))
p + geom_point(size=4, aes(colour =disp)) +
  scale_color_viridis(discrete=F) +
  
  theme_bw()


library(ggplot2)
highcol = 'red4'
lowcol = 'pink'


p <- ggplot(mtcars, aes(wt, mpg))
p + geom_point(size=4, aes(colour =disp)) +
  scale_colour_gradient(high = highcol, low = lowcol)  


mtcars$disp[1] = NA



p <- ggplot(mtcars, aes(wt, mpg))
p + geom_point(size=4, aes(colour =disp)) +
  scale_colour_gradient(high = highcol, low = lowcol) + labs(colour="Displac") + theme_bw()
  

 

wc_example <- ggplot(worldcup, aes(x = Time, y = Passes,
                                   color = Position, size = Shots)) + 
  geom_point(alpha = 0.5) 
