library(Cairo) 

# if not installed, install Cairo by typing install.library("Cairo"), then run this example again

# Example plot
plot(1:10, main = "Cairo Plot Example")

# Saving as PNG with Cairo
png("cairo_plot.png", type = "cairo")
plot(1:10, main = "Cairo PNG")
dev.off()

# Saving as PDF with Cairo (with ggsave)
library(ggplot2)
p <- ggplot(data.frame(x = 1:10, y = 1:10), aes(x, y)) + geom_point()

p # does it look better on screen?

ggsave("cairo_plot.pdf", plot = p, device = cairo_pdf)



