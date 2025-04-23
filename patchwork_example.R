library(data.table)
library(ggplot2)
library(tidyr)
library(patchwork) # better than gridExtra to align plots in a layout

fake <- data.table(ploidy=rnorm(1000,mean = c(0,4)),log_mut=rnorm(1000,mean = c(0,4)) %>% sample)

# simple example


g1 <- ggplot(fake) + aes(ploidy) +
  geom_histogram(bins=50,fill="lightblue") + 
  scale_x_continuous(position="top") 

g2 <- ggplot(fake) + aes(y=log_mut) + 
  geom_histogram(bins=50,fill="pink") + 
  scale_y_continuous(position="right") 

g3 <- ggplot(fake) + geom_point(aes(ploidy,log_mut))  

# plot g1 g2 g3 using patchwork
g1 + plot_spacer() + g3 + g2 + plot_layout(nrow=2,widths = c(3,1),heights=c(1,3))


# a bit nicer

bin1 <- .25 # adjust the width of bins

g1 <- ggplot(fake) + aes(ploidy) +
  geom_histogram(binwidth=binwidth,fill="lightblue") + 
  geom_density(aes(y=after_stat(density * bin1 * nrow(g1$data))),color="blue") +
  scale_x_continuous(position="top") +
  ylab("Count")

bin2 <- .25 # might be different than bin1

g2 <- g2 <- ggplot(fake) + aes(y=log_mut) + 
  geom_histogram(binwidth=bin2,fill="pink") + 
  geom_density(aes(x=after_stat(density * bin2 * nrow(g2$data))),color="red") +
  scale_y_continuous(position="right") + 
  xlab("Count")

# plot g1 g2 g3 using patchwork
g1 + plot_spacer() + g3 + g2 + plot_layout(nrow=2,widths = c(3,1),heights=c(1,3))



# visualize and manually determine the cutoffs (it can be done automatically too, but it's complex)

g1 <- g1 + geom_vline(xintercept = 2,lty=2)
g2 <- g2 + geom_hline(yintercept = 2,lty=2)
g3 <- g3 + geom_vline(xintercept = 2,lty=2) + geom_hline(yintercept = 2,lty=2)

g1 + plot_spacer() + g3 + g2 + plot_layout(nrow=2,widths = c(3,1),heights=c(1,3))






