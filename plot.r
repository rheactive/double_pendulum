library(ggplot2)

flnm <- "results/double_pendulum.dat"

d <- read.table(flnm)

#Column height
len <- length(d$V1)

#Create a new data frame with "long form" like this:

#___x1___y1___"y"
#___x2___y2___"y"
#___x3___y3___"y"
#___x4___y4___"y"
#___x5___y5___"y"
#___x1___z1___"z"
#___x2___z2___"z"
#___x3___z3___"z"
#___x4___z4___"z"
#___x5___z5___"z"

#This is how ggplot2 works

d1 <- data.frame(x = c(d$V1, d$V1), y = c(d$V2, d$V3), group = c(rep("Angle 1", len), rep("Angle 2",len)))

#Create the plot
g <- ggplot(d1, aes(x = x, y = y, colour = group)) + 
#Draw lines
#geom_line(size=1.5) + 
#Draw points
geom_point(size = 0.7) +
#Plot style parameters
theme_bw(base_size = 28) + 
theme(aspect.ratio = 0.618, axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black", angle = 90, hjust = 0.5)) + 
#Legend parameters
theme(legend.title = element_blank(), legend.text = element_text(size = 28), legend.background = element_rect(colour = "black", fill = "white", size = 0.7)) + 
theme(legend.position = "right", legend.spacing.y = unit(10, "pt"), legend.margin=margin(t=0,l=20,b=15,r=20, unit='pt')) + 
guides(colour = guide_legend(byrow = TRUE)) +
#Axes labels
labs(x = "Time, seconds", y = "Angle, degrees")

#Show the plot on screen
g

#Export the plot

flnm <- paste("Angles_vs_time.png",sep="")

ggsave(flnm,scale=1,width=16,height=10)


d1 <- data.frame(x = d$V2, y = d$V3)

#Create the plot
g <- ggplot(d1, aes(x = x, y = y)) + 
#Draw lines
#geom_line(size=1.5) + 
#Draw points
geom_point(size = 0.7, colour = "purple") +
#Plot style parameters
theme_bw(base_size = 28) + 
theme(aspect.ratio = 1, axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black", angle = 90, hjust = 0.5)) + 
#Legend parameters
#theme(legend.title = element_blank(), legend.text = element_text(size = 28), legend.background = element_rect(colour = "black", fill = "white", size = 0.7)) + 
#theme(legend.position = "right", legend.spacing.y = unit(10, "pt"), legend.margin=margin(t=0,l=20,b=15,r=20, unit='pt')) + 
#guides(colour = guide_legend(byrow = TRUE)) +
#Axes labels
labs(x = "Angle 1, degrees", y = "Angle 2, degrees")

#Show the plot on screen
g

#Export the plot

flnm <- paste("Angles_phase_diagram.png",sep="")

ggsave(flnm,scale=1,width=13,height=13)


