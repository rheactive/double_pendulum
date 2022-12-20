library(ggplot2)

flnm <- "dp_results/double_pendulum.dat"

dir.create("pictures")

d <- read.table(flnm)

t <- as.numeric(d$V1[-1]);
theta <- as.numeric(d$V2[-1]);
phi <- as.numeric(d$V3[-1]);
P <- as.numeric(d$V4[-1]);
Q <- as.numeric(d$V5[-1]);
cf <- as.numeric(d$V6[-1]);
H <- as.numeric(d$V7[-1]);

#Column height
len <- length(t)

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

d1 <- data.frame(x = c(t, t), y = c(theta, phi), group = c(rep("θ", len), rep("ϕ", len)))

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

flnm <- paste("pictures/Angles_vs_time.png",sep="")

ggsave(flnm,scale=1,width=16,height=10)


d1 <- data.frame(x = theta, y = phi)

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
labs(x = "theta, degrees", y = "phi, degrees")

#Show the plot on screen
g

#Export the plot

flnm <- paste("pictures/Angles_phase_diagram.png",sep="")

ggsave(flnm,scale=1,width=13,height=13)



d1 <- data.frame(x = c(t, t), y = c(H/max(H), cf), group = c(rep("H, a.u.", len), rep("c, a.u.", len)))

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
labs(x = "Time, seconds", y = "Total energy, normalized")

#Show the plot on screen
g

#Export the plot

flnm <- paste("pictures/Energy_and_friction.png",sep="")

ggsave(flnm,scale=1,width=16,height=10)


