# Figure 1 creation
library(ppso)
library(tidyverse)
library(ggpubr)

#Load in 
#Source files
setwd('C:/Users/Admin/Desktop/Github/JRSS_Lithium2/Publish2/src')
file_list = list.files()
sapply(file_list, source)

# # Long term run of PSO
# formSS3<-"dose * (KA/V * ((CLSE/V) - KA)/((((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - KA) * (((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - KA)) * exp(-KA * (t - 0))/(1 - exp(-KA * 24)) + KA/V * ((CLSE/V) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))/((KA - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2)) * (((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))) * exp(-((((ClS/V) + 
# (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * (t - 0))/(1 - exp(-((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * 24)) + KA/V * ((CLSE/V) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))/((KA - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2)) * (((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))) * exp(-((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * (t - 
# 0))/(1 - exp(-((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * 24)))" 
# ###one group 5 samp times
# Number_of_subjects<-100
# Dopt_lith_1grp_SS3<-function(time){
#   res<-funFIMem(formSS3,c("KA","V","ClS", "CLES", "CLSE"),
#                 c(0.93,22.30,1.24,11.1,4.15),c(0.72,0.30,0.20,0.27,0),c(0, 0.137),
#                 list(time[1:5]),2,c(36),c(time[6]/sum(time[6])),Number_of_subjects)
#   return(-res)
# }
# 
# ###evaluation of a design only (-D-criterion of the FIM)
# time_test<-c(0,2,4,6,8,Number_of_subjects)
# Dopt_lith_1grp_SS3(time = time_test)
# #####
# 
# number_of_parameters<-6 
# parameter_bounds = cbind(rep(0, number_of_parameters),rep(8, number_of_parameters))
# ptm <- proc.time()
# pso_long<-optim_pso(objective_function = Dopt_lith_1grp_SS3, number_of_parameters = number_of_parameters, number_of_particles = 40, 
#                             max_number_of_iterations = 100,max_number_function_calls= 1000, 
#                             parameter_bounds = parameter_bounds, 
#                             initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
#                             lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
# proc.time() - ptm
# #
# The log file created from the above code
# Is used to create Figure 1
# The long term convergence of each design point

setwd('C:/Users/Admin/Desktop/Github/JRSS_Lithium2/Publish2/data')
d1 = read.delim("ppso_long.log")
d1$iteration = rep(c(1:100), each = 40)

Number_of_subjects<-100

number_of_parameters<-6 
parameter_bounds = cbind(rep(0, number_of_parameters),rep(8, number_of_parameters))
d2 = d1 %>%
  mutate(part_no = rep(c(1:40), 100)) %>%
  group_by(iteration) %>%
  mutate(min = min(objective_function),best = ifelse(objective_function == min(objective_function), 1, 0))

# Place initial starting point in log of 0, 2, 4, 6, 8h.
init = matrix(rep(c(0,2,4,6,8,1), 40), byrow = T,nrow = 40)
init2 = init %>%
  data.frame() %>%
  mutate(objective_function = D1cov(c(0,2,4,6,8,1)))
names(init2)[1:6] = names(d2)[2:7]
init2$iteration = 0
init2$best = 1
d3 = plyr::rbind.fill(init2, d2)
p1 = ggplot(data = d3) +
  geom_point(aes(x = iteration, y = parameter_1), color = 'grey',alpha = 1) +
  geom_line(data = d3[d3$best == 1,], 
            aes(x = iteration, y = parameter_1), lwd = 3) +
  scale_y_continuous(limits = c(0, 8)) +
  labs(title = 'Time Point 4',x = 'Iteration', y = '') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 24),
        axis.title = element_text(size = 22),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = 'none',
        axis.text.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank()
  )

p2 = ggplot(data = d3) +
  geom_point(aes(x = iteration, y = parameter_2), color = 'grey',alpha = 1) +
  geom_line(data = d3[d3$best == 1,], 
            aes(x = iteration, y = parameter_2), lwd = 3) +
  scale_y_continuous(limits = c(0, 8)) +
  labs(title = 'Time Point 3',x = 'Iteration', y = '') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 24),
        axis.title = element_text(size = 22),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = 'none',
        axis.text.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank()
  )

p3 = ggplot(data = d3) +
  geom_point(aes(x = iteration, y = parameter_3), color = 'grey',alpha = 1) +
  geom_line(data = d3[d3$best == 1,], 
            aes(x = iteration, y = parameter_3), lwd = 3) +
  scale_y_continuous(limits = c(0, 8)) +
  labs(title = 'Time Point 2', x = 'Iteration', y = '') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 24),
        axis.title = element_text(size = 22),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = 'none',
        axis.text.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank()
  )

p4 = ggplot(data = d3) +
  geom_point(aes(x = iteration, y = parameter_4), color = 'grey',alpha = 1) +
  geom_line(data = d3[d3$best == 1,], 
            aes(x = iteration, y = parameter_4), lwd = 3) +
  scale_y_continuous(limits = c(0, 8)) +
  labs(title = 'Time Point 1', x = 'Iteration', y = 'Time of Measurement (Hour)') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 24),
        axis.title = element_text(size = 22),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = 'none',
        axis.text.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank()
  )

p5 = ggplot(data = d3) +
  geom_point(aes(x = iteration, y = parameter_5), color = 'grey',alpha = 1) +
  geom_line(data = d3[d3$best == 1,], 
            aes(x = iteration, y = parameter_5), lwd = 3) +
  scale_y_continuous(limits = c(0, 8)) +
  labs(title = 'Time Point 5', x = 'Iteration', y = '') +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 24),
        axis.title = element_text(size = 22),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = 'none',
        axis.text.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank()
  )

t1 = ggarrange(p4, p3, p2, p1, p5, ncol = 5, nrow = 1)
t1
setwd('C:/Users/Admin/Desktop/Github/JRSS_Lithium2/Publish2/')

ggsave('figure_2.png', plot = t1,
       dpi = 1200, width = 18, height = 5)

