library(ecr)
library(tidyverse)
library(scales)  # Ensure the scales package is loaded for the alpha() function

#Load in source files
setwd('C:/Users/Admin/Desktop/Github/JRSS_Lithium/Publish/src')
file_list = list.files()
sapply(file_list, source)

setwd('C:/Users/Admin/Desktop/Github/JRSS_Lithium/Publish/data')
multi_results = read.csv('mo_results.csv')

multi_results[,'Fitness1'] = multi_results[,'Fitness1'] / max(multi_results[,'Fitness1']) 
multi_results[,'Fitness2'] = multi_results[,'Fitness2'] / max(multi_results[,'Fitness2']) 

pareto_front = ecr::which.nondominated(t(multi_results[,c('Fitness1','Fitness2')]) * -1)

multi_results$pareto = 'Dominated'
multi_results[pareto_front,'pareto'] = 'Nondominated'
multi_results$par_size = ifelse(multi_results$pareto, 1.5, 1)
multi_results = data.frame(multi_results)# %>%
  # mutate(Fitness1 = Fitness1 * 100,
         # Fitness2 = Fitness2 * 100)

multi_plot = ggplot(multi_results, aes(x = Fitness1, y = Fitness2)) + 
  geom_point(data = subset(multi_results, pareto == "Dominated"), 
             aes(color = pareto, shape = pareto), size = 2, alpha = 0.6) +  # Smaller size for Dominated points
  geom_point(data = subset(multi_results, pareto == "Nondominated"), 
             aes(color = pareto, shape = pareto), size = 3, alpha = 0.9) +
  # geom_point(aes(color = pareto, shape = pareto)) +  # Pareto points
  geom_point(aes(x = 1, y = 1, color = 'Ideal Point', shape = 'Ideal Point'), size = 3) +  # Ideal Point
  labs(x = expression(D[s3]~"-"~efficiency),
       y = expression(D[s8]~"-"~efficiency),
       color = 'Pareto Optimal',  # Single title for the legend
       shape = 'Pareto Optimal',
       size = NULL) +  # Single title for the shape
  theme_bw() +
  scale_color_manual(
    values = c('Ideal Point' = 'black', 'Dominated' = 'grey', 'Nondominated' = 'black'),  # Reorder colors
    breaks = c('Ideal Point', 'Dominated', 'Nondominated')  # Order in the legend
  ) +
  scale_shape_manual(
    values = c('Ideal Point' = 8, 'Dominated' = 16, 'Nondominated' = 16),  # Reorder shapes
    breaks = c('Ideal Point', 'Dominated', 'Nondominated')  # Order in the legend
  ) +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +  # Show x-axis labels as percentages
  scale_y_continuous(labels = scales::percent)
multi_plot

setwd('C:/Users/Admin/Desktop/Github/JRSS_Lithium/Publish/')
ggsave('figure_3.png',plot = multi_plot, dpi = 1200,
       width = 10, height = 5)







