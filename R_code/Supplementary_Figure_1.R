library(tidyverse)

#HCD
#Load in 
#Source files
# setwd("~/Optimal_Lithium/src")
setwd('~/JRSS_Lithium-main/Publish/src')
file_list = list.files()
sapply(file_list, source)


# Design parameters within the 95\% confidence bound limits
# Fixed effects
f_ka   = c(0.93 - 1.96*0.35449651, 0.93 + 1.96*0.35449651)
f_V    = c(22.3 - 1.96*7.59238787, 22.3 + 1.96*7.59238787)
f_CLS  = c(1.24 - 1.96*0.05837381, 1.24 + 1.96*0.05837381)
f_CLES = c(4.15 - 1.96*3.07473842, 4.15 + 1.96*3.07473842)
f_CLSE = c(11.1 - 1.96*8.49259113, 11.1 + 1.96*8.49259113)
# Random effects
# r_ka   = seq(0.3, 1.5, length.out = 2)
# r_V    = seq(0.1, 0.5, length.out = 2)
# r_ClS  = seq(0.1, 0.4, length.out = 2)
# r_CLES = seq(0.1, 0.5, length.out = 2)
# r_CLSE = 0
# f_sigma = 0
# r_sigma  = seq(0.1, 0.3, length.out = 2)names_test = as.character(c('f_ka', 'f_V', 'f_CLS',' f_CLES', 'f_CLSE'))
grid_search = expand.grid(f_ka, f_V, f_CLS, f_CLES, f_CLSE)
dim(grid_search)
names(grid_search) = names_test
grid_df = data.frame(grid_search)
# 
# ans = c()
# ans2 = matrix(nrow = 32, ncol = 5)
# 
# ans = data.frame()
# # ln_sum = 0
# # HCD <- function(time){
# for(jj in 1:32){
#   # print(jj)
#   beta  = as.numeric(grid_df[jj, 1:5])
#   o     = c(0.72,0.30,0.20,0.27,0)#as.numeric(grid_search[jj, 6:10])
#   sigma = c(0, 0.137)#as.numeric(grid_search[jj, 11:12])
#   Number_of_subjects = 100
#   #Function
#   Dopt_lith_1grp_SS3<-function(time){
#     res<-funFIMem(formSS3,c("KA","V","ClS", "CLES", "CLSE"),
#                   beta,o,sigma,
#                   list(time[1:5]),2,c(36),c(time[6]/sum(time[6])),Number_of_subjects)
#     return(-res)
#   }
#   # ln_sum = ln_sum + log(-Dopt_lith_1grp_SS3(time))
#   # ans = 
#   number_of_parameters = 6
#   parameter_bounds = cbind(rep(0, number_of_parameters),rep(8, number_of_parameters))
#   library(tictoc)
#   tic()
#   result6_cov_1grp<-optim_pso(objective_function = Dopt_lith_1grp_SS3, number_of_parameters = number_of_parameters, number_of_particles = 40, 
#                               max_number_of_iterations = 20,max_number_function_calls= 10000,
#                               abstol = 0.0001, reltol = 0.0001,
#                               parameter_bounds = parameter_bounds, 
#                               initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
#                               lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
#   toc()
#   vals = result6_cov_1grp$par %>%
#     data.frame()
#   if(jj == 1){
#     ans = vals
#   }
#   ans = cbind(ans, vals)
#   print(jj)
# }
# 
# work = ans[,-1]
# 
# work[6,] = 1:ncol(work)  
# work[1:5,] = apply(work[1:5,], 2, sort)
# work2 = work[-6,] %>%
#   data.frame()
# names(work2) = work[6,]
# work3 = work2 %>% 
#   pivot_longer(cols = everything(),
#                names_to = 'test')
# work3$test = rep(c(1:5), each = 32)
# setwd('..')
# setwd('Publish/data/')
# write.csv(work3, 'HCD_results.csv', row.names = F )

# setwd('data/')
setwd('~/JRSS_Lithium-main/Publish/data')
hcd_results = read.csv('HCD_results.csv')

hcd_plot = ggplot(data = hcd_results) +
  geom_count(aes(x = test, y = value)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24)) +
  labs(x = "Design point",
       y = "Time in hours") +
  guides(size = F)
hcd_plot


setwd('~/JRSS_Lithium-main/Publish')
ggsave('Web_appendix_figure_1.png', plot = hcd_plot, dpi = 1200, width = 8, height = 8)

