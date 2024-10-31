#Load in 
#Source files
# setwd("~/Optimal_Lithium/src")
setwd('~/JRSS_Lithium-main/Publish/src')
file_list = list.files()
sapply(file_list, source)

# setwd("~/Optimal_Lithium/data/PFIM_files")
setwd('~/JRSS_Lithium-main/Publish/data/PFIM_files')

## Table 1 analyses
start1 = proc.time()
t3 = PFIM('stdin_cov_3_points.r')
end1 = proc.time() - start1
end1

start2 = proc.time()
t4 = PFIM('stdin_cov_4_points.r')
end2 = proc.time() - start2
end2

start3 = proc.time()
t5 = PFIM('stdin_cov.r')
end3 = proc.time() - start3
end3

start4 = proc.time()
t6 = PFIM('stdin_cov_6_points.r')
end4 = proc.time() - start4
end4

start5 = proc.time()
t7 = PFIM('stdin_cov_7_points.r')
end5 = proc.time() - start5
end5


# Web Appendix B
## Table 2 analyses
start6 = proc.time()
pfim_dc = PFIM('stdin_cov.r')
end6 = proc.time() - start6
end6

start7 = proc.time()
pfim_dc_2 = PFIM('stdin_cov_2.r')
end7 = proc.time() - start7
end7

start8 = proc.time()
pfim_d = PFIM('stdin.r')
end8 = proc.time() - start8
end8

start9 = proc.time()
pfim_d_2 = PFIM('stdin_2.r')
end9 = proc.time() - start9
end9


# Web Appendix B
# Five group analyses
start10 = proc.time()
pfim_dc_5 = PFIM('stdin_cov_5.r')
end10 = proc.time() - start10
end10

start11 = proc.time()
pfim_dc_6 = PFIM('stdin_5.r')
end11 = proc.time() - start11
end11

