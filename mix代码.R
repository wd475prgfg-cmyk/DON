getwd() #得到目前的路径
setwd("D:/研究生材料/数据/FT数据/溯源FT/转换溯源/贝叶斯/去除异常点/mix") #设置当前路径
.libPaths("D:/Rbao")
library(MixSIAR)
# Load mix data
mix <- load_mix_data(filename="consumer_DOM.csv",
                     iso_names=c("Intensity",  "N_C"), # 每个列名作为单独的同位素
                     factors="id",
                     fac_random=FALSE,
                     fac_nested=FALSE,
                     cont_effects=NULL)

# Load source data
source <- load_source_data(filename="source_DOM.csv",
                           source_factors=NULL,
                           conc_dep=FALSE,
                           data_type="means",
                           mix)
library(ggplot2)
# Load discrimination/TDF data
discr <- load_discr_data(filename="discirmination_DOM.csv", mix)

# Make isospace plot
plot_data(filename="isospace_plot",
          plot_save_pdf=TRUE,
          plot_save_png=TRUE,
          mix,source,discr)

# Calculate standardized convex hull area
library(splancs)
if(mix$n.iso==2) calc_area(source=source,mix=mix,discr=discr)

# Plot your prior
#plot_prior(alpha.prior=1,source)


# Define model structure and write JAGS model file
model_filename <- "MixSIAR_model.txt"
resid_err <- FALSE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

set.seed(12345) # 你可以使用任何数字作为种子，保证每一次的输出结果一致
# Run the JAGS model ("test" first, then "short")
# jags.1 <- run_model(run="test", mix, source, discr, model_filename, alpha.prior=1,resid_err,process_err)
jags.1 <- run_model(run="normal", mix, source, discr, model_filename, alpha.prior=1)


#设置最大输出行数，当输出过多数据的时候需要设置
options(max.print=10000)
# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.1, mix, source)

###清空环境中的图片
graphics.off()