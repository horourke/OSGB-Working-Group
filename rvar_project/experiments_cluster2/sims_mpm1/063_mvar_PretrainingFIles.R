
source("002_requirements_lite.R")

index  <- 1
main_folder <- "300_mvar/"
runindex <- "14_1"
load(paste0(main_folder, "pretrainings1/data/output14_1.RData"))
output <- get(paste0("output", runindex))

names(output)
output$mvar_adaptive$mats$total

par(mfrow = c(3,2), mar = c(5.1, 4.1, 4.1, 4.1))
for (i in 1:3) {
    B_true <- output$data$B_list[[i]]
    B_est  <- output$mvar_adaptive$mats$total[[i]]


    plot(B_true, breaks = 10)
    plot(as.matrix(B_est), breaks = 10)
}
dev.off()



par(mfrow = c(3,2), mar = c(5.1, 4.1, 4.1, 4.1))
for (i in 1:3) {
    B_true <- output$data$B_list[[i]]
    B_est  <- output$mvar_standard$mats$total[[i]]


    plot(B_true, breaks = 10)
    plot(as.matrix(B_est), breaks = 10)
}
dev.off()


