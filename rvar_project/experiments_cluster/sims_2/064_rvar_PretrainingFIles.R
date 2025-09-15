
source("002_requirements_lite.R")

index  <- 1
main_folder <- "400_rvar/"
runindex <- "14_1"
load(paste0(main_folder, "pretrainings1/data/output14_1.RData"))
output <- get(paste0("output", runindex))

attach(output)

names(output)
output$phi_list

phiest_list <- output$BICmodel$rvar_opt_coeffs %>%
    select(-lambda1, -lambda2, -var) %>% 
    {lapply(0:2, function(x_val, data) {
        print(colnames(data))
        print(paste0("x", x_val))
        which_keep <- str_detect(colnames(data), paste0("x", x_val))
        data %>% select_if(which_keep) %>% 
        as.matrix() %>%
        return()
    }, data = .)}

phiest_list

par(mfcol = c(3,2), mar = c(5.1, 4.1, 4.1, 4.1))
lapply(output$phi_list, plot, breaks = 10)
lapply(phiest_list, plot, breaks = 10)
dev.off()


phiest_list <- output$CVmodel$rvar_opt_coeffs %>%
    select(-lambda1, -lambda2, -var) %>% 
    {lapply(0:2, function(x_val, data) {
        print(colnames(data))
        print(paste0("x", x_val))
        which_keep <- str_detect(colnames(data), paste0("x", x_val))
        data %>% select_if(which_keep) %>% 
        as.matrix() %>%
        return()
    }, data = .)}

phiest_list

par(mfcol = c(3,2))
lapply(output$phi_list, plot)
lapply(phiest_list, plot)
dev.off()
