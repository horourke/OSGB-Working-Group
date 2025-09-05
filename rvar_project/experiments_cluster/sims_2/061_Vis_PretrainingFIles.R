
source("002_requirements_lite.R")

index  <- 1
main_folder <- "400_rvar/"
load(paste0(main_folder, "pretrainings1/data/output1_4.RData"))


attach(output1_4)

names(output1_4)
output1_4$phi_list

phiest_list <- output1_4$BICmodel$rvar_opt_coeffs[c(3:10, 1:2), ] %>%
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
lapply(output1_4$phi_list, plot)
lapply(phiest_list, plot)
dev.off()


phiest_list <- output1_4$CVmodel$rvar_opt_coeffs[c(3:10, 1:2), ] %>%
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
lapply(output1_4$phi_list, plot)
lapply(phiest_list, plot)
dev.off()
