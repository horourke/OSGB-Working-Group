################################################################
################################################################
################################################################
plot.pdm = function(data, main, breaks = 21, diagonal = FALSE,
                    cex = 0.05, entryrange = NULL, cex.main){
  .d = dim(data)[1]
  
  
  if( is.null(entryrange)){
    entryrange = 0
    if(diagonal){
      entryrange = max(abs(data))
    }
    else{
      for(.i in 1:.d){
        for(.j in (1:.d)[-.i]){
          entryrange = max(entryrange, abs(data[.i,.j]))
        }
      }
    }
  }
  
  rbPal <- colorRampPalette(c('red',"gray95",'blue'))
  color <- rbPal(breaks)
  plot(1,1,
       col = "white",
       xlim = c(1 - 0.15*cex, .d+ 0.15*cex),
       ylim = c(1 - 0.15*cex, .d+ 0.15*cex),
       xlab = "Row Index",
       ylab = "Column Index",
       xaxt = "none",
       yaxt = "none",
       main = main,
       cex.main = cex.main)
  #axis(2, at = c(seq(from = 1, to = .d, length.out = 5)), 
  #     labels = c(seq(from = .d, to = 1, length.out = 5)))
  #axis(1, at = seq(from = 1, to = .d, length.out = 5), 
  #     labels = c(seq(from = 1, to = .d, length.out = 5)))
  
  for(.i in 1:.d){
    if(diagonal){
      for(.j in 1:.d){
        colorindex = floor( breaks*(data[.i,.j] + entryrange)/(2*entryrange) + 0.99 )
        points(x= c(.j), y= c(.d-.i+1), col = color[colorindex], 
               pch = 19, cex = cex)
      }
    } else {
      for(.j in (1:.d)[-.i]){
        colorindex = floor( breaks*(data[.i,.j] + entryrange)/(2*entryrange) + 0.99 )
        points(x= c(.j), y= c(.d-.i+1), col = color[colorindex], 
               pch = 19, cex = cex)
      }
    }
  }
}


################################################################
################################################################

cbPalette <- c("#999999", "#E69F00", "#56B4E9", 
               "#009E73", "#F0E442", "#0072B2", 
               "#D55E00", "#CC79A7")
d <- 5
mag <- runif(d^2, 3,4) * (2 * rbinom(d^2,1,0.5) - 1)


MJ <- mag * matrix(
  c(0,1,1,0,0,
    1,0,0,0,0,
    0,1,0,1,0,
    0,0,0,0,0,
    0,0,1,0,0), ncol = 5, byrow = FALSE)
MM1 <- mag * matrix(
  c(0,0,0,0,0,
    0,0,1,1,0,
    0,0,0,0,0,
    0,1,0,0,0,
    0,0,0,0,0), ncol = 5, byrow = FALSE)
MM2 <- mag * matrix(
  c(0,0,0,0,0,
    0,0,0,0,0,
    0,0,0,0,0,
    1,0,0,0,1,
    1,0,0,1,0), ncol = 5, byrow = FALSE)


M1 <- (MJ + MM1) 
M2 <- (MJ + MM1 + MM2)
M3 <- (MJ) 
M4 <- (MJ + MM2)


#############################
#############################
#############################
## Networks of N subjects.
layout(matrix(c(1,2,3,4,5), nrow = 1), widths = c(2,2,2,1,2))
par(mar = c(3,2,5,2))
plot.pdm(
  M1, main = "GCN Subject 1", breaks = 21, diagonal = TRUE,
  cex = 3, entryrange = c(-4,4), cex.main = 2.5)
plot.pdm(
  M2, main = "GCN Subject 2", breaks = 21, diagonal = TRUE,
  cex = 3, entryrange = c(-4,4), cex.main = 2.5)
plot.pdm(
  M3, main = "GCN Subject 3", breaks = 21, diagonal = TRUE,
  cex = 3, entryrange = c(-4,4), cex.main = 2.5)

par(mar = c(0,0,0,0))
plot.new()
text(0.5, 0.5, "...", cex = 5)

par(mar = c(3,2,5,2))
plot.pdm(
  M4, main = "GCN Subject n", breaks = 21, diagonal = TRUE,
  cex = 3, entryrange = c(-4,4), cex.main = 2.5)




#############################
#############################
#############################
## Decomposition of Subject 2 network:
layout(matrix(c(1:7), nrow = 1), widths = c(3,1,3,1.5,3,1.5,3))
par(mar = c(3,2,5,2))
plot.pdm(
  M2, main = "Subject GCN", breaks = 21, diagonal = TRUE,
  cex = 3, entryrange = c(-4,4), cex.main = 2.5)
par(mar = c(0,0,0,0))
plot.new()
text(0.5, 0.5, "=", cex = 5)
par(mar = c(3,2,5,2))

plot.pdm(
  M3, main = "Shared Structure", breaks = 21, diagonal = TRUE,
  cex = 3, entryrange = c(-4,4), cex.main = 2.5)
par(mar = c(0,0,0,0))
plot.new()
text(0.5, 0.5, expression(paste("+ ",italic(x[1]))), cex = 5)
par(mar = c(3,2,5,2))

plot.pdm(
  MM1, main = "Moderator 1", breaks = 21, diagonal = TRUE,
  cex = 3, entryrange = c(-4,4), cex.main = 2.5)

par(mar = c(0,0,0,0))
plot.new()
text(0.5, 0.5, expression(paste("+ ",italic(x[2]))), cex = 5)
par(mar = c(3,2,5,2))

plot.pdm(
  MM2, main = "Moderator 2", breaks = 21, diagonal = TRUE,
  cex = 3, entryrange = c(-4,4), cex.main = 2.5)





#############################
#############################
#############################
## Subject 2 network
par(mfrow = c(1,1), mar = c(2,2,2,2))
plot.pdm(
  M2, main = "", breaks = 21, diagonal = TRUE,
  cex = 3, entryrange = c(-4,4), cex.main = 2.5)
