library(HMDHFDplus)
library(glue)
# Canada, USA, Norway, Republic of Korea

# Canada

Countries = c("CAN","USA","JPN","NOR","KOR","AUS")

color <- c(
  "CAN" = "#1b9e77",  # teal green
  "USA" = "#d95f02",  # orange
  "JPN" = "#7570b3",  # purple
  "NOR" = "#e7298a",  # pink/magenta
  "KOR" = "#66a61e",  # olive green
  "AUS" = "#e6ab02"   # mustard yellow
)

for (country in Countries){
  
  HMD_CAN = readHMDweb(CNTRY = country,item = "Deaths_1x10",username = "furkanb@my.yorku.ca",
                       password = "",fixup = TRUE)
  
  HMD_CAN_2010 = subset(HMD_CAN,Year == 2010)
  
  HMD_CAN_2010 = HMD_CAN_2010[1:(nrow(HMD_CAN_2010)-1),]
  
  HMD_CAN_2010_M = HMD_CAN_2010[,c("Year","Age","Male")]
  HMD_CAN_2010_F = HMD_CAN_2010[,c("Year","Age","Female")]
  
  HMD_CAN_2010_M$Male = round(HMD_CAN_2010_M$Male)
  HMD_CAN_2010_F$Female = round(HMD_CAN_2010_F$Female)
  
  HMD_CAN_2010_M_a_i = HMD_CAN_2010_M$Age[1:(nrow(HMD_CAN_2010_M))]
  
  grid = c(HMD_CAN_2010_M_a_i,110)
  n_male = HMD_CAN_2010_M$Male
  n_female = HMD_CAN_2010_F$Female
  
  lcd_results_male = get_fhatn_with_n(n_i = n_male,grid = grid)
  
  lcd_results_male_mu = round(EMemp_with_n(n_male,grid)$mu_hat)
  
  grid2 <- seq(min(grid),max(grid), length.out=1001)
  
  density_male = evaluateLogConDens(grid2,lcd_results_male$fhatn,which = c(4,5))
  
  cdf_male = density_male[,6]
  
  estimated_grid_density_male = density_male[,5]
  
  lcd_results_female = get_fhatn_with_n(n_female,grid)
  
  lcd_results_female_mu = round(EMemp_with_n(n_female,grid)$mu_hat)
  
  density_female = evaluateLogConDens(grid2,lcd_results_female$fhatn,which = c(4,5))
  
  cdf_female = density_female[,6]
  
  estimated_grid_density_female = density_female[,5]
  
  hazard_func_male = estimated_grid_density_male / (1 - cdf_male)
  hazard_func_female = estimated_grid_density_female / (1 - cdf_female)
    
  grDevices::pdf(
    glue("/Users/furkandanisman/Desktop/RA_Documents/RA-DOCUMENTS/Mortaility_{country}.pdf"),
    width = 12, height = 10
  )
  
  par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=6,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=3,cex.axis=3,
      font.axis=1,cex.main=3)
  
  plot(grid2, estimated_grid_density_female, type = "l",
       xlab = "", ylab = "",
       main = "",col="purple",lwd = 5)
  
  lines(grid2, estimated_grid_density_male,col="brown3",lwd = 5)

  legend("topleft", 
         legend = c(glue("Female ({lcd_results_female_mu})"),glue("Male ({lcd_results_male_mu})")), 
         col = c("purple","brown3"), 
         lty = c(1,1),
         lwd = c(5,5),
         bty = "n", 
         pt.cex = 2, 
         cex = 2, 
         text.col = "black", 
         horiz = F, 
         inset = c(0, 0),x.intersp = 0.2, y.intersp = 1,seg.len = 0.75) # -0.3
  
  dev.off()
  
  
  grDevices::pdf(
    glue("/Users/furkandanisman/Desktop/RA_Documents/RA-DOCUMENTS/Hazard_{country}.pdf"),
    width = 12, height = 10
  )
  
  par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=6,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=3,cex.axis=3,
      font.axis=1,cex.main=3)
  
  plot(grid2, hazard_func_male, type = "l",
       xlab = "", ylab = "",
       main = "",col="brown3",lwd = 5)
  
  lines(grid2, hazard_func_female,col="purple",lwd = 5)
  
  dev.off()
 
  print(country)
  
}

d = 1 

color <- c(
  "CAN" = "maroon",  # teal green
  "USA" = "#d95f02",  # orange
  "JPN" = "#7570b3",  # purple
  "NOR" = "#e7298a",  # pink/magenta
  "KOR" = "#66a61e",  # olive green
  "AUS" = "#e6ab02"   # mustard yellow
)


Countries = c("NOR","USA","JPN","CAN","KOR","AUS")


for (country in Countries){
  
  HMD_CAN = readHMDweb(CNTRY = country,item = "Deaths_1x10",username = "furkanb@my.yorku.ca",
                       password = "",fixup = TRUE)
  
  HMD_CAN_2010 = subset(HMD_CAN,Year == 2010)
  
  HMD_CAN_2010 = HMD_CAN_2010[1:(nrow(HMD_CAN_2010)-1),]
  
  HMD_CAN_2010$Total = round(HMD_CAN_2010$Total)

  HMD_CAN_2010_a_i = HMD_CAN_2010$Age[1:(nrow(HMD_CAN_2010))]
  
  grid = c(HMD_CAN_2010_a_i,110)
  n_all = HMD_CAN_2010$Total
  
  lcd_results = Log_Con_Discr_Algorithm_method3(n_all,grid)
  estimated_grid_density = evaluateLogConDens(grid,lcd_results$fhatn,which = 4)[,5]

  if(country == "NOR"){
    
    grDevices::pdf(
      glue("/Users/furkandanisman/Desktop/RA_Documents/RA-DOCUMENTS/Mortaility_all.pdf"),
      width = 15, height = 10
    )
    
    par(bg='white', mar = c(4, 5, 4, 1), xpd=T,cex=6,bty="l",font.lab=1,mfrow=c(1,1),cex.lab=3,cex.axis=3,
        font.axis=1,cex.main=3)
    
    plot(grid, estimated_grid_density, type = "l",
         xlab = "Age", ylab = "PDF",
         main = "",col=color[d],lwd = 5)
    
    d = d + 1 
    
  }else{
  
    lines(grid, estimated_grid_density,col=color[d],lwd = 5)
    
    d = d + 1
      
  }
  
}



legend("topleft", 
       legend = c("Norway","Japan","Australia","Canada","Republic of Korea","USA"), 
       col = color, 
       lty = c(1,1,1,1,1,1),
       lwd = c(3,3,3,3,3,3),
       bty = "n", 
       pt.cex = 2, 
       cex = 2, 
       text.col = "black", 
       horiz = F, 
       inset = c(0, 0),x.intersp = 0.2, y.intersp = 1,seg.len = 0.75) # -0.3


dev.off()



# Setup
library(grDevices)
library(glue)

d <- 1

color <- c(
  "CAN" = "maroon",
  "USA" = "#d95f02",
  "JPN" = "#7570b3",
  "NOR" = "#e7298a",
  "KOR" = "#66a61e",
  "AUS" = "#e6ab02"
)

# --- NOR ---
HMD_NOR <- readHMDweb(CNTRY = "NOR", item = "Deaths_1x10", username = "furkanb@my.yorku.ca",
                      password = "", fixup = TRUE)
HMD_NOR_2010 <- subset(HMD_NOR, Year == 2010)
HMD_NOR_2010 <- HMD_NOR_2010[1:(nrow(HMD_NOR_2010) - 1), ]
HMD_NOR_2010$Total <- round(HMD_NOR_2010$Total)
grid_NOR <- c(HMD_NOR_2010$Age, 110)
n_all_NOR <- HMD_NOR_2010$Total
lcd_results_NOR <- get_fhatn_with_n(n_all_NOR, grid_NOR)

grid2 <- seq(min(grid),max(grid), length.out=1001)
estimated_grid_density_NOR <- evaluateLogConDens(grid2,lcd_results_NOR$fhatn,which = 4)[,5]

lcd_results_NOR_mu <- EMemp_with_n(n_all_NOR, grid_NOR)

# Plot the first country
grDevices::pdf(
  glue("/Users/furkandanisman/Desktop/RA_Documents/RA-DOCUMENTS/Mortaility_all_2.pdf"),
  width = 15, height = 10
)

par(bg = 'white', mar = c(4, 5, 4, 1), xpd = TRUE, cex = 6, bty = "l", font.lab = 1,
    mfrow = c(1, 1), cex.lab = 3, cex.axis = 3, font.axis = 1, cex.main = 3)

plot(grid2, estimated_grid_density_NOR, type = "l",
     xlab = "Age", ylab = "PDF",
     main = "", col = color["NOR"], lwd = 5)
d <- d + 1


# --- USA ---
HMD_USA <- readHMDweb(CNTRY = "USA", item = "Deaths_1x10", username = "furkanb@my.yorku.ca",
                      password = "", fixup = TRUE)
HMD_USA_2010 <- subset(HMD_USA, Year == 2010)
HMD_USA_2010 <- HMD_USA_2010[1:(nrow(HMD_USA_2010) - 1), ]
HMD_USA_2010$Total <- round(HMD_USA_2010$Total)
grid_USA <- c(HMD_USA_2010$Age, 110)
n_all_USA <- HMD_USA_2010$Total
lcd_results_USA <- get_fhatn_with_n(n_all_USA, grid_USA)

grid2 <- seq(min(grid),max(grid), length.out=1001)

estimated_grid_density_USA <- evaluateLogConDens(grid2,lcd_results_USA$fhatn,which = 4)[,5]
lines(grid2, estimated_grid_density_USA, col = color["USA"], lwd = 5)

lcd_results_USA_mu <- EMemp_with_n(n_all_USA, grid_USA)

d <- d + 1

# --- JPN ---
HMD_JPN <- readHMDweb(CNTRY = "JPN", item = "Deaths_1x10", username = "furkanb@my.yorku.ca",
                      password = "", fixup = TRUE)
HMD_JPN_2010 <- subset(HMD_JPN, Year == 2010)
HMD_JPN_2010 <- HMD_JPN_2010[1:(nrow(HMD_JPN_2010) - 1), ]
HMD_JPN_2010$Total <- round(HMD_JPN_2010$Total)
grid_JPN <- c(HMD_JPN_2010$Age, 110)
n_all_JPN <- HMD_JPN_2010$Total
lcd_results_JPN <- get_fhatn_with_n(n_all_JPN, grid_JPN)

grid2 <- seq(min(grid),max(grid), length.out=1001)

estimated_grid_density_JPN <- evaluateLogConDens(grid2,lcd_results_JPN$fhatn,which = 4)[,5]
lines(grid2, estimated_grid_density_JPN, col = color["JPN"], lwd = 5)

lcd_results_JPN_mu <- EMemp_with_n(n_all_JPN, grid_JPN)

d <- d + 1

# --- CAN ---
HMD_CAN <- readHMDweb(CNTRY = "CAN", item = "Deaths_1x10", username = "furkanb@my.yorku.ca",
                      password = "", fixup = TRUE)
HMD_CAN_2010 <- subset(HMD_CAN, Year == 2010)
HMD_CAN_2010 <- HMD_CAN_2010[1:(nrow(HMD_CAN_2010) - 1), ]
HMD_CAN_2010$Total <- round(HMD_CAN_2010$Total)
grid_CAN <- c(HMD_CAN_2010$Age, 110)
n_all_CAN <- HMD_CAN_2010$Total
lcd_results_CAN <- get_fhatn_with_n(n_all_CAN, grid_CAN)

grid2 <- seq(min(grid),max(grid), length.out=1001)

estimated_grid_density_CAN <- evaluateLogConDens(grid2,lcd_results_CAN$fhatn,which = 4)[,5]
lines(grid2, estimated_grid_density_CAN, col = color["CAN"], lwd = 5)

lcd_results_CAN_mu <- EMemp_with_n(n_all_CAN, grid_CAN)

d <- d + 1

# --- KOR ---
HMD_KOR <- readHMDweb(CNTRY = "KOR", item = "Deaths_1x10", username = "furkanb@my.yorku.ca",
                      password = "", fixup = TRUE)
HMD_KOR_2010 <- subset(HMD_KOR, Year == 2010)
HMD_KOR_2010 <- HMD_KOR_2010[1:(nrow(HMD_KOR_2010) - 1), ]
HMD_KOR_2010$Total <- round(HMD_KOR_2010$Total)
grid_KOR <- c(HMD_KOR_2010$Age, 110)
n_all_KOR <- HMD_KOR_2010$Total
lcd_results_KOR <- get_fhatn_with_n(n_all_KOR, grid_KOR)

grid2 <- seq(min(grid),max(grid), length.out=1001)

estimated_grid_density_KOR <- evaluateLogConDens(grid2,lcd_results_KOR$fhatn,which = 4)[,5]
lines(grid2, estimated_grid_density_KOR, col = color["KOR"], lwd = 5)

lcd_results_KOR_mu <- EMemp_with_n(n_all_KOR, grid_KOR)

d <- d + 1

# --- AUS ---
HMD_AUS <- readHMDweb(CNTRY = "AUS", item = "Deaths_1x10", username = "furkanb@my.yorku.ca",
                      password = "", fixup = TRUE)
HMD_AUS_2010 <- subset(HMD_AUS, Year == 2010)
HMD_AUS_2010 <- HMD_AUS_2010[1:(nrow(HMD_AUS_2010) - 1), ]
HMD_AUS_2010$Total <- round(HMD_AUS_2010$Total)
grid_AUS <- c(HMD_AUS_2010$Age, 110)
n_all_AUS <- HMD_AUS_2010$Total
lcd_results_AUS <- get_fhatn_with_n(n_all_AUS, grid_AUS)

grid2 <- seq(min(grid),max(grid), length.out=1001)

estimated_grid_density_AUS <- evaluateLogConDens(grid2,lcd_results_AUS$fhatn,which = 4)[,5]
lines(grid2, estimated_grid_density_AUS, col = color["AUS"], lwd = 5)

lcd_results_AUS_mu <- EMemp_with_n(n_all_AUS, grid_AUS)

d <- d + 1

round(lcd_results_NOR_mu)
round(lcd_results_JPN_mu)
round(lcd_results_AUS_mu)
round(lcd_results_CAN_mu)
round(lcd_results_KOR_mu)
round(lcd_results_USA_mu)

color_2 <- c(
  "NOR" = "#e7298a",
  "JPN" = "#7570b3",
  "AUS" = "#e6ab02",
  "CAN" = "maroon",
  "USA" = "#d95f02",
  "KOR" = "#66a61e"
)


legend("topleft", 
       legend = c("Norway (79)","Japan (80)","Australia (77)","Canada (76)","USA (74)","Republic of Korea (73)"), 
       col = color_2, 
       lty = c(1,1,1,1,1,1),
       lwd = c(3,3,3,3,3,3),
       bty = "n", 
       pt.cex = 2, 
       cex = 2, 
       text.col = "black", 
       horiz = F, 
       inset = c(0, 0),x.intersp = 0.2, y.intersp = 1,seg.len = 0.75) # -0.3

# Close the PDF device
dev.off()


