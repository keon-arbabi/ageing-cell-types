library(png)

setwd("C:/Users/arbabik/Downloads/CAMH R/ageing-cell-types/")
paths = c("./GSE30272/MGP-GSE30272_files/figure-html", 
          "./GSE36192/MGP-GSE36192_files/figure-html",
          "./GSE60862/MGP-GSE60862_files/figure-html")


par(mar=rep(0,4)) # no margins
layout(matrix(1:6, ncol=2, byrow=TRUE))

for(n in 1:3){
(files = list.files(path = paths[n], pattern = "*.png", all.files = T, full.names = T))
(files = files[(length(files)-1):(length(files))])

  for(i in 1:2){
    img = readPNG(files[i])
    
    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
    rasterImage(img,0,0,1,1)
  }
}

# write to PDF
dev.print(pdf, "summary comparisons.pdf")
