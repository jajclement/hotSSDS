satCurveHS <- function(fIN = 'satCurve.tab',
                       sampleName = 'all'){
  
  library(ggplot2)
  library(gridExtra)
  
  options(scipen=999)

  theme_set(theme_bw(base_family = 'mono', base_size=9))

  s <- read.table(fIN,header=TRUE)
  
  s$reads <- s$reads/1000000
  
  nMax <- max(s$reads*10)
  
  if (nMax > 1.5){
    nu   <- data.frame(reads = c(seq(0.05,.95,0.05),
                                 seq(1,nMax,0.5)))
  }else{
    nu   <- data.frame(reads = c(seq(0,15,0.05)))
  }
  
  pLog <- lm(hs ~ log(reads), data = s)
  prLog <- predict(pLog, nu, se=TRUE)
  
  nu$fit   <- prLog$fit
  nu$upper <- prLog$fit + qt(0.975,prLog$df)*prLog$se.fit
  nu$lower <- prLog$fit -qt(0.975,prLog$df)*prLog$se.fit
  
  gL <- ggplot(data = s,aes(x=reads,y=hs)) + 
    geom_ribbon(data = nu, 
                aes(x=reads,y=fit,ymin=lower,ymax=upper),
                fill='orange', alpha=.5) +
    geom_line(data = nu, 
              aes(x=reads,y=fit),
              linetype='dashed') +
    geom_point(shape=21,fill='darkorange',color='grey50',size=4) + 
    scale_x_log10() + 
    annotation_logticks(sides='b') + 
    xlab('ssDNA Fragments (Million)') + 
    ylab('Hotspots') + 
    ggtitle(paste0('Log saturation curve (',sampleName,')')) 
      
  gN <- ggplot(data = s,aes(x=reads,y=hs)) + 
    geom_ribbon(data = nu, 
                aes(x=reads,y=fit,ymin=lower,ymax=upper),
                fill='orange', alpha=.5) +
    geom_line(data = nu, 
              aes(x=reads,y=fit),
              linetype='dashed') +
    geom_point(shape=21,fill='darkorange',color='grey50',size=4) + 
    xlab('ssDNA Fragments (Million)') + 
    ylab('Hotspots') + 
    ggtitle(paste0('Saturation curve (',sampleName,')'))
  
  png(paste0(sampleName,'.saturationCurve.png'), 
      res = 300, 
      width = 8, 
      height = 8, 
      units='in')
  
  grid.arrange(gN,gL,nrow=2)
  dev.off()
  
  cairo_pdf(paste0(sampleName,'.saturationCurve.pdf'),
      family="Sans",
      width = 8,
      height = 8)

  grid.arrange(gN,gL,nrow=2)
  dev.off()

  nOut1 <- data.frame(reads=as.numeric(s$reads*1000000),HS=s$hs,type='actual')
  nOut2 <- data.frame(reads=as.numeric(nu$reads*1000000),HS=round(nu$fit),type='predicted')
  nOut <- rbind(nOut1,nOut2)
  
  write.table(nOut,file = paste0(sampleName,'.preditedHotspotsFromMoreSequencing.tab'), row.names = FALSE, quote = FALSE, sep="\t")
  
}
