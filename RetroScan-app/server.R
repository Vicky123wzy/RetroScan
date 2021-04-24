library(shiny)
library(DT)
library(shinydashboard)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(UpSetR)
library(Biostrings)
library(muscle)
library(ape)
library(ggmsa)
library(patchwork)
library(pheatmap)

options(shiny.maxRequestSize = 900*1024^2)

retrocopy_structure <- function(final,info,retrocopyID,utr_color,exon_color,intron_color,frame_color,text_size,title_size){
  retro <- final[final$Retrocopy_ID %in% retrocopyID,"Retrocopy_ID"]
  #first:parent_exon & parent_cds
  parent_exon <- info[info$mRNA %in% final[final$Retrocopy_ID %in% retrocopyID,"Parental_gene_ID"],"exon_site"]
  parent_exon <- sort(as.numeric(unlist(strsplit(as.character(parent_exon), split = ",|-"))))
  parent_start <- min(parent_exon)
  parent_exon <-parent_exon-parent_start+1
  parent_cds <- info[info$mRNA %in% final[final$Retrocopy_ID %in% retrocopyID,"Parental_gene_ID"],"cds_site"]
  parent_cds <- sort(as.numeric(unlist(strsplit(as.character(parent_cds), split = ",|-"))))
  parent_cds <-parent_cds-parent_start+1
  parent_len <- max(parent_exon)
  
  #second:pro_exon & pro_cds
  pro_exon <- 1
  for(j in 1:(length(parent_exon)/2)){
    part <- parent_exon[j*2] - parent_exon[j*2-1] + pro_exon[length(pro_exon)]
    pro_exon <- c(pro_exon,part)
  }
  pro_start <- 0
  for(j in 1:(length(parent_exon)/2)){
    if(min(parent_cds) > parent_exon[j*2]){pro_start <- parent_exon[j*2] - parent_exon[j*2-1] + 1 + pro_start}
    if(min(parent_cds) >= parent_exon[j*2-1] & min(parent_cds) <= parent_exon[j*2]){pro_start <- min(parent_cds) - parent_exon[j*2-1] + 1 + pro_start}}
  pro_cds <- pro_start
  for(j in 1:(length(parent_cds)/2)){
    part <- parent_cds[j*2] - parent_cds[j*2-1] + pro_cds[length(pro_cds)]
    pro_cds <- c(pro_cds,part)
  }
  
  #third:retro_exon & retro_cds
  if(is.na(final[final$Retrocopy_ID %in% retrocopyID,"Host_gene_ID"])==F){
    host_exon <- info[info$mRNA %in% final[final$Retrocopy_ID %in% retrocopyID,"Host_gene_ID"],"exon_site"]
    host_exon <- sort(as.numeric(unlist(strsplit(as.character(host_exon), split = ",|-"))))
    host_cds <- info[info$mRNA %in% final[final$Retrocopy_ID %in% retrocopyID,"Host_gene_ID"],"cds_site"]
    host_cds <- sort(as.numeric(unlist(strsplit(as.character(host_cds), split = ",|-"))))
    all_start <- min(host_exon,final[final$Retrocopy_ID %in% retrocopyID,"Retro_start"])
    retro_start <- final[final$Retrocopy_ID %in% retrocopyID,"Retro_start"]-all_start+1
    retro_end <- final[final$Retrocopy_ID %in% retrocopyID,"Retro_end"]-all_start+1
    host_exon <- host_exon-all_start+1
    host_cds <- host_cds-all_start+1
    all_len <- max(host_exon,retro_end)
  }
  if(is.na(final[final$Retrocopy_ID %in% retrocopyID,"Host_gene_ID"])==T){
    retro_start <- 1
    retro_end <- final[final$Retrocopy_ID %in% retrocopyID,"Retro_end"]-final[final$Retrocopy_ID %in% retrocopyID,"Retro_start"]+1
    all_len <- retro_end
  }
  
  #plot
  plot(c(0, max(parent_len,all_len)+500), c(280, 600), axes=F, type= "n", xlab = "", ylab = "", bty="l",
       main=list(paste0("The Structure of ",final[final$Retrocopy_ID %in% retrocopyID,"Retrocopy_ID"]," and parental gene ",
                        final[final$Retrocopy_ID %in% retrocopyID,"Parental_gene_ID"]),cex=title_size))
  #first
  abline(h=550,lwd=4,col="gray")#添加一条水平直线y=3，线宽为4，颜色蓝色
  rect(min(parent_exon),545,max(parent_exon),555,col=intron_color, border = "grey")
  for(j in 1:(length(parent_exon)/2)){rect(parent_exon[j*2-1], 545, parent_exon[j*2], 555, col = utr_color, border = "grey")}
  for(j in 1:(length(parent_cds)/2)){rect(parent_cds[j*2-1], 540, parent_cds[j*2], 560, col = exon_color, border = "grey")}
  text(0,570,paste0("Parental gene: ",final[final$Retrocopy_ID %in% retrocopyID,"Parental_gene_ID"],"(",
                    final[final$Retrocopy_ID %in% retrocopyID,"Parent_chr"],":",
                    final[final$Retrocopy_ID %in% retrocopyID,"Parent_start"],"-",
                    final[final$Retrocopy_ID %in% retrocopyID,"Parent_end"],")(",
                    info[info$mRNA %in% final[final$Retrocopy_ID %in% retrocopyID,"Parental_gene_ID"],"mRNA_strand"],")")
       ,cex=text_size,adj=c(0,0))
  
  #second
  #abline(h=450,lwd=4,col="gray")
  for(j in 2:length(pro_exon)){rect(pro_exon[j], 445, pro_exon[j-1], 455, col = utr_color, border = "grey")}
  for(j in 2:length(pro_cds)){rect(pro_cds[j], 440, pro_cds[j-1], 460, col = exon_color, border = "grey")}
  text(0,470,"mRNA",cex=text_size,adj=c(0,0))
  
  #third
  abline(h=350,lwd=4,col="gray")
  if(is.na(final[final$Retrocopy_ID %in% retrocopyID,"Host_gene_ID"])==F){
    rect(min(host_exon),345,max(host_exon),355,col=intron_color, border = "grey")
    for(j in 1:(length(host_exon)/2)){rect(host_exon[j*2-1], 345, host_exon[j*2], 355, col = utr_color, border = "grey")}
    for(j in 1:(length(host_cds)/2)){rect(host_cds[j*2-1], 340, host_cds[j*2], 360, col = exon_color, border = "grey")}}
  text(0,325,paste0("Retrocopy: ",retrocopyID,"(",
                    final[final$Retrocopy_ID %in% retrocopyID,"Retro_chr"],":",
                    final[final$Retrocopy_ID %in% retrocopyID,"Retro_start"],"-",
                    final[final$Retrocopy_ID %in% retrocopyID,"Retro_end"],")(",
                    final[final$Retrocopy_ID %in% retrocopyID,"Retro_strand"],")")
       ,cex=text_size,adj=c(0,1))
  
  rect(retro_start, 330, retro_end, 370, border = frame_color, lwd=2)
 # text(retro_start,375,paste0((retro_end-retro_start+1)," bp"),cex=text_size,col=frame_color,,adj=c(0,0))
  arrows(min(parent_exon), 550,min(pro_exon) , 450, col = frame_color, lwd=2,angle = 15)
  arrows(max(parent_exon), 550,max(pro_exon) , 450, col = frame_color, lwd=2,angle = 15)
  
  parent_strand <- info[info$mRNA %in% final[final$Retrocopy_ID %in% retrocopyID,"Parental_gene_ID"],"mRNA_strand"]
  if(parent_strand == "+"){
    align_start <- min(pro_cds) + (final[final$Retrocopy_ID %in% retrocopyID,"Pro_start"]-1)*3
    align_end <- min(pro_cds) + final[final$Retrocopy_ID %in% retrocopyID,"Pro_end"]*3 - 1
  }
  if(parent_strand == "-"){
    align_start <- max(pro_cds) - final[final$Retrocopy_ID %in% retrocopyID,"Pro_end"]*3 + 1
    align_end <- max(pro_cds) - (final[final$Retrocopy_ID %in% retrocopyID,"Pro_start"]-1)*3
  }
  rect(align_start, 430, align_end, 470, border = frame_color, lwd=2)
  text(align_start,475,paste0((retro_end-retro_start+1)," bp"),cex=text_size,col=frame_color,,adj=c(0,0))
  arrows(align_start, 430,retro_start , 370, col = frame_color, lwd=2,angle = 15)
  arrows(align_end, 430,retro_end , 370, col = frame_color, lwd=2,angle = 15)
  legend("bottomright", c("CDS","UTR","intron","Retrocopy"),bty="n", cex=1, ncol=2,adj=0,x.intersp=0.4,
         fill=c(utr_color,exon_color,intron_color,"white"),border=c("grey","grey","grey",frame_color)) 
}

shinyServer(function(input, output){
    shinyjs::addClass(id = "navBar", class = "navbar-right")
  
  final <- reactiveValues(data = NULL)
  observe({
    req(input$file1)
    final$data <- read.table(input$file1$datapath,header = T, sep="\t",quote = "",stringsAsFactors = F)
  })
  
  info <- reactiveValues(data = NULL)
  observe({
    req(input$file2)
    info$data <- read.table(input$file2$datapath,header = T, sep="\t",quote = "",stringsAsFactors = F)
  })
  
  cds_fasta <- reactiveValues(data = NULL)
  observe({ 
    req(input$cds_file)
    temp <- readDNAStringSet(input$cds_file$datapath, format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
    names(temp) <- gsub(" .*.*","",names(temp))
    cds_fasta$data <- temp
  })
  
  fpkm_value <- reactiveValues(data = NULL)
  observe({
    req(input$fpkm_file)
    fpkm_value$data <- read.table(input$fpkm_file$datapath,header = T, sep="\t",quote = "",stringsAsFactors = F)
  })
  
  phenotype <- reactiveValues(data = NULL)
  observe({   
    req(input$phenotype_file) 
    phenotype$data <- read.table(input$phenotype_file$datapath,header = F, sep="\t",quote = "",stringsAsFactors = F)
  })  
    
  #example  
  observeEvent(input$loadexample, {
    final$data <- read.table("data/final.out",header = T, sep="\t",stringsAsFactors = F)
    info$data <- read.table("data/gene_trans_info.txt",header = T, sep="\t",stringsAsFactors = F)
    cds_fasta$data <- readDNAStringSet("data/cds.fasta", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
    names(cds_fasta$data) <- gsub(" .*.*","",names(cds_fasta$data))
    fpkm_value$data <- read.table("data/all_samples.counts.txt",header = T, sep="\t",stringsAsFactors = F)
    phenotype$data <- read.table("data/phenotype.txt",header = F, sep="\t",stringsAsFactors = F)
    reset('file1')
    reset('file2')
    reset('cds_file')
    reset('fpkm_file')
    reset('phenotype_file')
    })
  
  observeEvent(input$reset, {
    final$data <- NULL
    info$data <- NULL
    cds_fasta$data <- NULL
    fpkm_value$data <- NULL
    phenotype$data <- NULL
    reset('file1')
    reset('file2')
    reset('cds_file')
    reset('fpkm_file')
    reset('phenotype_file')
  })
  
  #all_retrocopy
  summary_txt <- reactive({
    final <- final$data
    final <- as.data.frame(final)
    a <- data.frame(name="A",num=1:6)
    a$name <- as.character(a$name)
    a[1,1] <- "Retrocopy number:   "
    a[2,1] <- "Parental gene number:   "
    a[3,1] <- "Retrocopy average length:   "
    a[4,1] <- "Retrocopy average identity:   "
    a[5,1] <- "Retrocopy average coverage:   "
    a[6,1] <- "Parental gene average intron loss number:   "
    a[1,2] <- length(unique(final[,"Retrocopy_ID"]))
    a[2,2] <- length(unique(final[,"Parental_gene_ID"]))
    a[3,2] <- sum(final[,"Retro_end"] - final[,"Retro_start"] + 1) / dim(final)[1]
    a[4,2] <- sum(final[,"Identity"]) / (dim(final)[1]*100)
    a[5,2] <- sum(final[,"Coverage"]) / (dim(final)[1]*100)
    a[6,2] <- sum(final[,"Lost_intron"]) / dim(final)[1]
    a})
  
  output$summary_info <- renderTable(summary_txt(),rownames = F, colnames = F,spacing="xs", width="auto", digits=3, align="l",hover = T, bordered = F)
  
  output$dynamic <- DT::renderDataTable({
    final <- final$data
    final <- as.data.frame(final)
    datatable(final[,c("Retrocopy_ID","Retro_chr","Retro_strand","Retro_start","Retro_end","Parental_gene_ID","Identity","Coverage","Description","Host_gene_ID")],
              rownames = FALSE, filter = 'top', 
              options = list(lengthMenu = c(5, 10, 25, 50), pageLength = 10, autoWidth = TRUE,scrollX = TRUE,scrollY = TRUE))})
  
  output$Retrocopy_result<-downloadHandler(
    filename=function(){paste0("Retrocopy_result.txt")}, 
    content=function(file){write.table(final$data, file, row.names = F, sep="\t", quote=F)})
  
  chr_plot <- reactive({
    final <- final$data
    final <- as.data.frame(final)
    info <- info$data
    info <- as.data.frame(info)
    final <- final[final$Retro_chr %in% unique(info[info$mRNA_end >= 1000000,"mRNA_chr"]) & 
                          final$Parent_chr %in% unique(info[info$mRNA_end >= 1000000,"mRNA_chr"]),]
    ggplot(final,aes(x=Retro_chr)) + 
      geom_bar(aes(fill=factor(Parent_chr)),position="stack",alpha=0.75) +
      labs(title="The chromosome distribution of retrocopies and parental genes",x="Chromosome",y="Number of retrocopies",fill="Chr") + 
      theme_classic() + 
      scale_fill_manual(values = colorRampPalette(brewer.pal(brewer.pal.info[row.names(brewer.pal.info) %in% input$chr_color,1], input$chr_color))(length(unique(final$Parent_chr)))) + 
      theme(plot.title = element_text(hjust = 0.5, size=input$chr_title,face = "bold"),
            plot.margin = margin(10, 10, 10, 10, "pt"),
            axis.title=element_text(size=input$chr_xy),axis.text.y=element_text(size=input$chr_xy_text),
            axis.text.x=element_text(size=input$chr_xy_text,angle=input$chr_x_angle,vjust=0.5),
            legend.position=input$chr_legend)})
  
  output$retro_parent_chr <- renderPlot({chr_plot()})
  
  output$chrplotdown <- downloadHandler(
    filename=function(){paste0('Chr_Distribution.',input$chrplot)},
    content = function(file){
      if(input$chrplot == 'pdf'){pdf(file)}
      else if(input$chrplot == 'png'){png(file)}
      else if(input$chrplot == 'jpeg'){jpeg(file)}
      else if(input$chrplot == 'svg'){svg(file)}
      print(chr_plot())
      dev.off()})
  
  num_plot <- reactive({
    parent_num <- as.data.frame(table(as.data.frame(table(final$data$Parental_gene_ID))$Freq))
    ggplot(parent_num,aes(x=Var1,y=Freq)) + 
      geom_bar(stat="identity",fill=input$num_color) +
      labs(title="The number distribution of retrocopies owned by each parent gene",x="Number of parental genes",y="Number of retrocopies") + 
      theme_classic() + 
      theme(plot.title = element_text(hjust = 0.5, size=input$num_title,face = "bold"),
            plot.margin = margin(10, 10, 10, 10, "pt"),
            axis.title=element_text(size=input$num_xy),axis.text=element_text(size=input$num_text))})
  
  output$retro_parent_num <- renderPlot({num_plot()})
  
  output$numplotdown <- downloadHandler(
    filename=function(){paste0('Num_Distribution_between_retrocopies_parental_genes.',input$numplot)},
    content = function(file){
      if(input$numplot == 'pdf'){pdf(file)}
      else if(input$numplot == 'png'){png(file)}
      else if(input$numplot == 'jpeg'){jpeg(file)}
      else if(input$numplot == 'svg'){svg(file)}
      print(num_plot())
      dev.off()})
  
  len_plot <- reactive({
    final <- final$data
    final$Retro_len <- final$Retro_end - final$Retro_start + 1
    ggplot(final,aes(x=Retro_len)) + 
      geom_histogram(position = 'dodge',binwidth=input$len_bin,fill=input$len_color)+
      labs(title="Distribution of retrocopy length",x="Length of retrocopies",y="Number of retrocopies") + 
      theme_classic() + 
      theme(plot.title = element_text(hjust = 0.5, size=input$len_title,face = "bold"),
            plot.margin = margin(10, 10, 10, 10, "pt"),
            axis.title=element_text(size=input$len_xy),axis.text=element_text(size=input$len_text))})
  
  output$retro_len <- renderPlot({len_plot()})
  
  output$lenplotdown <- downloadHandler(
    filename=function(){paste0('Distribution_of_retrocopy_length.',input$lenplot)},
    content = function(file){
      if(input$lenplot == 'pdf'){pdf(file)}
      else if(input$lenplot == 'png'){png(file)}
      else if(input$lenplot == 'jpeg'){jpeg(file)}
      else if(input$lenplot == 'svg'){svg(file)}
      print(len_plot())
      dev.off()})
  
  iden_plot <- reactive({
    final <- final$data
    table_pie <- data.frame(type=c("50% <= Identity < 60%","60% <= Identity < 70%","70% <= Identity < 80%","80% <= Identity < 90%","90% <= Identity <= 100%"),num=0)
    table_pie[1,2] <- dim(final[final$Identity>=50 & final$Identity<60,])[1]
    table_pie[2,2] <- dim(final[final$Identity>=60 & final$Identity<70,])[1]
    table_pie[3,2] <- dim(final[final$Identity>=70 & final$Identity<80,])[1]
    table_pie[4,2] <- dim(final[final$Identity>=80 & final$Identity<90,])[1]
    table_pie[5,2] <- dim(final[final$Identity>=90,])[1]
    ggplot(data=table_pie, aes(x="type",y=num,fill=type))+
      geom_bar(stat = 'identity', position = 'stack',alpha=0.75)+
      coord_polar("y", start=0)+
      geom_text(aes(label = scales::percent(num/sum(num))), size=input$iden_text, position=position_stack(vjust = 0.5)) + 
      scale_fill_brewer(palette=input$iden_color)+
      theme_minimal()+
      labs(title="Percentage of identity",fill="")+
      theme(
        axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),
        panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank(),
        legend.position=input$iden_legend,legend.text=element_text(size=input$iden_legend_size),
        plot.title=element_text(hjust = 0.5,size=input$iden_title, face="bold"))})
  
  output$identity_plot <- renderPlot({iden_plot()})
  
  output$idenplotdown <- downloadHandler(
    filename=function(){paste0('Percentage_of_identity.',input$idenplot)},
    content = function(file){
      if(input$idenplot == 'pdf'){pdf(file)}
      else if(input$idenplot == 'png'){png(file)}
      else if(input$idenplot == 'jpeg'){jpeg(file)}
      else if(input$idenplot == 'svg'){svg(file)}
      print(iden_plot())
      dev.off()})
  
  cov_plot <- reactive({
    final <- final$data
    final$Coverage <- final$Coverage
    table_pie <- data.frame(type=c("50% <= Coverage < 60%","60% <= Coverage < 70%","70% <= Coverage < 80%","80% <= Coverage < 90%","90% <= Coverage <= 100%"),num=0)
    table_pie[1,2] <- dim(final[final$Coverage>=50 & final$Coverage<60,])[1]
    table_pie[2,2] <- dim(final[final$Coverage>=60 & final$Coverage<70,])[1]
    table_pie[3,2] <- dim(final[final$Coverage>=70 & final$Coverage<80,])[1]
    table_pie[4,2] <- dim(final[final$Coverage>=80 & final$Coverage<90,])[1]
    table_pie[5,2] <- dim(final[final$Coverage>=90,])[1]
    ggplot(data=table_pie, aes(x="type",y=num,fill=type))+
      geom_bar(stat = 'identity', position = 'stack',alpha=0.75)+
      coord_polar("y", start=0)+
      geom_text(aes(label = scales::percent(num/sum(num))), size=input$cov_text, position=position_stack(vjust = 0.5)) + 
      scale_fill_brewer(palette=input$cov_color)+
      theme_minimal()+
      labs(title="Percentage of coverage",fill="")+
      theme(
        axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),
        panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank(),
        legend.position=input$cov_legend,legend.text=element_text(size=input$cov_legend_size),
        plot.title=element_text(hjust = 0.5,size=input$cov_title, face="bold"))})
  
  output$coverage_plot <- renderPlot({cov_plot()})
  
  output$covplotdown <- downloadHandler(
    filename=function(){paste0('Percentage_of_coverage.',input$covplot)},
    content = function(file){
      if(input$covplot == 'pdf'){pdf(file)}
      else if(input$covplot == 'png'){png(file)}
      else if(input$covplot == 'jpeg'){jpeg(file)}
      else if(input$covplot == 'svg'){svg(file)}
      print(cov_plot())
      dev.off()})
  
  des_plot <- reactive({
    table_pie <- as.data.frame(table(final$data$Description))
    colnames(table_pie) <- c("type","num")
    ggplot(data=table_pie, aes(x="type",y=num,fill=type))+
      geom_bar(stat = 'identity', position = 'stack',alpha=0.75)+
      coord_polar("y", start=0)+
      geom_text(aes(label = scales::percent(num/sum(num))), size=input$des_text, position=position_stack(vjust = 0.5)) + 
      scale_fill_brewer(palette=input$des_color)+
      theme_minimal()+
      labs(title="Percentage of category",fill="")+
      theme(
        axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),
        panel.border = element_blank(),panel.grid=element_blank(),axis.ticks = element_blank(),
        legend.position=input$des_legend,legend.text=element_text(size=input$des_legend_size),
        plot.title=element_text(hjust = 0.5,size=input$des_title, face="bold"))})
  
  output$description_plot <- renderPlot({des_plot()})
  
  output$desplotdown <- downloadHandler(
    filename=function(){paste0('Percentage_of_description.',input$desplot)},
    content = function(file){
      if(input$desplot == 'pdf'){pdf(file)}
      else if(input$desplot == 'png'){png(file)}
      else if(input$desplot == 'jpeg'){jpeg(file)}
      else if(input$desplot == 'svg'){svg(file)}
      print(des_plot())
      dev.off()})
  
  venn_plot <- reactive({
    final=final$data
   part_venn <- data.frame("retro"=final$Retrocopy_ID,"Identity"=0,"Coverage"=0,"KaKs"=0,"Host_gene"=0,"Lost_intron"=0)
      rownames(part_venn) <- part_venn$retro
      part_venn[part_venn$retro %in% final[final$Identity >= input$iden_cutoff,"Retrocopy_ID"],"Identity"] <- 1
      part_venn[part_venn$retro %in% final[final$Coverage >= input$cov_cutoff,"Retrocopy_ID"],"Coverage"] <- 1
      part_venn[part_venn$retro %in% final[final$Ka.Ks <= input$kaks_cutoff,"Retrocopy_ID"],"KaKs"] <- 1
      part_venn[part_venn$retro %in% final[is.na(final$Host_gene_ID)==F,"Retrocopy_ID"],"Host_gene"] <- 1
      part_venn[part_venn$retro %in% final[final$Lost_intron >= input$intron_cutoff,"Retrocopy_ID"],"Lost_intron"] <- 1
      upset(part_venn[,2:6],text.scale = 3,sets.x.label = "",point.size = 3, line.size = 1,mainbar.y.label = "",
            matrix.color = input$venn_matrix_color, main.bar.color = input$venn_main_bar_color,
            sets.bar.color = input$venn_sets_bar_color,shade.color = input$venn_shade_color)
      
   })
  
  output$venn <- renderPlot({venn_plot()})
  
  output$vennplotdown <- downloadHandler(
    filename=function(){paste0('Percentage_of_description.',input$vennplot)},
    content = function(file){
      if(input$vennplot== 'pdf'){pdf(file)}
      else if(input$vennplot == 'png'){png(file)}
      else if(input$vennplot == 'jpeg'){jpeg(file)}
      else if(input$vennplot == 'svg'){svg(file)}
      print(venn_plot())
      dev.off()})
  
  #retrocopy_sample
  retro_txt <- reactive({
    final <- final$data
    a <- data.frame(name="A",num=1:9)
    a$name <- as.character(a$name)
    a[1,1] <- "Retrocopy ID:   "
    a[2,1] <- "Coordinates:   "
    a[3,1] <- "Strand:   "
    a[4,1] <- "Lost intron:   "
    a[5,1] <- "Coverage:   "
    a[6,1] <- "Identity:   "
    a[7,1] <- "Ka:   "
    a[8,1] <- "Ks:   "
    a[9,1] <- "Ka/ks:   "
    a[1,2] <- input$sample1
    a[2,2] <- paste0(final[final$Retrocopy_ID %in% input$sample1,"Retro_chr"],":",
               final[final$Retrocopy_ID %in% input$sample1,"Retro_start"],"-",
               final[final$Retrocopy_ID %in% input$sample1,"Retro_end"])
    a[3,2] <- final[final$Retrocopy_ID %in% input$sample1,"Retro_strand"]
    a[4,2] <- final[final$Retrocopy_ID %in% input$sample1,"Lost_intron"]
    a[5,2] <- paste0(final[final$Retrocopy_ID %in% input$sample1,"Coverage"],"%")
    a[6,2] <- paste0(final[final$Retrocopy_ID %in% input$sample1,"Identity"],"%")
    a[7,2] <- final[final$Retrocopy_ID %in% input$sample1,"Ka"]
    a[8,2] <- final[final$Retrocopy_ID %in% input$sample1,"Ks"]
    a[9,2] <- final[final$Retrocopy_ID %in% input$sample1,"Ka.Ks"]
    a})
  output$retro_info <- renderTable(retro_txt(),rownames = F, colnames = F,spacing="xs", width="auto", digits=3, align="l",hover = T, bordered = F)
  
  parent_txt <- reactive({
    final <- final$data
    info <- info$data
    a <- data.frame(name=c("Parental gene ID:   ","Parental coordinates:   ","Parental strand:   "),num=1:3)
    a$name <- as.character(a$name)
    a[1,2] <- final[final$Retrocopy_ID %in% input$sample1,"Parental_gene_ID"]
    a[2,2] <- paste0(final[final$Retrocopy_ID %in% input$sample1,"Parent_chr"],":",
                     final[final$Retrocopy_ID %in% input$sample1,"Parent_start"],"-",
                     final[final$Retrocopy_ID %in% input$sample1,"Parent_end"]) 
    a[3,2] <- info[info$mRNA %in% final[final$Retrocopy_ID %in% input$sample1,"Parental_gene_ID"],"mRNA_strand"]
    a})
  output$parent_info <- renderTable(parent_txt(),rownames = F, colnames = F,spacing="xs", width="auto", digits=3, align="l",hover = T, bordered = F)
  
  host_txt <- reactive({
    final <- final$data
    info <- info$data
    a <- data.frame(name=c("Host gene ID:   ","Host gene coordinates:   ","Host gene strand:   "), num=c("N/A","N/A","N/A"))
    a$name <- as.character(a$name)
    a$num <- as.character(a$num)
    if(is.na(final[final$Retrocopy_ID %in% input$sample1,"Host_gene_ID"])==F){
      a[1,2] <- final[final$Retrocopy_ID %in% input$sample1,"Host_gene_ID"]
      a[2,2] <- paste0(info[info$mRNA %in% final[final$Retrocopy_ID %in% input$sample1,"Host_gene_ID"],"mRNA_chr"],
                       ":",info[info$mRNA %in% final[final$Retrocopy_ID %in% input$sample1,"Host_gene_ID"],"mRNA_start"],
                       "-",info[info$mRNA %in% final[final$Retrocopy_ID %in% input$sample1,"Host_gene_ID"],"mRNA_end"]) 
      a[3,2] <- info[info$mRNA %in% final[final$Retrocopy_ID %in% input$sample1,"Host_gene_ID"],"mRNA_strand"]}
    a})
  output$host_info <- renderTable(host_txt(),rownames = F, colnames = F,spacing="xs", width="auto", digits=3, align="l",hover = T, bordered = F)

  output$retro_plot <- renderPlot({
    retrocopy_structure(final$data,info$data,input$sample1,input$utr_color,input$exon_color,input$intron_color,
                        input$frame_color,input$text_size,input$title_size)})
  
  output$retroplotdown <- downloadHandler(
    filename=function(){paste0('Structure_of_retrocopy.',input$retroplot)},
    content = function(file){
      if(input$retroplot== 'pdf'){pdf(file)}
      else if(input$retroplot == 'png'){png(file)}
      else if(input$retroplot == 'jpeg'){jpeg(file)}
      else if(input$retroplot == 'svg'){svg(file)}
      print(retrocopy_structure(final$data,info$data,input$sample1,input$utr_color,input$exon_color,input$intron_color,
                                input$frame_color,input$text_size,input$title_size))
      dev.off()})
  
  output$retroseq <- reactive({
    str = paste(">",input$sample1,sep="")
    sequence <- as.character(final$data[final$data$Retrocopy_ID %in% input$sample1,"Retro_sequence"])
    for(i in 1:ceiling(nchar(sequence)/60)){
      str <- paste(str, substr(sequence,(i-1)*60+1,i*60),sep = '\n')
    }
    str
  })
  
  output$retropep <- reactive({
      str = paste(">",input$sample1,sep="")
      sequence <- as.character(final$data[final$data$Retrocopy_ID %in% input$sample1,"Retro_protein"])
      for(i in 1:ceiling(nchar(sequence)/60)){
          str <- paste(str, substr(sequence,(i-1)*60+1,i*60),sep = '\n')
      }
      str
  })
  

  output$alignment <- renderPlot({
    final <- final$data
    info <- info$data
    cds_fasta <- cds_fasta$data
    retro_seq <- final[final$Retrocopy_ID %in% input$sample1,"Retro_sequence"]
    pep_id <- final[final$Retrocopy_ID %in% input$sample1,"Parental_gene_ID"]
    pep_seq <- cds_fasta[names(cds_fasta) %in% pep_id,]
    pep_seq <- subseq(pep_seq,(final[final$Retrocopy_ID %in% input$sample1,"Pro_start"]*3 + 1),
                      (final[final$Retrocopy_ID %in% input$sample1,"Pro_end"]*3))
    if(info[info$mRNA %in% pep_id,"mRNA_strand"] == "-"){ pep_seq <- reverseComplement(pep_seq)}
    retro_pep <- unlist(DNAStringSetList(retro_seq,pep_seq))
    names(retro_pep) <- c(input$sample1,pep_id)
    aln <- muscle::muscle(retro_pep)
  #  show(detail(aln))
    p <- ggmsa(aln, start = 0, end = 60, font = "helvetical", color = input$aln_color,char_width = 0.5,seq_name = TRUE)
    for(i in 2:ceiling(nchar(aln)/60)){
      p2 <- ggmsa(aln, start = (i-1)*60+1, end = i*60, font = "helvetical", color = input$aln_color,char_width = 0.5,seq_name = TRUE)
      p <- p/p2
    }
    p
  })
  
  mean_fpkm <- reactive({
    final <- final$data
    info <- info$data
    phenotype <- phenotype$data
    fpkm_value <- fpkm_value$data
    
    part_fpkm <- fpkm_value[fpkm_value$gene %in% c(final$Retrocopy_ID,final$Parental_gene_ID),]
    phenotype_list <- unique(phenotype$V2)
    mean_fpkm <- data.frame(gene=part_fpkm$gene)
    for(i in 1:length(phenotype_list)){
      mean_fpkm <- cbind(mean_fpkm,apply(as.matrix(part_fpkm[,colnames(part_fpkm) %in% phenotype[phenotype$V2 %in% phenotype_list[i],"V1"]]), 1, mean))
      colnames(mean_fpkm)[dim(mean_fpkm)[2]] <- phenotype_list[i]
    }
    mean_fpkm})
  
  retro_exp <- reactive({
    final <- final$data
    part_fpkm <- mean_fpkm()

    retrocopy_exp <- part_fpkm[part_fpkm$gene %in% input$sample1,][1,]
    parent_exp <- part_fpkm[part_fpkm$gene %in% final[final$Retrocopy_ID %in% input$sample1,"Parental_gene_ID"],][1,]
    if(dim(retrocopy_exp)[1]==1 && dim(parent_exp)[1]==1){
      a <- data.frame(retro_parent=c(rep("Retrocopy",dim(retrocopy_exp)[2]-1),rep("Parental gene",dim(parent_exp)[2]-1)),
                      tissue=rep(colnames(retrocopy_exp)[2:dim(retrocopy_exp)[2]],2),
                      fpkm=c(as.character(retrocopy_exp[1,2:dim(retrocopy_exp)[2]]),as.character(parent_exp[1,2:dim(parent_exp)[2]])))
      
      p=ggplot(a, aes(x=factor(tissue), y=fpkm, colour=retro_parent,group=retro_parent)) + geom_line(size=input$exp_linesize)+
        labs(title=paste0("The FPKM between retrocopy ", input$sample1, " and parental gene ",final[final$Retrocopy_ID %in% input$sample1,"Parental_gene_ID"]))+
        scale_colour_manual(values = c(input$exp_parentcolor,input$exp_retrocolor)) + 
        theme_classic() + 
        theme(plot.title = element_text(hjust = 0.5, size=input$exp_titlesize,face = "bold"),
              plot.margin = margin(10, 10, 10, 10, "pt"),
              axis.title=element_text(size=input$exp_textsize),axis.text.y=element_text(size=input$exp_textsize),
              axis.text.x=element_text(size=input$exp_textsize,angle=input$exp_xangle,vjust=0.5))
    }
   # else{
  #    P <- plot(c(0,500),c(0,600),axes=F, type= "n", xlab = "", ylab = "", bty="l")
  #    text(100,400,"mRNA",cex=3,adj=c(0,0))
   #   }
  p})
  
  output$retro_exp_plot <- renderPlot({retro_exp()})
  
  output$retroexpdown <- downloadHandler(
    filename=function(){paste0('Expression_bwteen_retrocopy_and_parental_gene.',input$retroexpplot)},
    content = function(file){
      if(input$retroexpplot == 'pdf'){pdf(file)}
      else if(input$retroexpplot == 'png'){png(file)}
      else if(input$retroexpplot == 'jpeg'){jpeg(file)}
      else if(input$retroexpplot == 'svg'){svg(file)}
      print(retro_exp())
      dev.off()})
  
  #kaks
  kaks_txt <- reactive({
    final <- final$data
    a <- data.frame(name="A",num=1:3)
    a$name <- as.character(a$name)
    a[1,1] <- "Average Ka:   "
    a[2,1] <- "Average Ks:   "
    a[3,1] <- "Average Ka/Ks:   "
    a[1,2] <- sum(final[is.na(final$Ka)==F,]$Ka)/dim(final[is.na(final$Ka)==F,])[1]
    a[2,2] <- sum(final[is.na(final$Ks)==F,]$Ks)/dim(final[is.na(final$Ks)==F,])[1]
    a[3,2] <- sum(final[is.na(final$Ka.Ks)==F,]$Ka.Ks)/dim(final[is.na(final$Ks)==F,])[1]
    a})
  
  output$kaks_info <- renderTable(kaks_txt(),rownames = F, colnames = F,spacing="xs", width="auto", digits=3, align="l",hover = T, bordered = F)
  
  output$kaks_result<-DT::renderDataTable({
    datatable(final$data[,c("Retrocopy_ID","Ka","Ks","Ka.Ks")],rownames = FALSE,filter = 'top',
              options = list(lengthMenu = c(5, 10, 25, 50),pageLength = 10, autoWidth = TRUE,scrollX = TRUE))
  })
  
  output$kaks_result_download<-downloadHandler(
    filename=function(){paste0("Retrocopy_kaks.txt")}, 
    content=function(file){write.table(final$data[,c("Retrocopy_ID","Ka","Ks","Ka.Ks")], file, row.names = F, sep="\t", quote=F)})
  
  ks_plot <- reactive({
    final <- final$data[,c("Retrocopy_ID","Ks")]
    ggplot(final,aes(x=Ks)) + 
      geom_histogram(position = 'dodge',binwidth=input$ka_bin,fill=input$ka_col)+
      labs(title="Distribution of Ks",x="Ks",y="Number of retrocopies") + 
      theme_classic() + 
      theme(plot.title = element_text(hjust = 0.5, size=input$ka_title,face = "bold"),
            plot.margin = margin(10, 10, 10, 10, "pt"),
            axis.title=element_text(size=input$ka_xy),axis.text=element_text(size=input$ka_text))})
  
  output$ks_len <- renderPlot({ks_plot()})
  
  output$ksplotdown <- downloadHandler(
    filename=function(){paste0('Distribution_of_ks.',input$kaplot)},
    content = function(file){
      if(input$kaplot == 'pdf'){pdf(file)}
      else if(input$kaplot == 'png'){png(file)}
      else if(input$kaplot == 'jpeg'){jpeg(file)}
      else if(input$kaplot == 'svg'){svg(file)}
      print(ks_plot())
      dev.off()})
  
  kaks_plot <- reactive({
    final <- final$data[,c("Retrocopy_ID","Ka.Ks")]
    ggplot(final,aes(x=Ka.Ks)) + 
      geom_histogram(position = 'dodge',binwidth=input$kaks_bin,fill=input$kaks_col)+
      labs(title="Distribution of Ka/Ks",x="Ka/Ks",y="Number of retrocopies") + 
      theme_classic() + 
      theme(plot.title = element_text(hjust = 0.5, size=input$kaks_title,face = "bold"),
            plot.margin = margin(10, 10, 10, 10, "pt"),
            axis.title=element_text(size=input$kaks_xy),axis.text=element_text(size=input$kaks_text))})
  
  output$kaks_len <- renderPlot({kaks_plot()})
  
  output$kaksplotdown <- downloadHandler(
    filename=function(){paste0('Distribution_of_kaks.',input$kaksplot)},
    content = function(file){
      if(input$kaksplot == 'pdf'){pdf(file)}
      else if(input$kaksplot == 'png'){png(file)}
      else if(input$kaksplot == 'jpeg'){jpeg(file)}
      else if(input$kaksplot == 'svg'){svg(file)}
      print(kaks_plot())
      dev.off()})
  
  ks_type <- reactive({
    final <- final$data[,c("Retrocopy_ID","Ks","Description")]
    ggplot(final,aes(x=Ks,fill=Description)) + 
      geom_bar(stat="bin",position=position_dodge(),binwidth = input$ka_type_bin) +
      labs(title="Distribution of Ks",x="Ks",y="Number of retrocopies",fill="Class") + 
      theme_classic() +
      scale_fill_manual(values = colorRampPalette(brewer.pal(4, input$ka_type_color))(4)) + 
      theme(plot.title = element_text(hjust = 0.5, size=input$ka_type_title,face = "bold"),
            plot.margin = margin(10, 10, 10, 10, "pt"),
            axis.title=element_text(size=input$ka_type_xy),axis.text.y=element_text(size=input$ka_type_xy_text),
            axis.text.x=element_text(size=input$ka_type_xy_text,vjust=0.5),
            legend.position=input$ka_type_legend)})
  
  output$ks_type_plot <- renderPlot({ks_type()})
  
  output$ks_typeplotdown <- downloadHandler(
    filename=function(){paste0('Distribution_of_ks_in_three_types.',input$ka_typeplot)},
    content = function(file){
      if(input$ka_typeplot == 'pdf'){pdf(file)}
      else if(input$ka_typeplot == 'png'){png(file)}
      else if(input$ka_typeplot == 'jpeg'){jpeg(file)}
      else if(input$ka_typeplot == 'svg'){svg(file)}
      print(ks_type())
      dev.off()})
  
  kaks_type <- reactive({
    final <- final$data[,c("Retrocopy_ID","Ka.Ks","Description")]
    ggplot(final,aes(x=Ka.Ks,fill=Description)) + 
      geom_bar(stat="bin",position=position_dodge(),binwidth = input$kaks_type_bin) +
      labs(title="Distribution of Ka/Ks",x="Ka/Ks",y="Number of retrocopies",fill="Class") + 
      theme_classic() +
      scale_fill_manual(values = colorRampPalette(brewer.pal(4, input$kaks_type_color))(4)) + 
      theme(plot.title = element_text(hjust = 0.5, size=input$kaks_type_title,face = "bold"),
            plot.margin = margin(10, 10, 10, 10, "pt"),
            axis.title=element_text(size=input$kaks_type_xy),axis.text.y=element_text(size=input$kaks_type_xy_text),
            axis.text.x=element_text(size=input$kaks_type_xy_text,vjust=0.5),
            legend.position=input$kaks_type_legend)})
  
  output$kaks_type_plot <- renderPlot({kaks_type()})
  
  output$kaks_typeplotdown <- downloadHandler(
    filename=function(){paste0('Distribution_of_kaks_in_three_types.',input$kaks_typeplot)},
    content = function(file){
      if(input$kaks_typeplot == 'pdf'){pdf(file)}
      else if(input$kaks_typeplot == 'png'){png(file)}
      else if(input$kaks_typeplot == 'jpeg'){jpeg(file)}
      else if(input$kaks_typeplot == 'svg'){svg(file)}
      print(kaks_type())
      dev.off()})
    
  
  #expression
  expr_txt <- reactive({
    mean_fpkm <- mean_fpkm()
    final <- final$data
    
    retro_fpkm <- mean_fpkm[mean_fpkm$gene %in% final$Retrocopy_ID,]
    retro_static <- data.frame(gene=retro_fpkm$gene)
    
    over_zero_num <- function(x){x <- as.numeric(x)
      length(x[x>0])}
    retro_static$num <- apply(retro_fpkm[,2:dim(retro_fpkm)[2]],1,over_zero_num)
    
    over_one_num <- function(x){x <- as.numeric(x)
    length(x[x>1])}
    retro_static$high_num <- apply(retro_fpkm[,2:dim(retro_fpkm)[2]],1,over_one_num)
    
    a <- data.frame(name="A",num=1:2)
    a$name <- as.character(a$name)
    a[1,1] <- "The number of expressed retrocopies (FPKM > 0):   "
    a[2,1] <- "The number of expressed retrocopies (FPKM > 0) in all tissues:   "
    a[3,1] <- "The number of robust expressed retrocopies (FPKM > 1):   "
    a[4,1] <- "The number of robust expressed retrocopies (FPKM > 1) in all tissues:   "
    a[1,2] <- dim(retro_static[retro_static$num > 0,])[1]
    a[2,2] <- dim(retro_static[retro_static$num == (dim(retro_fpkm)[2]-1),])[1]
    a[3,2] <- dim(retro_static[retro_static$high_num > 0,])[1]
    a[4,2] <- dim(retro_static[retro_static$high_num == (dim(retro_fpkm)[2]-1),])[1]
    a})
  
  output$expr_info <- renderTable(expr_txt(),rownames = F, colnames = F,spacing="xs", width="auto", digits=3, align="l",hover = T, bordered = F)
  
  output$expr_result<-DT::renderDataTable({
    datatable(mean_fpkm(),rownames = FALSE,filter = 'top',
              options = list(lengthMenu = c(5, 10, 25, 50),pageLength = 10, autoWidth = TRUE,scrollX = TRUE))
  })
  
  output$expr_result_download<-downloadHandler(
    filename=function(){paste0("FPKM.txt")}, 
    content=function(file){write.table(mean_fpkm(), file, row.names = F, sep="\t", quote=F)})
  
  expr_plot <- reactive({
    part_fpkm <- mean_fpkm()
    final <- final$data
    
    a<-data.frame()
    for(i in 2:dim(part_fpkm)[2]){
      b<-data.frame(class=colnames(part_fpkm)[i],factor="parent",num=1)
      b$num<-mean(part_fpkm[part_fpkm$gene %in% final$Parental_gene_ID,i])
      a<-rbind(a,b)
      b<-data.frame(class=colnames(part_fpkm)[i],factor="retrocopy",num=1)
      b$num<-mean(part_fpkm[part_fpkm$gene %in% final$Retrocopy_ID,i])
      a<-rbind(a,b)
    }
    
    ggplot(a,aes(x=class,y=num,fill=factor))+geom_col( position="dodge") +
      labs(title="The number of expressed retrocopies and parental genes",x="",y="Number of retrocopies",fill="Class") + 
      theme_classic() +
      scale_fill_manual(values = c(input$expr_parentcolor,input$expr_retrocolor)) + 
      theme(plot.title = element_text(hjust = 0.5, size=input$expr_titlesize,face = "bold"),
            plot.margin = margin(10, 10, 10, 10, "pt"),legend.position=input$expr_legend,
            axis.title=element_text(size=input$expr_text_size),axis.text.y=element_text(size=input$expr_text_size),
            axis.text.x=element_text(size=input$expr_text_size,angle=input$expr_xangle,vjust=0.5))
    })
  
  output$exprplot <- renderPlot({expr_plot()})
  
  output$expr_plotdown <- downloadHandler(
    filename=function(){paste0('The_expressed_number.',input$expr_typeplot)},
    content = function(file){
      if(input$expr_typeplot == 'pdf'){pdf(file)}
      else if(input$expr_typeplot == 'png'){png(file)}
      else if(input$expr_typeplot == 'jpeg'){jpeg(file)}
      else if(input$expr_typeplot == 'svg'){svg(file)}
      print(expr_plot())
      dev.off()})
  
  expr_heatmap <- reactive({pheatmap(t(log(mean_fpkm()[,2:dim(mean_fpkm())[2]]+1)))})
  
  output$exprheatmap <- renderPlot({expr_heatmap()})
  
  output$expr_heatdown <- downloadHandler(
    filename=function(){paste0('The_expressed_number.',input$expr_heat)},
    content = function(file){
      if(input$expr_heat == 'pdf'){pdf(file)}
      else if(input$expr_heat == 'png'){png(file)}
      else if(input$expr_heat == 'jpeg'){jpeg(file)}
      else if(input$expr_heat == 'svg'){svg(file)}
      print(expr_heatmap())
      dev.off()})
  
})