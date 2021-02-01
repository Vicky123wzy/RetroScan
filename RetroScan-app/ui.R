# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
# Developed with R version 3.3.2 (64-bit)
library(dplyr)
library(stringr)
library(png)
library(shinyjs)
library(DT)
library(visNetwork)
library(rintrojs)
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
library(colourpicker)

#source("carouselPanel.R")

# Panel div for visualization
# override the currently broken definition in shinyLP version 1.1.0
panel_div <- function(class_type, content) {
    div(class = sprintf("panel panel-%s", class_type),
        div(class = "panel-body", content)
    )
}

shinyUI(navbarPage(title = img(src="title_retro.png", height = "40px"), id = "navBar",
                   theme = "paper.css",
                   collapsible = TRUE,
                   inverse = TRUE,
                   windowTitle = "RetroScan",
                   position = "fixed-top",
                   footer = includeHTML("./www/include_footer.html"),
                   header = tags$style(
                       ".navbar-right {
                       float: right !important;
                       }",
                       "body {padding-top: 75px;}"),
                   
                   tabPanel("HOME", value = "home",
                            
                            shinyjs::useShinyjs(),
                            
                            tags$head(tags$script(HTML('
                                                       var fakeClick = function(tabName) {
                                                       var dropdownList = document.getElementsByTagName("a");
                                                       for (var i = 0; i < dropdownList.length; i++) {
                                                       var link = dropdownList[i];
                                                       if(link.getAttribute("data-value") == tabName) {
                                                       link.click();
                                                       };
                                                       }
                                                       };
                                                       '))),
                            fluidRow(
                                HTML("
                                     <section class='banner'>
                                     <h2 class='parallax'>RetroScan</h2>
                                     <p class='parallax_description'>An easy-to-use systematic pipeline for retrocopy annotation and visualization.</p>
                                     </section>
                                     ")
                                ),

                            fluidRow(shiny::HTML("<br><br><center> <h1>Ready to Get Started?</h1> </center>
                                                 <br>")
                            ),
                            fluidRow(
                                column(3),
                                column(6,
                                       tags$div(align = "center", 
                                                tags$a("Start", 
                                                       onclick="fakeClick('careerPF')", 
                                                       class="btn btn-primary btn-lg")
                                       )
                                ),
                                column(3)
                            ),
                            fluidRow(style = "height:25px;"
                            ),                                                        
                            
                            # PAGE BREAK
                            tags$hr(),
                            
                            # WHAT
                            fluidRow(
                                column(3),
                                column(6,
                                       shiny::HTML("<br><br><center> <h1>What is RetroScan</h1> </center><br>"),
                                       shiny::HTML("<h5>RetroScan is a publicly available and easy-to-use tool to scan, 
                                                   annotate and display retrocopies. The pipeline integrates a series of 
                                                   bioinformatics softwares and scripts to find retrocopies. This website 
                                                   is mainly to display the results of retrocopies, which shows the 
                                                   overall statistical information, the structure diagram of retrocopies, 
                                                   the distribution of ka/ks and the heatmap of FPKM using the Shiny package 
                                                   in R.</h5>")
                                       ),
                                column(3)
                                       ),
                            
                            fluidRow(
                                
                                style = "height:50px;"),
                            
                            # PAGE BREAK
                            tags$hr(),
                            
                            # HOW TO START
                            fluidRow(
                                shiny::HTML("<br><br><center> <h1>How to get started.</h1> </center>
                                            <br>")
                                ),
                            
                            fluidRow(
                                column(3),
                                
                                column(2,
                                       div(class="panel panel-default", 
                                           div(class="panel-body",  width = "600px",
                                               align = "center",
                                               div(
                                                   tags$img(src = "one.svg", 
                                                            width = "50px", height = "50px")
                                               ),
                                               div(
                                                   h5(
                                                       "Upload the results of RetroScan, including \"final.out\", 
                                                       \"gene_trans_info.txt\", \"cds.fasta\" and \"all_samples.counts.txt\". 
                                                       And \"phenotype.txt\" is provided by users. "
                                                   )
                                               )
                                           )
                                       )
                                ),
                                column(2,
                                       div(class="panel panel-default",
                                           div(class="panel-body",  width = "600px", 
                                               align = "center",
                                               div(
                                                   tags$img(src = "two.svg", 
                                                            width = "50px", height = "50px")
                                               ),
                                               div(
                                                   h5(
                                                       "Then this website will count the results of each part and draw 
                                                       statistics, structure and other pictures."
                                                   )
                                               )
                                           )
                                       )
                                ),
                                column(2,
                                       div(class="panel panel-default",
                                           div(class="panel-body",  width = "600px", 
                                               align = "center",
                                               div(
                                                   tags$img(src = "three.svg", 
                                                            width = "50px", height = "50px")),
                                               div(
                                                   h5(
                                                       "You can filter the information and change the parameters 
                                                       reasonably and download the pictures you need."
                                                   )
                                               )
                                           )
                                       )
                                ),
                                column(3)
                                
                            ),
                            

                            fluidRow(
                                
                                style = "height:50px;")
                            
                            ), # Closes the first tabPanel called "Home"
                   
                   tabPanel("RETROSCAN", value = "careerPF",
                            tabsetPanel(
                                type = "tabs",
                                #input
                                tabPanel("input",style = "height:600px;",
                                         shiny::HTML("<br><br>"),
                                         fluidRow(
                                             column(2),
                                             column(9,
                                                    tags$li(shiny::HTML("<h5 style='line-height:30px'>These three files are necessary and can be found in 
                                                                        the retroscan_results folder of RetroScan results.<br><br></h5>"))
                                             ),
                                             column(1)
                                         ),
                                         fluidRow(
                                             column(2),
                                         column(3,
                                                fileInput('file1',
                                                          shiny::HTML("<h4 style='line-height:30px'>Final out File<br>(final.out)</h4>"),
                                                          accept = c('text/csv','text/comma-separated-values,text/plain','.csv','.out','.txt'))
                                         ),
                                         column(3,
                                                fileInput('file2',
                                                          shiny::HTML("<h4 style='line-height:30px'>Gene annotation File<br>(gene_trans_info.txt)</h4>"),
                                                          accept = c('text/csv','text/comma-separated-values,text/plain','.csv','.out','.txt'))
                                                ),
                                         column(3,
                                                fileInput('cds_file',shiny::HTML("<h4 style='line-height:30px'>CDS fasta File<br>(cds.fasta)</h4>"))
                                                ),
                                         column(1)
                                         ),
                                         shiny::HTML("<br><br>"),
                                         fluidRow(
                                             column(2),
                                             column(9,
                                                    tags$li(shiny::HTML("<h5 style='line-height:30px'>If you provide RNASeq data and get the FPKM values, 
                                                                        please input the FPKM (in the the retroscan_results folder of RetroScan results) and 
                                                                        phenotype file (provided by the user and need to correspond to RNASeq 
                                                                        data).<br><br></h5>"))
                                             ),
                                             column(1)
                                         ),
                                         fluidRow(
                                             column(2),
                                             column(3,
                                                    fileInput('fpkm_file',
                                                              shiny::HTML("<h4 style='line-height:30px'>FPKM File<br>(all_samples.counts.txt)</h4>"),
                                                              accept = c('text/csv','text/comma-separated-values,text/plain','.csv','.out','.txt'))
                                             ),
                                             column(3,
                                                    fileInput('phenotype_file',
                                                              shiny::HTML("<h4 style='line-height:30px'>Phenotype File<br>(Corresponds to RNASeq)</h4>"),
                                                              accept = c('text/csv','text/comma-separated-values,text/plain','.csv','.out','.txt'))
                                             ),
                                             column(4)
                                         )
                                ),

                                #summary
                                tabPanel("summary",
                                        fluidRow( 
                                            shiny::HTML("<center> <h2>Summary</h2> </center>"),
                                            column(3,
                                                   wellPanel(
                                        h4("Brief Information"),
                                        h6(tableOutput("summary_info")),
                                        style = "background: #EEEEEE"
                                        )),
                                        column(9,
                                             
                                            tabsetPanel(
                                                type = "tabs",
                                                tabPanel("List",
                                                         br(),
                                                         DTOutput("dynamic"),
                                                         br(),
                                                         downloadButton('Retrocopy_result','Download')),
                                                tabPanel("Chr-distribution",
                                                         box(width = 10,plotOutput("retro_parent_chr",width = "100%",height = "600px"),
                                                             p("The figure focuses on chromosomes/scaffolds larger than 1000kbp.")),
                                                         box(width = 2,
                                                             selectInput("chr_color", label = "Select fill colors", 
                                                                         choices = list("Accent","Set1","Set2","Set3","Pastel2","Pastel1","Paired","Dark2"), 
                                                                         selected = 'Accent'),
                                                             numericInput("chr_title", label = "Title size", value = 12),
                                                             numericInput("chr_xy", label = "X&Y title size", value = 10),
                                                             numericInput("chr_xy_text", label = "X&Y text size", value = 10),
                                                             numericInput("chr_x_angle", label = "X text angle", value = 60),
                                                             selectInput("chr_legend",label="Legend position",choices =list("none","left","right","bottom","top"),selected = 'bottom'),
                                                             br(),
                                                             radioButtons('chrplot', 'Output format', choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg','SVG'='svg'),inline = T),
                                                             downloadButton('chrplotdown',label = 'Download plot'))),
                                                tabPanel("Retro-parental-gene",
                                                         box(width = 10,plotOutput("retro_parent_num",width = "100%",height = "600px")),
                                                         box(width = 2,
                                                             numericInput("num_title", label = "Title size", value = 12),
                                                             numericInput("num_xy", label = "X&Y title size", value = 10),
                                                             numericInput("num_text", label = "X&Y text size", value = 10),
                                                             colourInput("num_color", "Select colour", "lightblue"),
                                                             radioButtons('numplot', 'Output format', choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg','SVG'='svg'),inline = T),
                                                             downloadButton('numplotdown',label = 'Download plot'))),
                                                tabPanel("Retro-length",
                                                         box(width = 10,plotOutput("retro_len",width = "100%",height = "600px")),
                                                         box(width = 2,
                                                             numericInput("len_title", label = "Title size", value = 12),
                                                             numericInput("len_xy", label = "X&Y title size", value = 10),
                                                             numericInput("len_text", label = "X&Y text size", value = 10),
                                                             numericInput("len_bin", label = "Bin size", value = 100),
                                                             colourInput("len_color", "Select colour", "lightblue"),
                                                             radioButtons('lenplot', 'Output format', choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg','SVG'='svg'),inline = T),
                                                             downloadButton('lenplotdown',label = 'Download plot'))),
                                                tabPanel("Identity",
                                                         box(width = 10,plotOutput("identity_plot",width = "100%",height = "600px")),
                                                         box(width = 2,
                                                             selectInput("iden_color", label = "Select fill colors", 
                                                                         choices = list("Accent","Set1","Set2","Set3","Pastel2","Pastel1","Paired","Dark2"), 
                                                                         selected = 'Accent'),
                                                             numericInput("iden_title", label = "Title size", value = 12),
                                                             numericInput("iden_text", label = "Text size", value = 6),
                                                             selectInput("iden_legend",label="Legend position",choices =list("none","left","right","bottom","top"),selected = 'right'),
                                                             numericInput("iden_legend_size", label = "Legend size", value = 10),
                                                             radioButtons('idenplot', 'Output format', choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg','SVG'='svg'),inline = T),
                                                             downloadButton('idenplotdown',label = 'Download'))),
                                                tabPanel("Coverage",
                                                         box(width = 10,plotOutput("coverage_plot",width = "100%",height = "600px")),
                                                         box(width = 2,
                                                             selectInput("cov_color", label = "Select fill colors", 
                                                                         choices = list("Accent","Set1","Set2","Set3","Pastel2","Pastel1","Paired","Dark2"), 
                                                                         selected = 'Accent'),
                                                             numericInput("cov_title", label = "Title size", value = 12),
                                                             numericInput("cov_text", label = "Text size", value = 6),
                                                             selectInput("cov_legend",label="Legend position",choices =list("none","left","right","bottom","top"),selected = 'right'),
                                                             numericInput("cov_legend_size", label = "Legend size", value = 10),
                                                             radioButtons('covplot', 'Output format', choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg','SVG'='svg'),inline = T),
                                                             downloadButton('covplotdown',label = 'Download'))),
                                                tabPanel("Category",
                                                         box(width = 10,plotOutput("description_plot",width = "100%",height = "600px")),
                                                         box(width = 2,
                                                             selectInput("des_color", label = "Select fill colors", 
                                                                         choices = list("Accent","Set1","Set2","Set3","Pastel2","Pastel1","Paired","Dark2"), 
                                                                         selected = 'Accent'),
                                                             numericInput("des_title", label = "Title size", value = 12),
                                                             numericInput("des_text", label = "Text size", value = 6),
                                                             selectInput("des_legend",label="Legend position",choices =list("none","left","right","bottom","top"),selected = 'right'),
                                                             numericInput("des_legend_size", label = "Legend size", value = 10),
                                                             radioButtons('desplot', 'Output format', choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg','SVG'='svg'),inline = T),
                                                             downloadButton('desplotdown',label = 'Download'))),
                                                tabPanel("UpSetR",
                                                         box(width = 10,plotOutput("venn",width = "100%",height = "600px")),
                                                         box(width = 2,
                                                          #   selectInput("venn_method", label = "Select venn method",choices = list("upSetR","VennDiagram"),selected = 'upSetR'),
                                                             numericInput("iden_cutoff",label="Identity Threshold",min=50,max=100,value=90),
                                                             numericInput("cov_cutoff",label="Coverage Threshold",min=50,max=100,value=90),
                                                             numericInput("intron_cutoff",label="Lost intron Threshold",min=2,max=100,value=3),
                                                             numericInput("kaks_cutoff",label="Ka/Ks Threshold",min=0,max=1,step=0.01,value=0.1),
                                                             colourInput("venn_matrix_color", "Select matrix colour", "gray23"),
                                                             colourInput("venn_main_bar_color", "Select main bar colour", "gray23"),
                                                             colourInput("venn_sets_bar_color", "Select sets bar colour", "gray23"),
                                                             colourInput("venn_shade_color", "Select shade colour", "gray88"),
                                                             radioButtons('vennplot', 'Output format', choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg','SVG'='svg'),inline = T),
                                                             downloadButton('vennplotdown',label = 'Download')))
                                            ))
                                        ) 
                                ),
                                #Retrocopy Sample
                                tabPanel("retrocopy", 
                                         fluidRow(
                                         shiny::HTML("<center> <h2>Retrocopy Information</h2> </center>"),
                                         column(3,
                                                wellPanel(
                                              #  h4("Brief Information"),
                                                textInput("sample1",label=h4("Retrocopy ID search"),value=("Retrocopy1")),
                                                tabsetPanel(
                                            type = "tabs",
                                            tabPanel("Retrocopy",h6(tableOutput("retro_info"))),
                                            tabPanel("Parental gene",h6(tableOutput("parent_info"))),
                                            tabPanel("Host gene",h6(tableOutput("host_info")))),
                                            style = "background: #EEEEEE"
                                                )),
                                         column(9,
                                        tabsetPanel(
                                            type = "tabs",
                                            #utr_color,exon_color,intron_color,frame_color,text_size,title_size
                                            tabPanel("Structure",
                                                     box(width = 10,plotOutput("retro_plot",width = "100%",height = "600px")),
                                                     box(width = 2,
                                                         colourInput("exon_color", "Select exon color", "#377EB8"),
                                                         colourInput("utr_color", "Select UTR color", "#4DAF4A"),
                                                         colourInput("intron_color", "Select intron color", "#999999"),
                                                         colourInput("frame_color", "Select retrocopy color", "#E41A1C"),
                                                         numericInput("title_size", label = "Title size", value = 2),
                                                         numericInput("text_size", label = "Text size", value = 1),
                                                         radioButtons('retroplot', 'Output format', choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg','SVG'='svg'),inline=T),
                                                         downloadButton('retroplotdown',label = 'Download')
                                                     )),
                                            tabPanel("Sequence",
                                                     br(),
                                                tabsetPanel(
                                                    
                                                    tabPanel("Retrocopy DNA",verbatimTextOutput("retroseq", placeholder = TRUE)),
                                                    tabPanel("Retrocopy PEP",verbatimTextOutput("retropep", placeholder = TRUE))
                                                 #   tabPanel("Parental Gene CDS",verbatimTextOutput("retropep", placeholder = TRUE)),
                                                 #   tabPanel("Parental Gene PEP",verbatimTextOutput("retropep", placeholder = TRUE))
                                                )
                                                ),
                                            tabPanel("Alignment",
                                                     verbatimTextOutput("alignment")),
                                                     #box(width = 10,plotOutput("alignment",width = "100%",height = "1000px")),
                                                   #  box(width = 2,
                                                     #    selectInput("aln_color", label = "Select alignment colors", 
                                                     #                choices = list("Zappo_AA","Taylor_AA","Shapely_AA","Chemistry_AA","Clustal"), 
                                                      #               selected = 'Chemistry_AA'))),
                                            tabPanel("Expression",
                                                     box(width = 10,plotOutput("retro_exp_plot",width = "100%",height = "600px")),
                                                     box(width = 2,
                                                         colourInput("exp_parentcolor", "Select parental gene color", "#377EB8"),
                                                         colourInput("exp_retrocolor", "Select retrocopy color", "#4DAF4A"),
                                                         numericInput("exp_titlesize", label = "Title size", value = 15),
                                                         numericInput("exp_textsize", label = "Text size", value = 10),
                                                         numericInput("exp_linesize", label = "Line size", value = 2),
                                                         numericInput("exp_xangle", label = "X text angle", value = 60),
                                                         radioButtons('retroexpplot', 'Output format', choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg','SVG'='svg'),inline=T),
                                                         downloadButton('retroexpdown',label = 'Download')
                                                     ))
                                        ))
                                )),
                                #kaks
                                tabPanel("kaks",
                                         fluidRow(
                                             shiny::HTML("<center> <h2>KaKs Information</h2> </center>"),
                                             column(3,
                                                    wellPanel(
                                                        h4("Brief Information"),
                                                        h6(tableOutput("kaks_info")),
                                                        style = "background: #EEEEEE"
                                                    )),
                                             column(9,
                                        tabsetPanel(
                                            type = "tabs",
                                            tabPanel("KaKs Table",
                                                     br(),
                                                     DTOutput("kaks_result"),
                                                     br(),
                                                     downloadButton('kaks_result_download','Download')
                                            ),
                                            tabPanel("Ks",
                                                     box(width = 10,plotOutput("ks_len",width = "100%",height = "600px")),
                                                     box(width = 2,
                                                         numericInput("ka_title", label = "Title size", value = 12),
                                                         numericInput("ka_xy", label = "X&Y title size", value = 10),
                                                         numericInput("ka_text", label = "X&Y text size", value = 10),
                                                         numericInput("ka_bin", label = "Bin size", value = 0.01,step=0.01),
                                                         colourInput("ka_col", "Select color", "lightblue"),
                                                         radioButtons('kaplot', 'Output format', choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg','SVG'='svg'),inline = T),
                                                         downloadButton('ksplotdown',label = 'Download plot'))),
                                            tabPanel("KaKs",
                                                     box(width = 10,plotOutput("kaks_len",width = "100%",height = "600px")),
                                                     box(width = 2,
                                                         numericInput("kaks_title", label = "Title size", value = 12),
                                                         numericInput("kaks_xy", label = "X&Y title size", value = 10),
                                                         numericInput("kaks_text", label = "X&Y text size", value = 10),
                                                         numericInput("kaks_bin", label = "Bin size", value = 0.01,step=0.01),
                                                         colourInput("kaks_col", "Select color", "lightblue"),
                                                         radioButtons('kaksplot', 'Output format', choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg','SVG'='svg'),inline = T),
                                                         downloadButton('kaksplotdown',label = 'Download plot'))),
                                            tabPanel("Ks Distribution",
                                                     box(width = 10,plotOutput("ks_type_plot",width = "100%",height = "600px")),
                                                     box(width = 2,
                                                         selectInput("ka_type_color", label = "Select fill colors", 
                                                                     choices = list("Accent","Set1","Set2","Set3","Pastel2","Pastel1","Paired","Dark2"), 
                                                                     selected = 'Accent'),
                                                         numericInput("ka_type_title", label = "Title size", value = 12),
                                                         numericInput("ka_type_xy", label = "X&Y title size", value = 10),
                                                         numericInput("ka_type_xy_text", label = "X&Y text size", value = 10),
                                                         numericInput("ka_type_bin", label = "Bin size", value = 0.05,step=0.01),
                                                         selectInput("ka_type_legend",label="Legend position",choices =list("none","left","right","bottom","top"),selected = 'bottom'),
                                                         br(),
                                                         radioButtons('ka_typeplot', 'Output format', choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg','SVG'='svg'),inline = T),
                                                         downloadButton('ks_typeplotdown',label = 'Download plot'))),
                                            tabPanel("KaKs Distribution",
                                                     box(width = 10,plotOutput("kaks_type_plot",width = "100%",height = "600px")),
                                                     box(width = 2,
                                                         selectInput("kaks_type_color", label = "Select fill colors", 
                                                                     choices = list("Accent","Set1","Set2","Set3","Pastel2","Pastel1","Paired","Dark2"), 
                                                                     selected = 'Accent'),
                                                         numericInput("kaks_type_title", label = "Title size", value = 12),
                                                         numericInput("kaks_type_xy", label = "X&Y title size", value = 10),
                                                         numericInput("kaks_type_xy_text", label = "X&Y text size", value = 10),
                                                         numericInput("kaks_type_bin", label = "Bin size", value = 0.05,step=0.01),
                                                         selectInput("kaks_type_legend",label="Legend position",choices =list("none","left","right","bottom","top"),selected = 'bottom'),
                                                         br(),
                                                         radioButtons('kaks_typeplot', 'Output format', choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg','SVG'='svg'),inline = T),
                                                         downloadButton('kaks_typeplotdown',label = 'Download plot')))
                                        ))
                                )),
                                #expression
                                tabPanel("expression",
                                         fluidRow(
                                             shiny::HTML("<center> <h2>Expression Information</h2> </center>"),
                                             column(3,
                                                    wellPanel(
                                                        h4("Brief Information"),
                                                        h6(tableOutput("expr_info")),
                                                        style = "background: #EEEEEE"
                                                    )),
                                             column(9,
                                        tabsetPanel(
                                            type = "tabs",
                                            tabPanel("Expression table",
                                                     br(),
                                                     DTOutput("expr_result"),
                                                     br(),
                                                     downloadButton('expr_result_download','Download')
                                            ),
                                            tabPanel("Expression number",
                                                     box(width = 10,plotOutput("exprplot",width = "100%",height = "600px")),
                                                     box(width = 2,
                                                         colourInput("expr_parentcolor", "Select parental gene color", "#377EB8"),
                                                         colourInput("expr_retrocolor", "Select retrocopy color", "#4DAF4A"),
                                                         numericInput("expr_titlesize", label = "Title size", value = 12),
                                                         numericInput("expr_text_size", label = "Text size", value = 10),
                                                         numericInput("expr_xangle", label = "X text angle", value = 60),
                                                         selectInput("expr_legend",label="Legend position",choices =list("none","left","right","bottom","top"),selected = 'bottom'),
                                                         radioButtons('expr_typeplot', 'Output format', choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg','SVG'='svg'),inline = T),
                                                         downloadButton('expr_plotdown',label = 'Download plot'))),
                                            tabPanel("heatmap",
                                                     br(),
                                                     plotOutput("exprheatmap"),
                                                     br(),
                                                     radioButtons('expr_heat', 'Output format', choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg','SVG'='svg'),inline = T),
                                                     downloadButton('expr_heatdown',label = 'Download plot')
                                            )
                                        ))
                                ))
                            )# Closes tabsetPanel
                   ),  # Closes the second tabPanel called "RetroScan"
                   
                   tabPanel("ABOUT", value = "about",
                            
                            fluidRow(
                                shiny::HTML("<br><br><center> 
                                            <h1>About RetroScan</h1> 
                                            </center>
                                            <br>
                                            <br>"),
                                style = "height:200px;"),
 
                            fluidRow(
                                div(align = "center",
                                    tags$span(h3("The pipeline of RetroScan"), 
                                              style = "font-weight:bold"
                                    ))
                            ),
                            
                            fluidRow(
                                column(2),
                                column(8,
                                       tags$p(
                                           shiny::HTML("<h6 style='line-height:25px'>RetroScan is an easy-to-use tool for 
                                                       retrocopy identification, which integrates a series of bioinformatics 
                                                       tools (LAST, BEDtools, ClustalW2, KaKs_Calculator, HISAT2, StringTie, 
                                                       samtools and Shiny ) and scripts. It scans retrocopies based on the 
                                                       alignments between protein-coding genes and the whole genome sequences. 
                                                       This tool can also analyze heterosense substitution and synonymous 
                                                       substitution, compare the gene structure between parental genes and 
                                                       retrocopies, and calculate the expression values of them. Moreover, 
                                                       RetroScan has a user-friendly visualization interface, which showing 
                                                       the overall statistical information, the structure diagram of 
                                                       retrocopies, the distribution of ka/ks and the heatmap of FPKM using 
                                                       the Shiny package in R.</h6>")),
                                       tags$ul(
                                           tags$li(
                                               h5("Retrocopy identification"),
                                               shiny::HTML("<h6 style='line-height:25px'>RetroScan requires at least two input 
                                                           files : genome sequences file (fasta format) and corresponding 
                                                           annotation file (gff format), then it can identify the detailed 
                                                           information of retrocopies and parental genes on the genome. If users 
                                                           want to obtain the expression values of retrocopies, they need to 
                                                           submit additional RNA-Seq data. According to genome sequences and gff 
                                                           file, RetroScan first generate the peptide sequences which are used 
                                                           as queries in similarity searches against the complete genome 
                                                           sequences using LAST to find out candidate hits. To avoid duplicate 
                                                           results, the longest transcripts of each gene for alignment are 
                                                           remained for next step. Multi-exon proteins were selected for 
                                                           subsequent analysis because the parental genes must loss introns 
                                                           not less than two. According to the alignment results in the previous 
                                                           step, users can set the sequence identity, coverage and alignment 
                                                           length parameters to consider the specific conditions of the species. 
                                                           Multiple alignment hits to the same genomic locus were clustered 
                                                           using BEDTools. When the distance between the hits is less than a 
                                                           certain length that are not likely separated by introns (gap defaulted 
                                                           as 40bp by RetroScan), adjacent homology hits are merged using BEDTools. 
                                                           Next, the merged sequences are aligned back to multi-exon proteins 
                                                           using LAST and keep best hits as putative parental genes. At last, 
                                                           estimating the number of lost introns is applied to get reliable 
                                                           results according the alignment output. RetroScan only conserves 
                                                           parental genes (excluding the first and last 10 amino acids) which 
                                                           must span at least two introns and single-exon retrocopies. If there 
                                                           are a large numbers of retrocopies (max_parent_retrocopy_pair 
                                                           defaulted as 10) assigning to a single parental gene, filter these 
                                                           retrocopies in order to minimize the number of false positive findings. 
                                                           In addition, the retrocopies with either premature stop codons or frameshift 
                                                           mutations are defined as retropseudogenes otherwise as intact retrocopies. 
                                                           If one intact retrocopy could recruit novel regulatory elements or new 
                                                           protein-coding exons and evolve into a functional retrogene, it can be 
                                                           defined as a chimerical retrogene.</h6>")),
                                           tags$li(
                                               h5("Ka/Ks analysis"),
                                               shiny::HTML("<h6 style='line-height:25px'>The age distribution of the retrocopies is 
                                                           determined by calculating Ka (nonsynonymous substitutions), Ks (synonymous 
                                                           substitutions) and Ka/Ks radios between each retrocopy and its parental gene. 
                                                           The coding sequences (CDS) information of the retrocopies and their parental 
                                                           genes based on the annotation file are extracted to the Ka/Ks calculation. 
                                                           Then, RetroScan executes multiple alignment between their protein sequences 
                                                           using ClustalW2. Finally, calculate the Ka, Ks and Ka/Ks radios using 
                                                           KaKs_calculator_2.0.</h6>")),
                                           tags$li(
                                               h5("Retrocopy expression analysis"),
                                               shiny::HTML("<h6 style='line-height:25px'>In order to estimate expression values of 
                                                           retrogenes, RetroScan uses HISAT2, samtools and stringtie to analyze the 
                                                           RNA-Seq data based on the position information of retrocopies and their 
                                                           parental genes. After the reads mapped the corresponding annotated sequences 
                                                           using HISAT2, RetroScan converts sam files into bam files and sorts them 
                                                           using samtools. Finally, stringtie calculates FPKM (Fragments Per Kilobase 
                                                           per Million), which is helpful for analyzing differential expression. All 
                                                           programs are setting with the default parameters.</h6>")),
                                           tags$li(
                                               h5("Visualization"), 
                                               shiny::HTML("<h6 style='line-height:25px'>We develop a visual interface that can clearly 
                                                           display the structure of retrocopies, the distribution of ka/ks, expression 
                                                           heatmap, sequence alignment and statistical figures using the Shiny package 
                                                           in R . The interface layout is mainly divided into five parts: (1) Home; (2) 
                                                           Summary; (3) Retrocopy; (4) KaKs; (5) Expression. Users can filter data based 
                                                           on any column of tables in each page and directly search for keywords in the 
                                                           search box above the tables. All pictures colors and texts size can be 
                                                           adjusted according to users’ needs. All information tables and figures can be 
                                                           downloaded by clicking download tabs.</h6>"))
                                       )
                                       
                                ),
                                column(2)
                            ),
                            
                            
                            fluidRow(
                                div(align = "center",
                                    tags$span(h3("How to cite"), 
                                              style = "font-weight:bold"
                                    ))
                            ),
                            
                            fluidRow(
                                column(2),
                                column(8,
                                       tags$p(
                                           shiny::HTML("<h6 style='line-height:25px'>RetroScan: an easy-to-use systematic pipeline 
                                                       for retrocopy annotation and visualization.</h6>"))
                                ),
                                column(2)
                            ),
                      
                        
                            fluidRow(
                                div(align = "center",
                                    tags$span(h3("Contact"), 
                                              style = "font-weight:bold"
                                    ))
                            ),
                            
                            fluidRow(
                                column(2),
                                column(8,
                                       tags$p(
                                           shiny::HTML("<h6 style='line-height:25px'>wzy19951225@email.swu.edu</h6>"))
                                ),
                                column(2)
                            ),
                            
                            
                         
                            fluidRow(style = "height:150px;")
                                )  # Closes About tab
                   
                            )
        
                   )