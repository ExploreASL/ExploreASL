#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


## USERS PARAMETER  ##
#Load Libraries
library(shiny)
library(ggplot2)
library(readr)
library(ggplot2)
library(gridExtra)
library(RNifti)
library(oro.nifti)
library(shinyjs)
library(wesanderson)
library(plyr)
library(openxlsx)
library(plotly)
library(dplyr)

######### Functions ######### 
#Find outliers
isout <-  function(df){
  
  nums <- unlist(lapply(df, is.numeric))
  dfofint <- df[, nums]
  outliers = data.frame(matrix(NaN, nrow = nrow(dfofint), ncol = ncol(dfofint)))
  colnames(outliers)<- colnames(df[, nums])
  for (iobs in c(1:nrow(dfofint))) {
    for (ivar in c(1:ncol(dfofint))){
      if(is.na(dfofint[iobs, ivar])) {
        outliers[iobs, ivar] = NaN
      } else if (dfofint[iobs, ivar] >= (mean(dfofint[,ivar], na.rm = TRUE) + 1.96*sd(dfofint[,ivar], na.rm = TRUE))){
        outliers[iobs, ivar] = 1
      } else if (dfofint[iobs, ivar] <= (mean(dfofint[,ivar], na.rm = TRUE) - 1.96*sd(dfofint[,ivar], na.rm = TRUE))){
        outliers[iobs, ivar] = 1
      } else {
        outliers[iobs, ivar] = 0
      }
    }
  }
  return(outliers)
}

#Compute Z scores
isz <-  function(df){
  
  nums <- unlist(lapply(df, is.numeric))
  dfofint <- df[, nums]
  zed = data.frame(matrix(NaN, nrow = nrow(dfofint), ncol = ncol(dfofint)))
  colnames(zed)<- colnames(df[, nums])
  for (iobs in c(1:nrow(dfofint))) {
    for (ivar in c(1:ncol(dfofint))){
      if(is.na(dfofint[iobs, ivar])) {
        zed[iobs, ivar] = NaN
      } else{
        zed[iobs, ivar] = (dfofint[iobs, ivar]- mean(dfofint[, ivar]))/sd(dfofint[, ivar])
      
      }
    }
  }
  return(zed)
}

##############

#Read the Configuration File and set variables
configuration = openxlsx::read.xlsx("../ConfigFile.xlsx", 1)

analysisdir = configuration[which(configuration$Visualization.Module == 'AnalysisDir'), 'VIS_properties'] 
              
T1image = configuration[which(configuration$Visualization.Module == 'StructIm'), 'VIS_properties']  
FunctionalFold = configuration[which(configuration$Visualization.Module == 'FuncFold'), 'VIS_properties'] 
FunctionalImage = configuration[which(configuration$Visualization.Module == 'FuncIm'), 'VIS_properties'] 
DiffFolder = configuration[which(configuration$Visualization.Module == 'DiffFold'), 'VIS_properties'] 
DiffImage = configuration[which(configuration$Visualization.Module == 'DiffIm'), 'VIS_properties'] 



# Set directory (tmp) and read the output of xQC_Master
if (file.exists(file.path("dataframes", "QCed_data.csv"))) {
  qc_data_all = read.csv(file.path("dataframes", "QCed_data.csv"))
  
  
  print("One previous file has been found, starting from there. Delete the file to start over")
} else {
  qc_data_all = read.csv(file.path("dataframes", "QC.csv"))
  
  #Subset based on variables you want to visualize (from configuration file)
  includeparms = openxlsx::read.xlsx("../ConfigFile.xlsx", 3)
  includeparms = includeparms[,c(1:3,5)]
  includeparms  = includeparms[c(which(includeparms$Visualize == 1)),]
  parmstoinclude = paste(includeparms$Scantype, includeparms$Domain, includeparms$Parameter, sep = '_')
  parmstoinclude = c("Subject", "Site", parmstoinclude)
  qc_data_all <- qc_data_all[, names(qc_data_all) %in% parmstoinclude]
  
  qc_data_all$Structural = rep("Good",nrow(qc_data_all)) 
  qc_data_all$Functional = rep("Good",nrow(qc_data_all)) 
  qc_data_all$Diffusion = rep("Good",nrow(qc_data_all)) 
  qc_data_all$ASL = rep("Good",nrow(qc_data_all)) 
  
  #this won't be needed later
  qc_data_all[qc_data_all == 0] <- NaN
  
  # numbers of outliers and zscores
  qc_data_all$n_struct_outliers <- 0
  qc_data_all$z_score_struct <- 0
  
  # If it's the first time we need to flag scans 
  # Flag using defined function
  for (isit in unique(qc_data_all$Site)){
    sitdf <- qc_data_all[which(as.character(qc_data_all$Site) ==isit),]
    sitoutliers <- isout(sitdf)
    sitezed <- isz(sitdf)
    for (isub in c(1:nrow(sitdf))){
      
      # structural 
      xS = sum(sitoutliers[isub,startsWith(colnames(sitoutliers), 'Structural')], na.rm = TRUE)
      if (between(xS, 2, 3)){ # 2 or 3 outliers
        sitdf$Structural[isub] <- "Moderate"
      } else if (xS > 3){
        sitdf$Structural[isub] <- "Poor"
      }
      
      # save the total number of structural outliers for that subject
      sitdf$n_struct_outliers[isub] <- xS 
      sitdf$z_score_struct[isub] <- sum(abs(sitezed[isub, startsWith(colnames(sitoutliers), 'Structural')]))
      
      # functional 
      xF = sum(sitoutliers[isub,startsWith(colnames(sitoutliers), 'Functional')], na.rm = TRUE)
      if (between(xF, 2, 3)){ # 2 or 3 outliers
        sitdf$Functional[isub] <- "Moderate"
      } else if (xF > 3){
        sitdf$Functional[isub] <- "Poor"
      }
      
      # DTI
      xD = sum(sitoutliers[isub,startsWith(colnames(sitoutliers), 'Diffusion')], na.rm = TRUE)
      if (between(xD, 2, 3)){ # 1 or 2 outliers
        sitdf$Diffusion[isub] <- "Moderate"
      } else if (xD > 3){
        sitdf$Diffusion[isub] <- "Poor"
      }
      
    }
    
    # Put the results in the main dataframe 
    qc_data_all[which(as.character(qc_data_all$Site) ==isit),'Structural'] <- sitdf$Structural
    qc_data_all[which(as.character(qc_data_all$Site) ==isit),'n_struct_outliers'] <- sitdf$n_struct_outliers
    qc_data_all[which(as.character(qc_data_all$Site) ==isit),'z_score_struct'] <- sitdf$z_score_struct
    
    qc_data_all[which(as.character(qc_data_all$Site) ==isit),'Functional'] <- sitdf$Functional
    qc_data_all[which(as.character(qc_data_all$Site) ==isit),'Diffusion'] <- sitdf$Diffusion
    
  }
}


#Manage Site and patient name, This is temporary since is specific for EPAD study
qc_data_all$Site <- as.factor(qc_data_all$Site)


#This should already be present in the dataframe, will be left out later
qc_data_all$patient = qc_data_all$Subject




####################

# Create variable with number of scan per site to 
for (i in 1:nrow(qc_data_all)) {
  qc_data_all$numberofscanspersite[i]=nrow(qc_data_all[which(qc_data_all$Site == qc_data_all$Site[i]),])
}


# Define UI for application that draws a histogram
ui <- fluidPage(
  useShinyjs(),
  
  wellPanel(
    # Application title
    titlePanel("Explore QC visualization Tool"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
      
      
      sidebarPanel(
        
        selectInput("modality", 
                    "Choose The MRI modality",
                    choices = c("Structural", "Functional", "Diffusion", "ASL"), 
                    selected = "Structural"),
        
        sliderInput("nscans",
                    "Exclude Sites with fewer scans than:",
                    min = 1,
                    max = 50,
                    value = 0),
        
        uiOutput("inputUI")
        
      ),
      
      
      # Show a plot of the generated distribution
      mainPanel(
        plotlyOutput("violinplot"),
        
        plotlyOutput("scatterplot"),#, click = "wt"), 
        
        textOutput("Text1"),
        
        #imageOutput("myImage"),
        fluidRow( 
          column(5, 
                 sliderInput("sliceX", 
                             label = "Axial",
                             min = 1,
                             max = 250,
                             value = 170),
                 
                 plotOutput("axial"),
                 includeScript("https://cdnjs.cloudflare.com/ajax/libs/jquery-mousewheel/3.1.13/jquery.mousewheel.min.js")
          ), 
          column(5, offset = 1,
                 sliderInput("sliceY", 
                             label = "Sagittal",
                             min = 1,
                             max = 170,
                             value = 170),
                 
                 plotOutput("sagittal"),
                 includeScript("https://cdnjs.cloudflare.com/ajax/libs/jquery-mousewheel/3.1.13/jquery.mousewheel.min.js")
          )
          
          
        ),
        
        fluidRow(
          column(5, 
                 sliderInput("sliceZ", 
                             label = "Coronal",
                             min = 1,
                             max = 250,
                             value = 250),
                 
                 plotOutput("coronal"),
                 includeScript("https://cdnjs.cloudflare.com/ajax/libs/jquery-mousewheel/3.1.13/jquery.mousewheel.min.js")
          ), 
          column(4, offset = 1,
                 
                 uiOutput("mygui") ,
                 
                 "\n  \n  \n  Do you want to save the results of the visual QC? \n The results will be stored in a dataframen and you can continue editing later",
                 actionButton("save", "Save QCed dataframe")
                 
                 
          )
          
          
        )
      )
      
    )
    ,style = "overflow-y:scroll; max-height: 100%"
  )
)
# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  ## ADMINISTRATION  ##
  # select columns based on inputs
  qc_data <- reactive({
    data.frame(qc_data_all[,which(startsWith(colnames(qc_data_all), input$modality))], 
               qc_data_all$Site, 
               qc_data_all$patient,
               qc_data_all$numberofscanspersite,
               qc_data_all[,input$modality]
    )
  })
  
  #qc_data_TH -> qc_data but with just sites with more scans than threshold (input$nscans)
  qc_dataTH <- reactive({qc_data()[which(qc_data()$qc_data_all.numberofscanspersite>=input$nscans),] })
  
  # Select inputs based on modality (just parameters of one modality)
  
  output$inputUI<- renderUI({
    tagList(
      
      selectInput("QCpar",
                  label = "Select QC parameter:", 
                  choices = colnames(qc_data()),
                  selected = "Structural_Noise_SNR_GM_Ratio" ),
      
      selectInput("site",
                  label = "Select Site",
                  choices = unique(qc_dataTH()[,which(colnames(qc_data())== "qc_data_all.Site")]),
                  selected =  "040")
    )
   
  })
  
  
  # SiteData --> subset of thresholded dataframe includin only the selected site
  SiteData = reactive({qc_dataTH()[which(qc_dataTH()$qc_data_all.Site==input$site),]})
  
  
  
  
  ##   VISUALIZATION  ##
  pal <- c("green", "yellow" ,"red")
  #Violin Plots for Between sites distributions
  output$violinplot <- renderPlotly({
    violin <- plot_ly(qc_dataTH(),
    x = ~qc_dataTH()$qc_data_all.Site,
    y = ~qc_dataTH()[,input$QCpar],
    split = ~qc_dataTH()$qc_data_all.Site,
    type = 'violin',
    box = list(
      visible = T
    ))
    violin %>% layout(title = 'Between-Site Distribution',
                      xaxis = list(title = 'Sites', zeroline = FALSE),
                      yaxis = list(title = input$QCpar, zeroline = FALSE ), showlegend = FALSE)
  })
  
  # Scatter Plot for Within site distributions 
  output$scatterplot <-  renderPlotly({
    # color is SiteData()$includecolumn
   scatter1 <- plot_ly(data = SiteData(), x = ~SiteData()$qc_data_all.patient, y = ~SiteData()[,input$QCpar], color = ~qc_data_all[which(qc_data_all$Site == input$site), input$modality], 
                       type = "scatter",
                       mode = "markers", 
                       colors = pal, 
                       marker = list(size =10, line = list(color = "black", width = 2) ))
   scatter1 %>% layout(title = paste('Within-Site Distribution Site', input$site),
                       xaxis = list(title = 'Subjects' ,
                                    zeroline = FALSE),
                       yaxis = list(title = input$QCpar, zeroline = FALSE ), showlegend = TRUE) 
                       #shapes=list(type='line', 
                                   #x0= 0, 
                                   #x1=length(qc_data_all$patient), 
                                   #y0 = mean(SiteData()[,input$QCpar]),
                                   #y1=mean(SiteData()[,input$QCpar]), 
                                   #line=list(dash='dot', width=1, color = "green")))
  })
  
  output$Text1 <- renderText({d <- event_data("plotly_click")$x
  })

  
  
  
  
  # Set image Path and Load Image based on PatientNum   
  
  Patientnum <- reactive({event_data("plotly_click")$x})
  
  ImagePath<- NULL
  im <- NULL
  
  makeReactiveBinding("ImagePath")
  makeReactiveBinding("im")
  observeEvent(Patientnum(), { 
      if (input$modality == 'Structural'){
        ImagePath <<- reactive({file.path(analysisdir, Patientnum(), T1image)})
        im <<- reactive({oro.nifti::readNIfTI(ImagePath(), reorient = FALSE)})
      } else if (input$modality == 'Functional'){
        ImagePath <<- reactive({file.path(analysisdir, Patientnum(), FunctionalFold, FunctionalImage)})
        im <<- reactive({oro.nifti::readNIfTI(ImagePath(), reorient = FALSE) })
      }
      else if (input$modality == 'Diffusion'){
        ImagePath <<- reactive({file.path(analysisdir, Patientnum(), DiffFolder, DiffImage)})
        im <<- reactive({oro.nifti::readNIfTI(ImagePath(), reorient = FALSE)})
      }
    })
  
  
  
  ## MRI IMAGE Visualization ##
  
  # Plot the three views and allow the scrolling, unfortunately the scrolling has some problems but this is the best we can do so far
  #Axial
  output$axial <- renderPlot({
    if (is.null(Patientnum)) return(NULL)
    oro.nifti::slice(im(), z=input$sliceX)
  })
  onevent("mousewheel", "axial", {
    updateSliderInput(session, "sliceX", value = (input$sliceX-5))
  })
  #Sagittal
  output$sagittal <-renderPlot({
    if (is.null(Patientnum)) return(NULL)
    oro.nifti::slice(im(), z=input$sliceY, plane = "sagittal")
    
  })
  onevent("mousewheel", "sagittal", {
    updateSliderInput(session, "sliceY", value = input$sliceY-5)
  })
  #Coronal
  output$coronal <-renderPlot({
    if (is.null(Patientnum)) return(NULL)
    oro.nifti::slice(im(), z=input$sliceZ, plane = "coronal")
    
  })
  onevent("mousewheel", "coronal", {
    updateSliderInput(session, "sliceZ", value = input$sliceZ-5)
  })
  
  
  
  
  
  
  
  ## DECISION PART ## here we set up the user interface to include exlude scans
  
  
  # This is a piece of code that respond to the clicking of CONFIRM botton and:
  observeEvent(input$update,{  
    # 1. Save the output (Pass/Fail) in the qc_data_all dataframe in the column called as the modality
    if (input$include == "Good") {
      qc_data_all[which(qc_data_all$patient == Patientnum()),input$modality] <<-"Good"
      #SiteData()[which(SiteData()$qc_data_all.patient == Patientnum()),input$modality]<<-"Passed"
      #SiteData()$include[which(SiteData()$patient == Patientnum()),]<-1
      
    }  else if (input$include == "Poor"){
      qc_data_all[which(qc_data_all$patient == Patientnum()),input$modality] <<- "Poor"
      #SiteData()[which(SiteData()$qc_data_all.patient == Patientnum()),input$modality]<<-"Passed"
      
      
    } else if (input$include == "Moderate"){
      qc_data_all[which(qc_data_all$patient == Patientnum()),input$modality] <<- "Moderate"
      #SiteData()[which(SiteData()$qc_data_all.patient == Patientnum()),input$modality]<<-"Passed"

    }
    

    # 2. Update Plot
    output$scatterplot <-  renderPlotly({
      scatter1 <- plot_ly(data = SiteData(), x = ~SiteData()$qc_data_all.patient, y = ~SiteData()[,input$QCpar], color = ~qc_data_all[which(qc_data_all$Site == input$site), input$modality], 
                          type = "scatter",
                          mode = "markers", 
                          colors = pal, 
                          marker = list(size =10, line = list(color = "black", width = 2) ))
      scatter1 %>% layout(title = paste('Within-Site Distribution Site', input$site),
                          xaxis = list(title = 'Subjects' ,
                                       zeroline = FALSE),
                          yaxis = list(title = input$QCpar, zeroline = FALSE ), showlegend = TRUE)
    })
    
    
    
  })
  
  
  # This is a gui to print the Fail/Pass choice and to select the option if the scan has already been seen, 
  #(if you already excluded the subject the gui wil appear with the Fail option selected) 
  output$mygui = renderUI({
    if (is.null(Patientnum)) return(NULL)
    if (qc_data_all[which(qc_data_all$patient == Patientnum()),input$modality] == "Poor"){
      tagList(
        (radioButtons("include", 
                      label = "How do you judge the Scan's quality?",
                      choices = c("Good","Moderate", "Poor"),
                      selected = "Poor",
                      inline = TRUE
        )),
        (actionButton("update", "Confirm")))
    } else if (qc_data_all[which(qc_data_all$patient == Patientnum()),input$modality] == "Good") {
      tagList(
        (radioButtons("include", 
                      label = "How do you judge the Scan's quality?",
                      choices = c("Good","Moderate", "Poor"),
                      selected = "Good",
                      inline = TRUE
        )),
        (actionButton("update", "Confirm")))
    } else if (qc_data_all[which(qc_data_all$patient == Patientnum()),input$modality] == "Moderate") {
      tagList(
        (radioButtons("include", 
                      label = "How do you judge the Scan's quality?",
                      choices = c("Good","Moderate", "Poor"),
                      selected = "Moderate",
                      inline = TRUE
        )),
        (actionButton("update", "Confirm")))
    }
  })
  
  #  Finally, here we respond to the clicking of the button SAVE and we save a CSV in the dataframe with the new columns,
  #this csv will be read autoatically in future runs. (CHANGE COLUMNS NAME? so that is not just structural...)
  observeEvent(input$save, {write.csv(qc_data_all, file = "dataframes/QCed_data.csv")})
}
# Run the application 
shinyApp(ui = ui, server = server, enableBookmarking = "url")
