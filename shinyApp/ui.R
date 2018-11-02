#GUI of the shiny app
shinyUI(fluidPage(
  #Title
  titlePanel("OncoScan vs. cancer genes"),
  verticalLayout(
    mainPanel(
      #Short description
	  p('Crosses the results from an OncoScan run with the Cancer Gene Census (COSMIC database), a list of known 
        amplified/deleted cancer genes (Zack et al., Nature 2013), and the ', 
        a('CIVIC database',href = 'https://civic.genome.wustl.edu', target = '_blank'), '(downloaded on 11.07.2016).'),
      #File selection button
      fileInput('file1', 'Uploadez le resultat OncoScan (.txt)',
                accept = c(
                  'text/csv',
                  'text/comma-separated-values',
                  'text/tab-separated-values',
                  'text/plain',
                  '.csv',
                  '.tsv'
                )
      ),
	  #???
      tags$hr(),
	  #Print the main table
      tableOutput('contents')
    )
  )
))