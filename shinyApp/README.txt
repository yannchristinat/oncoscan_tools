App est accessible ici: https://hug-spc.shinyapps.io/shinyApp/

Admin interface:
https://www.shinyapps.io/admin
email: yann@hug, pwd: shinyHUG2018

Setup in R:
install.packages('rsconnect')
library(rsconnect)
options(RCurlOptions = list(proxy = "http://user:pwd@129.195.0.195:8080")) #replace user:pwd by your actual HUG initials and password
options(shinyapps.http = "rcurl")
rsconnect::setAccountInfo(name='hug-spc',
			  token='7C2FCBC52D558DAFC55FA9FB9136C46A',
			  secret='+kWtFp5w3g0hNI3+A4IMbjMShvLphR0xEcSWdvQ5')

Test:
library(shiny)
setwd('Q:/GitCentralRepo/CNV-tools/CNV-tools/oncoscan/shinyApp')
runApp()

Deploy:
library(rsconnect)
setwd('path/to/app/directory')
rsconnect::deployApp() #Proxy issues... TODO!
