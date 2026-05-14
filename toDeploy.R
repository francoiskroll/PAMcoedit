# install.packages('rsconnect')
library(rsconnect)
 
# to configure rsconnect so it links to your online repo: in shinyapps.io dashboard, go below your avatar > Tokens
# then Show > Show Secret
# then copy command paste/run in R Console
deployApp()
