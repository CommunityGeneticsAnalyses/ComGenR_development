###Script to install the development version of ComGenR

.instpacks <- as.character(installed.packages()[,1])
if (all(.instpacks!='devtools')){install.packages('devtools')}
library('devtools')
install_github('ComGenR',user='MKLau')
library('ComGenR')
