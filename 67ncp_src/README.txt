README.txt 13 Jan 2020
Jeffrey R. Stevens (jeffrey.r.stevens@gmail.com)

This file describes how to reproduce the results from the following paper: Duque, J.F., Rasmussen, T., Rodriguez, A., & Stevens, J.R. (2020). Effects of mesotocin on social bonding in pinyon jays. Ethology. 126, 165-175. doi:10.1111/eth.12990

All materials presented here are released under the Creative Commons Attribution 4.0 International Public License (CC BY 4.0). You are free to:
    Share — copy and redistribute the material in any medium or format
    Adapt — remix, transform, and build upon the material for any purpose, even commercially.
Under the following terms:
    Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.
    No additional restrictions — You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.

To reproduce these results, first unzip duque_etal_2020_rr.zip into a folder.  Then, ensure that a subfolder named "figures" is in the folder with all other files.  Next, open duque_etal_2020_rcode.R and ensure that all packages mentioned at the top of the script are installed.  Once all packages are installed, run the script in R using "source("duque_etal_2020_rcode.R")".  

Once the script runs without errors, you can compile the RMarkdown document duque_etal_2020.Rmd.  Open this file in RStudio and ensure that you have packages 'knitr' and 'rmarkdown' installed.  Once installed, use knitr to compile the document (control-shift-k).  Use the same process to compile duque_etal_2020_SM.Rmd.
