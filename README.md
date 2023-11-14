## Disease-gene IMmune cell Expression (DIME)

### About

DIME is a R Shinyapp, that can be used to assess the impact of disease genes on immune cells. It uses non-negative matrix factorization (NMF) technique to identify clusters/network of disease associated cells (DAC) that express the disease associated genes (DAG). The DIME network in this tool has been implemented on two disease networks, the DisGeNet and the EBI GWAS network. It can also be used on custom set of genes to assess the impact of these genes on the immune cells.

---

### Installation

To use the DIME app, please install all of the following dependencies first:

**Software**

1. [R](https://www.r-project.org/)
2. [Rstudio](https://rstudio.com/)

Once you have installed R and Rstudio, proceed to install the following packages in R:

**R Packages**

1. shiny
2. shinyWidgets
3. reshape2
4. ggplot2
5. DT
6. NMF
7. igraph
8. qgraph
9. pheatmap
10. plyr
11. RColorBrewer
12. grid

```
# Use the following code to install the packages in R
install.packages(c("shiny", "shinyWidgets", "reshape2", "ggplot2", "DT", "NMF", "igraph", "qgraph", "pheatmap", "plyr", "RColorBrewer", "grid"), dependencies = T)
```
**Installing DIME**

To install the DIME app open a shell terminal and type the following command:

```
cd ~
git clone https://bitbucket.org/systemsimmunology/dime.git

```
The installation is now complete!

Now, To run the DIME app, type the following commands
```
cd ~/dime/
R -e "shiny::runApp('DIME.R', launch.browser=T)"
```
DIME should open in your default browser. The app takes a minute to load. Once done you can begin your analysis. 

Please see the about page in the app to get more information about the analysis and the results.

DIME has been used to analysis the disease associated cell types and the key genes involved in immune-mediated inflammatory diseases (IMIDs). Read more about our work in our pre-print: 

```
Integration of immunome with disease-gene network reveals common cellular mechanisms between IMIDs and drug repurposing strategies
Abhinandan Devaprasad, Timothy RDJ Radstake, Aridaman Pandit
doi: https://doi.org/10.1101/2019.12.12.874321
```
