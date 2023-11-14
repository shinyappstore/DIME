###########################################################################################

# App Name: DIME
# Author : Abhinandan Devaprasad
# Contact: a.devaprasad@umcutrecht.nl

###########################################################################################
###################################### INSTRUCTIONS #######################################
###########################################################################################

# Please make sure all of the following libraries are installed first.                                                         
# loading libraries
library(shiny)
library(shinyWidgets)
library(reshape2)
library(ggplot2)
library(DT)
library(NMF)
library(igraph)
library(qgraph)
library(pheatmap)
library(plyr)
library(RColorBrewer)
library(grid)

# To run the app from Rstudio, open this file (DIME.R) in Rstudio and simply click on the 
# 'Run App' option.                                                                      
#
# Make sure the accompanying file "DIME_shiny.Rdata" is in the same folder as the 'DIME.R'
# file.                                                                                   

###########################################################################################

# loading requisite data
load(file="DIME_shiny.rdata")

# UI page here
ui<-navbarPage("DIME",
               
               tabPanel("DisGeNet",
                        
                        sidebarLayout(
                          sidebarPanel(
                            helpText("To profile diseases based on their disease associated genes and multiple cell types, please enter the desired disease name in the box below and click Go!"),
                            helpText("Input is case insensitive."),
                            selectInput("ds", "Disease:", choices= names(table(dgaf_app$diseaseName)[table(dgaf_app$diseaseName)>10]),selected="Systemic scleroderma",multiple=FALSE),
                            actionButton("goButton", "Go!"),
                            progressBar(id = "dime", value = 0, total = 100, display_pct =T,striped = T,title = "Progress bar"),
                            helpText("The analysis can take about 5-10 minutes to run. In the meantime, please see the 'About' page to learn more on details of the result"),
                            helpText("Authors: ",a("Abhinandan Devaprasad,", href="mailto:a.devaprasad@umcutrecht.nl"),"Timothy Radstake,", a("Aridaman Pandit", href="mailto:a.pandit@umcutrecht.nl")),
                            helpText(a("Copyright (c) 2019, COSI Group", href="https://bitbucket.org/aridaman/cosi/wiki/Home")),
                            helpText("All rights reserved")
                            
                            
                          ),
                          
                          mainPanel(
                            h1("DIME Plot"),
                            plotOutput(outputId = "dime", height= "100%", width = "100%"),
                            h1("Download DIME Network"),
                            h5("Download the cytoscape-friendly DIME Network"),
                            downloadButton(outputId = "dl_disgenet_net_node", label="Download Node attributes"),
                            downloadButton(outputId = "dl_disgenet_net_edge", label="Download Edge attributes"),
                            h1("DIME Drug Network"),
                            plotOutput(outputId = "drug_net", height= "100%", width = "100%"),
                            helpText(a("Copyright (c) 2019, COSI Group", href="https://bitbucket.org/aridaman/cosi/wiki/Home")),
                            helpText("Authors: ",a("Abhinandan Devaprasad,", href="mailto:a.devaprasad@umcutrecht.nl"),"Timothy Radstake,", a("Aridaman Pandit", href="mailto:a.pandit@umcutrecht.nl")),
                            helpText("All rights reserved")
                          )
                        )
               ),
               
               tabPanel("EBI-GWAS",
                        
                        sidebarLayout(
                          sidebarPanel(
                            helpText("To profile diseases based on their disease associated genes and multiple cell types, please enter the desired disease name in the box below and click Go!"),
                            helpText("Input is case insensitive."),
                            selectInput("gwas_ds", "Disease:", choices= names(table(clean_dgnet2$d)[table(clean_dgnet2$d)>10]), selected="Systemic_sclerosis",multiple=FALSE),
                            actionButton("goButton2", "Go!"),
                            progressBar(id = "dime2", value = 0, total = 100, display_pct = T,striped = T,title = "Progress bar"),
                            helpText("The analysis can take about 5-10 minutes to run. In the meantime, please see the 'About' page to learn more on details of the result"),
                            helpText("Authors: ",a("Abhinandan Devaprasad,", href="mailto:a.devaprasad@umcutrecht.nl"),"Timothy Radstake,", a("Aridaman Pandit", href="mailto:a.pandit@umcutrecht.nl")),
                            helpText(a("Copyright (c) 2019, COSI Group", href="https://bitbucket.org/aridaman/cosi/wiki/Home")),
                            helpText("All rights reserved")
                            
                            
                          ),
                          
                          mainPanel(
                            h1("DIME Plot"),
                            plotOutput(outputId = "gwas_dime", height= "100%", width = "100%"),
                            h1("Download DIME Network"),
                            h5("Download the cytoscape-friendly DIME Network"),
                            downloadButton(outputId = "dl_gwas_net_node", label="Download Node attributes"),
                            downloadButton(outputId = "dl_gwas_net_edge", label="Download Edge attributes"),
                            h1("DIME Drug Network"),
                            plotOutput(outputId = "gwas_drug_net", height= "100%", width = "100%"),
                            helpText(a("Copyright (c) 2019, COSI Group", href="https://bitbucket.org/aridaman/cosi/wiki/Home")),
                            helpText("Authors: ",a("Abhinandan Devaprasad,", href="mailto:a.devaprasad@umcutrecht.nl"),"Timothy Radstake,", a("Aridaman Pandit", href="mailto:a.pandit@umcutrecht.nl")),
                            helpText("All rights reserved")
                          )
                        )
               ),
               
               
               tabPanel("Custom Gene Analysis",
                        
                        sidebarLayout(
                          sidebarPanel(
                            helpText("To do a custom gene analysis, please input your genes below and click Go!"),
                            helpText("NOTE: Minimum of 10 genes required to compute DIME"), 
                            helpText("Input is case insensitive."),
                            selectInput("cgi", "Genes:", choices = sort(unique(rownames(annot_mcpmk5fpval))), selected = NA, multiple = T),
                            actionButton("goButton1", "Go!"),
                            progressBar(id = "cga", value = 0, total = 100, display_pct = T,striped = T,title = "Progress bar"),
                            helpText("The analysis can take about 5-10 minutes to run. In the meantime, please see the 'About' page to learn more on details of the result"),
                            helpText(a("Copyright (c) 2019, COSI Group", href="https://bitbucket.org/aridaman/cosi/wiki/Home")),
                            helpText("Authors: ",a("Abhinandan Devaprasad,", href="mailto:a.devaprasad@umcutrecht.nl"),"Timothy Radstake,", a("Aridaman Pandit", href="mailto:a.pandit@umcutrecht.nl")),
                            helpText("All rights reserved")
                            
                          ),
                          mainPanel(
                            h1("DIME Plot"),
                            plotOutput(outputId = "cg_dime", height= "100%", width = "100%"),
                            h1("Download DIME Network"),
                            h5("Download the cytoscape-friendly DIME Network"),
                            downloadButton(outputId = "dl_cg_net_node", label="Download Node attributes"),
                            downloadButton(outputId = "dl_cg_net_edge", label="Download Edge attributes"),
                            h1("DIME Drug Network"),
                            plotOutput(outputId = "cg_drug_net", height= "100%", width = "100%"),
                            helpText(a("Copyright (c) 2019, COSI Group", href="https://bitbucket.org/aridaman/cosi/wiki/Home")),
                            helpText("Authors: ",a("Abhinandan Devaprasad,", href="mailto:a.devaprasad@umcutrecht.nl"),"Timothy Radstake,", a("Aridaman Pandit", href="mailto:a.pandit@umcutrecht.nl")),
                            helpText("All rights reserved")
                            
                          )
                        )
                        
               ),
               tabPanel("DisGeNet-explorer",
                        dataTableOutput("DiG")),
               
               tabPanel("EBI-GWAS-explorer",
                        dataTableOutput("gwas_DiG")),
               
               tabPanel("About",
                          mainPanel(
                            h2("Disease-gene Immune cell Expression (DIME) Network"),
                            br(),
                            h4("About"),
                            h5("The DIME network is a tool to assess the impact of disease genes on immune cells. It uses non-negative matrix factorization (NMF) technique to identify clusters of disease associated cells (DACs) that express the disease associated genes (DAGs). The DIME network in this tool has been implemented on two disease networks, the DisGeNet and the EBI GWAS network. It can also be used on custom set of genes to assess the impact of these genes on the immune cells."),
                            br(),
                            h4("The NMF model - DIME"),
                            h5("The NMF is performed on the input gene (DAGs) expression data of the immunome. The NMF splits the gene expression data into ‘k’ clusters. The estimation of the ‘k’ value is the most time-consuming process of the tool. This takes about 4-5 minutes to run, depending on the configuration of the machine. The top cluster (among the 'k' clusters) is identified by calculating the Frobenius norm of all the clusters; the cluster with the highest Frobenius norm is top cluster or rank1, subsequent clusters are rank2, rank3, etc. The rank1 cluster represents the most important cluster that maximally represents the data. Hence genes and cells of this cluster are regarded the most important for the given input genes (DAGs). Within each cluster, the genes and cell types are scored using the basis and coefficient matrices computed by the NMF. These scores range between 0 and 1. With 1 being the highest score. The genes and cell types in the top 25 percentile of the score are regarded as the key genes and cell types."),
                            br(),
                            h4("DIME Plot"),
                            h5("This informative heatmap comprises of the gene expression of all the input genes (DAGs) across all cell types of the immunome. The cell types are ordered based on their score in the rank1 (top) cluster. The genes are ranked according to the corresponding rank's score. The top 10 genes of each cluster are displayed in this plot."),
                            br(),
                            h4("DIME Drug Network"),
                            h5("The DIME drug network comprises of the key genes (as identified by DIME) and their corresponding drugs from the drug-gene network."),
                            br(),
                            h4("The DisGeNet Network"),
                            h5("The DisGeNet comprises of expert curated disease gene interactions from 16 different databases such as UNIPROT, CGI, ClinGen, Genomics England, CTD, PsyGeNET, Orphanet, RGD, MGD, CTD, Human Phenotype Ontology, Clinvar, GWAS catalogue, GWAS DB, LHDGN and BeFree. The full disease gene association network from DisGeNet was downloaded from the DisGeNet database (www.disgenet.org/downloads). All HLA associated genes was removed from the network, this was done to ensure that bias towards myeloid cells and B cells are removed, since the HLA genes are largely expressed by these cells. The resulting network was further filtered to include only those genes that were present in the immunome. The final network comprised of 15367 diseases and 13510 DAG."),
                            br(),
                            h4("The EBI-GWAS Network"),
                            h5("The GWAS catalogue of Version 1.0.2, e96, was downloaded from the EBI website. The catalogue also provided the corresponding mapped gene information for all the SNP which was used to construct the disease to gene association network. We further filtered this network using the same filtering criteria that was used for the DisGeNet. The final network comprised of 3670 diseases and 17810 DAG."),
                            br(),
                            h4("The Drug-gene Network"),
                            h5("The drug to gene target database from DGIdb was downloaded. The data was filtered to contain only the CHEMBL interactions and only those pertaining to the drugs approved by the food and drug administration (FDA) of USA. This FDA approved drug to gene target interaction serves as the drug-gene network in this study."),
                            br(),
                            h4("Application of DIME in IMIDs"),
                            h5("To know more on how we use DIME to identify the key cell types and genes of immune-mediate inflammatory diseases (IMIDs), read our pre-print article:"),
                            helpText(a("Integration of immunome with disease-gene network reveals common cellular mechanisms between IMIDs and drug repurposing strategies. Abhinandan Devaprasad, Timothy RDJ Radstake, Aridaman Pandit", href="https://doi.org/10.1101/2019.12.12.874321"))
                          )
                        )#)
)

# Server page here
server<-function(input, output, session) {
  
  # Functions
  # Scaling function
  scale01<-function(x){
    (x-min(x))/(max(x)-min(x))
  }
  
  # NMF Function
  nmfnet<-function(G,id){
    
    ID<-as.character(id)
    
    # Get immcell df
    if(ID!="maintenance"){updateProgressBar(session = session, id = ID, value = (1/5)*100, total = 100,title = "Get immunome")}
    #mtsetg<-as.data.frame(annot_mcpmk5fpval[row.names(annot_mcpmk5fpval) %in% G,cellorder])
    mtsetg<-as.data.frame(annot_mcpmk5fpval[match(G, row.names(annot_mcpmk5fpval)),cellorder])
    mtsetg<-as.matrix(2^mtsetg)
    
    # Estimating ideal Rank - 'K' clusters - Most time consuming step!
    if(ID!="maintenance"){updateProgressBar(session = session, id = ID, value = (2/5)*100, total = 100, title = "Estimating Rank,takes time!")}
    set.seed(123456) # Set seed to get reproducible results
    nmf.options(pbackend=NA)
    nmf.options(gc=5)
    estNMFR<-nmfEstimateRank(mtsetg, seq(2,6), method='brunet', nrun=30, seed=123456) # ideal runs should be over 30
    tmp_ranks<-c()
    for(i in 1:(length(estNMFR$measures$cophenetic)-1)){
      if(estNMFR$measures$cophenetic[i]>estNMFR$measures$cophenetic[i+1]){
        tmp_ranks<-c(tmp_ranks, (i+1))
      }
    }
    ideal_rank<-min(tmp_ranks)
    # ideal_rank=3 # For testing purposes
    
    # Calculate NMF
    if(ID!="maintenance"){updateProgressBar(session = session, id = ID, value = (3/5)*100, total = 100, title = "Compute NMF")}
    res<-nmf(mtsetg, rank = ideal_rank, method = "brunet", seed="nndsvd")
    
    # Get W and B matrix
    w<-basis(res) ; w<-as.matrix(apply(w, 2, scale01))
    h<-coef(res) ; h<-as.matrix(t(apply(h, 1, scale01)))

    # Choose feature vector (Top weighted cluster)
    if(ID!="maintenance"){updateProgressBar(session = session, id = ID, value = (4/5)*100, total = 100, title = "Identify key feature")}
    norms<-c()
    for(r in 1:ideal_rank){
      tmp_norm<-norm(as.matrix(h[r,])%*%t(as.matrix(w[,r])), type="F")
      norms<-c(norms, tmp_norm)
    }
    cfv<-which.max(norms)
    
    # Save results
    if(ID!="maintenance"){updateProgressBar(session = session, id = ID, value = (5/5)*100, total = 100, title = "Gather result + Plot")}
    out<-list(w=w, h=h, cfv=cfv, ideal_rank=ideal_rank)
    
  }
  
  # Plotting DIME SAUCEPLOT (informative plot)
  
  # addflag function for pheatmap
  add.flag <- function(pheatmap,
                       kept.labels,
                       repel.degree){
    
    heatmap <- pheatmap$gtable
    
    new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
    
    # keep only labels in kept.labels, replace the rest with ""
    new.label$label <- ifelse(new.label$label %in% kept.labels, 
                              new.label$label, "")
    
    # calculate evenly spaced out y-axis positions
    repelled.y <- function(d, d.select, k = repel.degree){
      # recursive function to get current label positions
      strip.npc <- function(dd){
        if(!"unit.arithmetic" %in% class(dd)) {
          return(as.numeric(dd))
        }
        
        d1 <- strip.npc(dd$arg1)
        d2 <- strip.npc(dd$arg2)
        fn <- dd$fname
        return(lazyeval::lazy_eval(paste(d1, fn, d2)))
      }
      
      full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
      selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
      
      return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                      to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                      length.out = sum(d.select)), 
                  "npc"))
    }
    new.y.positions <- repelled.y(new.label$y,
                                  d.select = new.label$label != "")
    new.flag <- segmentsGrob(x0 = new.label$x,
                             x1 = new.label$x + unit(0.15, "npc"),
                             y0 = new.label$y[new.label$label != ""],
                             y1 = new.y.positions)
    
    # shift position for selected labels
    new.label$x <- new.label$x + unit(0.2, "npc")
    new.label$y[new.label$label != ""] <- new.y.positions
    
    # add flag to heatmap
    heatmap <- gtable::gtable_add_grob(x = heatmap,
                                       grobs = new.flag,
                                       t = 4, 
                                       l = 4
    )
    
    # replace label positions in heatmap
    heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
    
    # plot result
    invisible(heatmap)
  }
  
  # DIME PLOT
  plot_dime<-function(w,h,cfv,idealranks){
    
    cutoff = 0
    ngenes = 10
    # Arrange based on ranks
    ws<-as.matrix(apply(w, 2, scale01))
    hs<-as.matrix(t(apply(h, 1, scale01)))
    
    # Calculate Frob norms
    rankassign<-data.frame()
    for(r in 1:idealranks){
      tmp_norm<-data.frame(norm=norm(as.matrix(hs[r,])%*%t(as.matrix(ws[,r])), type="F"), oldrank=r)
      rankassign<-rbind(rankassign, tmp_norm)
    }
    rankassign$newrank<-order(rankassign$norm, decreasing = T)
    
    # Rearrange based on top ranks
    w<-w[,rankassign$newrank]
    h<-h[rankassign$newrank,]
    
    # Identifying key DACs and DAGs
    w<-cbind(w, apply(w,1,max))
    h<-t(h)
    h<-as.data.frame(h)
    h$max<-as.data.frame(apply(h,1,max))[,1]
  
    dagdf<-data.frame()
    for(i in 1:idealranks){
      ind<-intersect(which(w[,i]==w[,idealranks+1]), which(w[,i]>=quantile(w[,i], c(cutoff))))
      tmp<-data.frame(dag=row.names(w)[ind], rank=i)
      dagdf<-rbind(dagdf, tmp)
    }
    
    # Extracting scores
    dagdf$score<-w[match(dagdf$dag, rownames(w)),idealranks+1]
    dagdf<-ddply(dagdf,"dag",numcolwise(min))
    
    # Adding annotations
    annot_row<-dagdf; rownames(annot_row)<-annot_row$dag; annot_row$dag<-NULL
    annot_row$rank<-as.character(annot_row$rank)
    
    annot_col<-h[match(cellorder, row.names(h)),rev(1:idealranks)]
    colnames(annot_col)<-paste("score_rank",rev(1:idealranks), sep="_")
  
    rankcol<-brewer.pal(6,"Set1")[c(1:idealranks)]
    names(rankcol)<-1:idealranks
    
    annot_row_clr<-list(rank=rankcol)
    
    # Labelling only top 'n' key genes
    keygenes<-c()
    for(i in 1:idealranks){
      tmp<-as.character(head(dagdf[dagdf$rank==i,][order(dagdf[dagdf$rank==i,"score"], decreasing = T),"dag"], ngenes))
      keygenes<-c(keygenes, tmp)
    }
    
    annot_row<-annot_row[order(annot_row$rank,annot_row$score, decreasing = T),]
    ar1<-data.frame()
    for(r in 1:idealranks){
      tmpar<-annot_row[annot_row$rank==r,]
      tmpar<-tmpar[order(tmpar$score, decreasing = T),]
      ar1<-rbind(ar1, tmpar)
    }
    
    annot_row<-ar1
    
    dgcdf<-as.data.frame(annot_mcpmk5fpval[match(row.names(annot_row), row.names(annot_mcpmk5fpval)), cellorder[order(annot_col$score_rank_1, decreasing = T)]])
    rownames(dgcdf)<-row.names(annot_row)
    
    # Plotting using pheatmap
    pheat<-pheatmap(dgcdf,
                    annotation_row = annot_row, annotation_col = annot_col,
                    border_color = NA, cluster_rows = F, cluster_cols = F, annotation_colors = annot_row_clr,
                    silent=T)
    # Adding labels
    add.flag(pheat,
             kept.labels = keygenes,
             repel.degree = 0.5)
  }
  
  # Net Extracter
  
  netext<-function(w,h,idealranks){
    cutoff=0.75
    # Arrange based on ranks
    ws<-as.matrix(apply(w, 2, scale01))
    hs<-as.matrix(t(apply(h, 1, scale01)))
    
    # Calculate Frob norms
    rankassign<-data.frame()
    for(r in 1:idealranks){
      tmp_norm<-data.frame(norm=norm(as.matrix(hs[r,])%*%t(as.matrix(ws[,r])), type="F"), oldrank=r)
      rankassign<-rbind(rankassign, tmp_norm)
    }
    rankassign$newrank<-order(rankassign$norm, decreasing = T)
    
    # Rearrange based on top ranks
    w<-w[,rankassign$newrank]
    h<-h[rankassign$newrank,]
    
    w<-cbind(w, apply(w,1,max))
    h<-t(h)
    h<-cbind(h, apply(h,1,max))  
    
    w<-as.data.frame(w); h<-as.data.frame(h)
    # net
    edges<-data.frame()
    nodes<-data.frame()
    for(i in 1:idealranks){
      dag=row.names(w)[intersect(which(w[,i]==w[,idealranks+1]), which(w[,i]>=quantile(w[,i], c(cutoff))))]
      dac=row.names(h)[intersect(which(h[,i]==h[,idealranks+1]), which(h[,i]>=quantile(h[,i], c(cutoff))))]
      tmp<-expand.grid(dag, dac)
      if(nrow(tmp)>0){
        tmp$rank<-as.character(i)
        edges<-rbind(edges, tmp)
      }
      
      for(g in dag){
        tmp1<-data.frame(node=g, class="Gene", rank=i, score=w[row.names(w)==g,i])
        nodes<-rbind(nodes,tmp1)
      }
      for(c in dac){
        tmp1<-data.frame(node=c, class="Cell", rank=i, score=h[row.names(h)==c,i])
        nodes<-rbind(nodes,tmp1)
      }
      
    }
    
    # Net
    out<-list(edges=edges, nodes=nodes)
  }
  

  # Plotting dime drug network
  plot_drug_net<-function(w,h,ideal_rank){
    cutoff=0.75
    w<-as.data.frame(w)
    w$max<-as.data.frame(apply(w,1,max))[,1]
    h<-t(h)
    h<-as.data.frame(h)
    h$max<-as.data.frame(apply(h,1,max))[,1]
    
    dag<-c()
    for(rank in 1:ideal_rank){
      dag<-c(row.names(w[intersect(which(w[,rank]==w$max), which(w[,rank]>=quantile(w[,rank], c(cutoff)))),]), dag)
    }
    
    dg<-fda_genes_net[fda_genes_net$V1 %in% dag,]
    if(nrow(dg)>=1){
      g_fda_d1d2_imm_net<-graph_from_edgelist(as.matrix(dg), directed = F)
      l<-qgraph.layout.fruchtermanreingold(as.data.frame(get.edgelist(g_fda_d1d2_imm_net, names = F)),vcount=vcount(g_fda_d1d2_imm_net),
                                           area=(2*vcount(g_fda_d1d2_imm_net)^2),repulse.rad=(2*vcount(g_fda_d1d2_imm_net)^2))
      
      plot(g_fda_d1d2_imm_net, layout=l,
           vertex.size=4,
           vertex.shape=ifelse(vertex.attributes(g_fda_d1d2_imm_net)$name %in% dg$V1, "square", "circle" ),
           vertex.color=ifelse(vertex.attributes(g_fda_d1d2_imm_net)$name %in% dg$V1, "greenyellow", "burlywood1" ),
           edge.arrow.size=0, vertex.label.color="black", vertex.frame.color=NA, edge.color="lightgrey",
           vertex.label.cex=0.65)
      
    } else {
      par(mar = c(0,0,0,0))
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, paste("No Drugs found"), 
           cex = 1.6, col = "black")
      
    }
  }
  
  # End of functions #

  # Output section
  ### Disease Network page - DisGeNet page
  observeEvent(input$goButton,{
    # Get Input
    d<-eventReactive(input$goButton,{paste0((input$ds))})
    genes<-eventReactive(input$goButton,{as.character(unique(dgaf_app[dgaf_app$diseaseName==d(),"geneSymbol"]))})
    
    # Run NMF
    nmf_out<-eventReactive(input$goButton,{
      nmfnet(G = genes(),id = "dime")
    })
    
    plotting<-reactive({
      req(input$goButton)
      plots<-plot_dime(w=nmf_out()$w, h=nmf_out()$h, idealranks = nmf_out()$ideal_rank, cfv = nmf_out()$cfv)
      return(plots)
    })
    
    # DIME plot
    output$dime <- renderPlot({
      req(input$goButton)
      grid.newpage()
      grid.draw(plotting())
    }, height = 1050, width=1050)
    
    # Download DIME Net
    
    net<-eventReactive(input$goButton,{
      netext(w=nmf_out()$w, h=nmf_out()$h, idealranks=nmf_out()$ideal_rank)
    })
    
    output$dl_disgenet_net_node<-downloadHandler(
      req(input$goButton),
      filename = function(){
        paste("DIME_disgenet_net_node_", gsub(" ", "_", d()), ".txt", sep="")
      },
      content = function(file){
        write.table(net()$nodes, file, row.names = F, quote = F)
      }
    )
    
    output$dl_disgenet_net_edge<-downloadHandler(
      req(input$goButton),
      filename = function(){
        paste("DIME_disgenet_net_edge_", gsub(" ", "_", d()), ".txt", sep="")
      },
      content = function(file){
        write.table(net()$edges, file, row.names = F, quote = F)
      }
    )

    # DIME drug net
    output$drug_net <- renderPlot({
      req(input$goButton)
      plots<-plot_drug_net(nmf_out()$w, nmf_out()$h, nmf_out()$ideal_rank)
      plots
    }, height = 1050, width=1050)
    
  })
  
  ### Disease Network page - GWAS page
  observeEvent(input$goButton2,{
    # Get Input
    gwas_d<-eventReactive(input$goButton2,{paste0((input$gwas_ds))})
    gwas_genes<-eventReactive(input$goButton2,{as.character(unique(clean_dgnet2[clean_dgnet2$d==gwas_d(),"g"]))})
    
    # Run NMF
    gwas_nmf_out<-eventReactive(input$goButton2,{
      nmfnet(G = gwas_genes(),id = "dime2")
    })
    
    gwas_plotting<-reactive({
      req(input$goButton2)
      plots<-plot_dime(w=gwas_nmf_out()$w, h=gwas_nmf_out()$h, idealranks = gwas_nmf_out()$ideal_rank, cfv = gwas_nmf_out()$cfv)
      return(plots)
    })
    
    # DIME plot
    output$gwas_dime <- renderPlot({
      req(input$goButton2)
      grid.newpage()
      grid.draw(gwas_plotting())
    }, height = 1050, width=1050)
    
    # Download DIME Net
    
    gwas_net<-eventReactive(input$goButton2,{
      netext(w=gwas_nmf_out()$w, h=gwas_nmf_out()$h, idealranks=gwas_nmf_out()$ideal_rank)
    })
    
    output$dl_gwas_net_node<-downloadHandler(
      req(input$goButton2),
      filename = function(){
        paste("DIME_gwas_net_node_", gsub(" ", "_", gwas_d()), ".txt", sep="")
      },
      content = function(file){
        write.table(gwas_net()$nodes, file, row.names = F, quote = F)
      }
    )
    
    output$dl_gwas_net_edge<-downloadHandler(
      req(input$goButton2),
      filename = function(){
        paste("DIME_gwas_net_edge_", gsub(" ", "_", gwas_d()), ".txt", sep="")
      },
      content = function(file){
        write.table(gwas_net()$edges, file, row.names = F, quote = F)
      }
    )
    
    
    # DIME Drug Net
    output$gwas_drug_net <- renderPlot({
      req(input$goButton2)
      plots<-plot_drug_net(gwas_nmf_out()$w, gwas_nmf_out()$h, gwas_nmf_out()$ideal_rank)
      plots
    }, height = 1050, width=1050)
    
  })

  ### Custon Gene Analysis Panel
  
  observeEvent(input$goButton1,{
    # Read input
    cg<-eventReactive(input$goButton1,{as.character(input$cgi)})
    cg_d<-eventReactive(input$goButton1,{paste0(("Custom_Net"))})
    
    # Run NMF
    cg_nmf_out<-eventReactive(input$goButton1,{
      nmfnet(G = cg(),id ="cga")
    })
    
    cg_plotting<-reactive({
      req(input$goButton1)
      plots<-plot_dime(w=cg_nmf_out()$w, h=cg_nmf_out()$h, idealranks = cg_nmf_out()$ideal_rank, cfv = cg_nmf_out()$cfv)
      return(plots)
    })
    
    # DIME plot
    output$cg_dime <- renderPlot({
      req(input$goButton1)
      grid.newpage()
      grid.draw(cg_plotting())
    }, height = 1050, width=1050)
    
    # Download DIME Net
    
    cg_net<-eventReactive(input$goButton1,{
      netext(w=cg_nmf_out()$w, h=cg_nmf_out()$h, idealranks=cg_nmf_out()$ideal_rank)
    })
    
    output$dl_cg_net_node<-downloadHandler(
      req(input$goButton1),
      filename = function(){
        paste("DIME_cg_net_node_", gsub(" ", "_", cg_d()), ".txt", sep="")
      },
      content = function(file){
        write.table(cg_net()$nodes, file, row.names = F, quote = F)
      }
    )
    
    output$dl_cg_net_edge<-downloadHandler(
      req(input$goButton1),
      filename = function(){
        paste("DIME_cg_net_edge_", gsub(" ", "_", cg_d()), ".txt", sep="")
      },
      content = function(file){
        write.table(cg_net()$edges, file, row.names = F, quote = F)
      }
    )
    
    # DIME Drug net
    output$cg_drug_net <- renderPlot({
      req(input$goButton1)
      plots<-plot_drug_net(cg_nmf_out()$w, cg_nmf_out()$h, cg_nmf_out()$ideal_rank)
      plots
    }, height = 1050, width=1050)
    
  })

  # DisGeNet Network Panel
  output$DiG<-renderDataTable({
    datatable(dgaf_app, rownames = F, filter = 'top', colnames = c("GENES", "DISEASES"), class='cell-border stripe')
  })

  # GWAS Network Panel
  output$gwas_DiG<-renderDataTable({
    datatable(clean_dgnet2[,c("g","d")], rownames = F, filter = 'top', colnames = c("GENES", "DISEASES"), class='cell-border stripe')
  })
  
}

# Compiling shiny app ui and server
shinyApp(ui = ui, server = server)