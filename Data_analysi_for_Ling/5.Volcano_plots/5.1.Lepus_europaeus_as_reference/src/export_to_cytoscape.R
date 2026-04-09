#function to generate network files for cytoscape

export_to_cytoscape <- function(network = bwnet,
                                modules,
                                edge_threshold = 0.02
                                ){
  #load TOM
  load(network$TOMFiles)
  TOM <- as.matrix(TOM)
  
  # Select module genes 
  inModule=is.finite(match(network$colors,modules)) 
  modGenes=names(network$colors)[inModule] 
  
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modGenes, modGenes)
  
  # Export the network into edge and node list files for Cytoscape 
  cyt = exportNetworkToCytoscape(modTOM, 
                                 edgeFile=paste("out/cytoscape/CytoEdge",paste(modules,collapse="-"),".txt",sep=""), 
                                 nodeFile=paste("out/cytoscape/CytoNode",paste(modules,collapse="-"),".txt",sep=""),
                                 weighted = TRUE, 
                                 threshold = edge_threshold,
                                 nodeNames = modGenes, 
                                 nodeAttr = bwnet$colors[inModule])
  
  return("edge and node files saved in out/cytoscape folder")
  }