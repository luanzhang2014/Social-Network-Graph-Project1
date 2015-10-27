library("igraph")
library("fit.models")
library("ggplot2")
library("hash")

edgelistFile <- "/Users/luanzhang/Downloads/facebook_combined.txt"

g <- read.graph(edgelistFile , format = "ncol" , directed=FALSE)

############   QUESTION_1     ################
is.connected(g)
is.directed(g)

diameter(g, directed = FALSE)
h = hist(degree(g), breaks=seq(0.0, by=1 , length.out=max(degree(g))+2))       
df = data.frame(x=h$mids, y=h$density)
plot(df,main="Degree Distribution of Facebook Graph", xlab="Nodes", ylab="Degree Distribution",type="o")

models <- list(
  nls(y ~ (1/x*a) + b*x, data = df, start = list(a = 0, b = 0),trace=T), 
  nls(y ~ (a + b*log(x)), data=df, start = list(a = 0, b = 0),trace=T),
  nls(y ~ (exp(a + b * x)), data=df, start = list(a=0,b=0),trace=T),
  nls(y ~ (1/x*a)+b, data=df, start = list(a=1,b=1)),trace=T)


ggplot(df, aes(x, y)) + geom_point(size = 2) +
  geom_line(aes(x,fitted(models[[1]])), size = 1,colour = "blue") + 
  geom_line(aes(x,fitted(models[[2]])), size = 1, colour = "yellow") +
  geom_line(aes(x,fitted(models[[3]])),size = 1,  colour = "red") +
  geom_line(aes(x,fitted(models[[4]])),size = 1,  colour = "purple")+
  ggtitle("Fitted curves for degree distribution")+ xlab("Nodes") +ylab("Degree Distribution")

summary(models[[3]])
models[[3]]

mean(degree(g))

############  END OF QUESTION_1   #############


############   QUESTION_2     ################

subGraphNodes <- neighborhood(g , 1 , 1)
subGraphNodes <- subGraphNodes[[1]]
subGraph <- induced.subgraph(g,which( ( (1:vcount(g)) %in% subGraphNodes)  ))
length(E(subGraph))  ##no_of_edges
length(V(subGraph))  ##no_of_nodes 


############  END OF QUESTION_2   #############


############   QUESTION_3     ###########

core_nodes = {}
core_nodes<-which(neighborhood.size(g, 1) > 200)
length(core_nodes)  #number of core_nodes
mean(degreeDist[core_nodes])  #average degree

core <- core_nodes[10]
subGraphCoreNodes <- neighborhood(g , 1 , nodes=core)
subGraphCoreNodes <- subGraphCoreNodes[[1]]
subGraphCoreNodes <- induced.subgraph(g,which( ( (1:vcount(g)) %in% subGraphCoreNodes)  ))
plot(subGraphCoreNodes,main="Sub Network of the core node 10")

comm1 = fastgreedy.community(subGraphCoreNodes)
plot(comm1,subGraphCoreNodes,main="Communities for Core Graph's Node 10, FAST-GREEDY",vertex.size=5,vertex.label=NA)
comm2 <- edge.betweenness.community(subGraphCoreNodes)
plot(comm2, subGraphCoreNodes,main="Communities for Core Graph's Node 10, EDGE-BETWEENNESS",,vertex.size=5,vertex.label=NA)
comm3 <- infomap.community(subGraphCoreNodes)
plot(comm3, subGraphCoreNodes,main="Communities for Core Graph's Node 10, INFOMAP",,vertex.size=5,vertex.label=NA)
hist(comm1$membership,col="dark red",main="Community Distribution of the Fast Greedy Algorithm",xlab="Community Number",ylab="No. of Nodes in a Community")
hist(comm2$membership,col="dark blue",main="Community Distribution of the Edge-Betweenness Algorithm",xlab="Community Number",ylab="No. of Nodes in a Community")
hist(comm3$membership,col="dark green",main="Community Distribution of the Infomap Algorithm",xlab="Community Number",ylab="No. of Nodes in a Community")


#plot(subGraph_core,vertex.color="black",edge.color="gray",vertex.size=4,vertex.label=NA,edge.arrow.size=0.01)


############  END OF QUESTION_3   #############


############   QUESTION_4     ###########
core <- core_nodes[10]
subGraphCoreNodes <- neighborhood(g , 1 , nodes=core)
subGraphCoreNodes <- subGraphCoreNodes[[1]]
subGraphCoreNodes1 <- subGraphCoreNodes[-1]
subGraphCoreNodes1 <- induced.subgraph(g,which( ( (1:vcount(g)) %in% subGraphCoreNodes1)  ))
plot(subGraphCoreNodes1,main="Sub Network of the core node 10 (removing the core node itself)")

comm1_1 = fastgreedy.community(subGraphCoreNodes1)
plot(comm1_1,subGraphCoreNodes1,main="Communities for Core Graph without core node 10, FAST-GREEDY")
comm2_1 <- edge.betweenness.community(subGraphCoreNodes1)
plot(comm2_1, subGraphCoreNodes1,main="Core Node 10's Personal Network Structure Without Node 10, EDGE-BETWEENNESS")
comm3_1 <- infomap.community(subGraphCoreNodes1)
plot(comm3_1, subGraphCoreNodes1,main="Communities for Core Graph without core node 10, INFOMAP")
hist(comm1_1$membership,col="dark red",main="Community Distribution of the Fast Greedy Algorithm",xlab="Community Number",ylab="No. of Nodes in a Community")
hist(comm2_1$membership,col="dark blue",main="Community Distribution of the Edge-Betweenness Algorithm",xlab="Community Number",ylab="No. of Nodes in a Community")
hist(comm3_1$membership,col="dark green",main="Community Distribution of the Infomap Algorithm",xlab="Community Number",ylab="No. of Nodes in a Community")
############  END OF QUESTION_4   #############


############   QUESTION_5     ###########

# defining initial functions
# ===========================
commNeib_find <- function(u,v,g)
{
  neighborsU <- neighborhood(g,1,u)[[1]][-1]
  neighborsV <- neighborhood(g,1,v)[[1]][-1]
  intersect(neighborsU, neighborsV)
}

# embeddedness calculations
# ===========================
embd_find <- function (u,v,g)
{
  embd = length(commNeib_find(u,v,g))
  embd
}

perNet_find <- function(u, g)
{
  pnNodes <- neighborhood(g , 1 , nodes=u)[[1]]
  nonPNNodes <- which( !( (1:vcount(g)) %in% pnNodes)  )
  perNetw <- delete.vertices(g , nonPNNodes)
  perNetw$name =  sort(pnNodes)
  perNetw
}

# node ids
# ========
nodeID_find <- function(g, vertex)
{
  temp <- which(g$name == vertex)
  temp
}

# dispersion
# ==========
disp_find <- function(u,v,g) 
{
  disp <- 0
  commonUV <- commNeib_find(u, v, g)
  gNoUV <- delete.vertices(g, c(u, v))
  
  for(s in commonUV) 
  {
    for(t in commonUV) 
    {
      if(s != t && s!=u && s!=v && t!=u && t!=v) 
      {
        short_d<-get.shortest.paths(gNoUV,from=s,to=t)
        if(length(short_d$vpath[[1]])>0){
        d<-length(short_d$vpath[[1]])-1
        disp <- disp + d}
        
      }
    }
  }
  disp=disp/2
}

# ratio, emb, disp
# =====
dispEmb_find <- function(g,coreNode){
  
  dHigh = 0;
  dNode = 0;
  rHigh = 0;
  rNode = 0;
  eHigh = 0;
  eNode = 0;
  
  pnOfU <- perNet_find(coreNode,g)
  u <- nodeID_find(pnOfU, coreNode)
  
  nodes <- V(pnOfU)
  for(v in nodes){
    if(v == u)
      next
    
    dip = disp_find(u,v,g)
    embd = embd_find(u,v,g)
    
    if (embd > 0)
    {
      rt = dip/embd
      if (rt > rHigh)
      {
        rNode = v;
        rHigh=rt;
      }
    }
    
    if (dip > dHigh)
    {
      dNode = v;
      dHigh=dip;
    }
    if (embd > eHigh)
    {
      eNode = v
      eHigh=embd;
    }
    
  }
  
#figure 1  
#=============
  if (dNode > 0)
  {
    # community detection
    fc = fastgreedy.community(pnOfU); sizes(fc)
    mfc = membership(fc)
    
    sizeVet = rep(3, length(V(pnOfU)));
    sizeVet[dNode] = 8;  
    colEd = rep(8, length(E(pnOfU)));
    colEd[which(get.edgelist(pnOfU,name=F)[,1] == dNode | get.edgelist(pnOfU,name=F)[,2] == dNode)] = 3;
    E(pnOfU)$color = colEd;
    widEd = rep(1, length(E(pnOfU)));
    widEd[which(get.edgelist(pnOfU,name=F)[,1] == dNode | get.edgelist(pnOfU,name=F)[,2] == dNode)] = 3;
    dev.new ();
    plot(pnOfU, vertex.label= NA, vertex.color=mfc,vertex.size=sizeVet, edge.width = widEd,mark.groups = by(seq_along(mfc), mfc, invisible),main="Max dispersion");
    }
  
  else
  {
    print (paste(c("No high Disp node", toString(coreNode)), collapse=" "));
  }

#figure 2  
#=============  
  if (eNode > 0)
  {
    # community detection
    fc = fastgreedy.community(pnOfU); sizes(fc)
    mfc = membership(fc)
    sizeVet = rep(3, length(V(pnOfU)));
    sizeVet[eNode] = 8;  
    colEd = rep(8, length(E(pnOfU)));
    colEd[which(get.edgelist(pnOfU,name=F)[,1] == eNode | get.edgelist(pnOfU,name=F)[,2] == eNode)] = 3;
    E(pnOfU)$color = colEd;
    widEd = rep(1, length(E(pnOfU)));
    widEd[which(get.edgelist(pnOfU,name=F)[,1] == eNode | get.edgelist(pnOfU,name=F)[,2] == eNode)] = 3;
    dev.new ();
    plot(pnOfU, vertex.label= NA, vertex.color=mfc,vertex.size=sizeVet, edge.width = widEd,mark.groups = by(seq_along(mfc), mfc, invisible),main="Max embeddedness");# ,mark.groups = by(seq_along(mfc), mfc) );
    }
  else
  {
    print (paste(c("No high Emb node", toString(coreNode)), collapse=" "));
  }

#figure 3  
#=============
  if (rNode > 0)
  {

    # community detection
    fc = fastgreedy.community(pnOfU); sizes(fc)
    mfc = membership(fc)
    
    sizeVet = rep(3, length(V(pnOfU)));
    sizeVet[rNode] = 8;  
    colEd = rep(8, length(E(pnOfU)));
    colEd[which(get.edgelist(pnOfU,name=F)[,1] == rNode | get.edgelist(pnOfU,name=F)[,2] == rNode)] = 3;
    E(pnOfU)$color = colEd;
    widEd = rep(1, length(E(pnOfU)));
    widEd[which(get.edgelist(pnOfU,name=F)[,1] == rNode | get.edgelist(pnOfU,name=F)[,2] == rNode)] = 3;
    dev.new ();
    plot(pnOfU, vertex.label= NA, vertex.color=mfc,vertex.size=sizeVet, edge.width = widEd,mark.groups = by(seq_along(mfc), mfc, invisible) , main="Max dispersion/embeddedness");
    }
  else
  {
    print (paste(c("No high Disp node", toString(coreNode)), collapse=" "));
  }
  
}

# =======================

dispVec <- c();
embVec <- c();


for(coreNode in core_nodes)
{
  pnOfU <- perNet_find(coreNode,g)
  u <- nodeID_find(pnOfU, coreNode)
  
  nodes <- V(pnOfU)
  for(v in nodes){
    if(v == u)
      next
    
    embd = embd_find(u,v,g)
    dip = disp_find(u,v,g)
    embVec <- c(embVec, embd);
    dispVec <- c(dispVec, dip);
    
  }
}


hist (embVec, breaks=seq (-0.5, by=1, length.out=max(embVec) +2), main ="embd Distribution", xlab="embd");
hist (dispVec, breaks=seq (-0.5, by=1, length.out=max(dispVec) +2), main="dip Distribution", xlab="dip");

dispEmb_find(g,core_nodes[1])
dispEmb_find(g,core_nodes[9])
dispEmb_find(g,core_nodes[10])


############  END OF QUESTION_5   #############


############   QUESTION_6     ###########

   get_Gu <- function(u, g){
   subGraphNodes <- neighborhood(g , 1 , nodes=u)
   subGraphNodes <- subGraphNodes[[1]]
   nonSubGraphNodes <- which( !( (1:vcount(g)) %in% subGraphNodes)  )
   subGraph <- delete.vertices(g , nonSubGraphNodes)
   subGraph$name =  sort(subGraphNodes)
   subGraph
 }

h_main <- hash(keys=1,values=1)
for(node in core_nodes)
{
  sub_graph <- get_Gu(node, g)
  neighborhood_fastgreedy <- fastgreedy.community(sub_graph)
  membership <- neighborhood_fastgreedy$membership
  comSize <- sizes(neighborhood_fastgreedy)
  h_temp <- hash(keys=1,values=1)
  for(idx in 1:10)
  {
    members <- which(membership == idx)
    
    if(length(members) > 10)
    {
      hash:::.set(h_temp,keys=idx,values=comSize[idx])
    }
  }
  print(h_temp)
  hash:::.set(h_main,keys=node,values=h_temp)
}

for(key in keys(h_main))
{
  cat("Node = ", key, "\n")
  
  for(k in keys(values(h_main,keys=key)[[1]]))
  {
    val <- values(values(h_main,keys=key)[[1]],keys=k)
    if(k == 1 && val == 1)
      next
    if(val < 100)
      type=1
    else
      type=2
    cat("membership = ", k , " Size = " , val , " Community type = ",  type ,"\n")
  }
}

#for(node in core_nodes)
#{
#  sub_graph <- get_Gu(node, g)
#  neighborhood_fastgreedy <- fastgreedy.community(sub_graph)
#  membership <- neighborhood_fastgreedy$membership
#  print(table(membership))
#}

############  END OF QUESTION_6   #############



############ Question 7##############
library("igraph")

#####
##### google+ ego networks
#####



filesPath <- "D://Rtemp/gplus/"

allfiles <- list.files(filesPath)
IDs <-sub("^([^.]*).*", "\\1", allfiles)
egoNodeIds <- unique(IDs) 

for(ii in 1: length(egoNodeIds)){
  egoNodeId <- egoNodeIds[ii]
  
  edgelistFile <- paste(filesPath , egoNodeId  , ".edges" , sep="")
  circlesFile <- paste(filesPath , egoNodeId , ".circles" , sep="")
  
  g2Raw <- read.graph(edgelistFile , format = "ncol" , directed=TRUE)
  
  
  ### adding ego node in the ego network
  
  nonEgoNodes = V(g2Raw)
  
  g2 <- add.vertices(g2Raw,1,name=egoNodeId)
  egoNodeIndex <- which(V(g2)$name==egoNodeId)
  
  edgeAppendList <- c()
  for (nodeIndex in 1:(vcount(g2)-1)) {
    edgeAppendList <- c(edgeAppendList , c(vcount(g2),nodeIndex))
  }
  
  g2 <- add.edges(g2,edgeAppendList)
  
  g2U <- as.undirected(g2)
  
  ### reading circles
  
  fileConnection <- file(circlesFile , open="r")
  lines <- readLines(fileConnection)
  circles <- list()
  
  for (i in 1:length(lines)) {
    sp <- strsplit(lines[i],"\t")
    circles[[i]] <- sp[[1]][-1]
    
  }
  close(fileConnection)
  
  if(length(circles)>2) {
    
    wc <- walktrap.community(g2U)
    plot(wc,g2U,vertex.label=NA,vertex.size=7,edge.arrow.size=0.2,main="Community Structure for ego node 1 - WALK TRAP")
    
    # wc <- infomap.community(g2U)
    #plot(wc,g2U,vertex.label=NA,vertex.size=7,edge.arrow.size=0.2,main="Community Structure for ego node 1 - INFO MAP")
    
    for(j in 1:max(wc$membership)){
      select <- vector()
      for(i in 1:length(wc$membership)){
        if(wc$membership[i]==j){
          select <- c(select,(wc$name[i]))
        }
      }
      percentage <- vector()
      percentage_circle <- vector()
      for(k in 1:length(circles)){
        intersect_id <- intersect(select,circles[[k]])
        temp_percentage <- length(intersect_id)/length(select)
        temp_percentage_circle <- length(intersect_id)/length(circles[[k]])
        percentage_circle<- c(percentage_circle,temp_percentage_circle)
        percentage <- c(percentage, temp_percentage)
      }
      write(paste(egoNodeIds[ii]," :percentage:", percentage, " ,percentage_circle:", percentage_circle, " intersect_id:", intersect_id,sep=""), file="D://Rtemp/out4.txt", append=TRUE)
      write('=================', file="D://Rtemp/out4.txt", append=TRUE)
      
      print('-----------------------')
      print(percentage)
      print(percentage_circle)  
    }
    
  }
  
}


