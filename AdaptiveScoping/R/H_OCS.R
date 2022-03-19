## Calculate H-OCS score based on allier et al.
#-------------------------------------------------------------------------
#Copyright (c) 2019 Antoine Allier
#-------------------------------------------------------------------------

OCS_H<-function(candidate.marker.genos.i,
                        candidate.haplotype=candidate.haploid.i,
                        Bhat=predictions.out$solve.out$u,
                        pheno.train,
                        p.sel=100,
                        presel=300,
                        A.sc,
                        parents.information=parent.selections.list,
                        trainingpanel=list(geno=TP.genos.i,pheno=TP.phenos.i),
                        map=map,
                        rep.iter,
                        genome=hv.genome){
  
  Parent.genos <- genotype.loci(haploid.genos = candidate.haplotype, 
                                genome = genome, 
                                include.QTL = F)
  
  
  GEBV <-   Parent.genos %*% Bhat
  
  
  
  # Assuming the Elites are the 10 best lines
  PopE <- names(GEBV[order(GEBV,decreasing = TRUE),])[seq(p.sel/2)]
  
  preselect <- names(GEBV[order(GEBV,decreasing = TRUE),])[seq(presel)]

  PopD <- sample(setdiff(preselect,PopE),(length(preselect)-p.sel/2),replace = FALSE)
  
 
  
  ###########
  ###########
  Haplo<-as.matrix(do.call("rbind", candidate.haplotype))
              
  
  positions<-c()
  begin=0
  for (i in seq(7)){
    positions<-c(positions, genome[[i]]@pos.add.qtl$ID+begin)
    begin<-begin+genome[[i]]@num.snp.chr+genome[[i]]@num.add.qtl.chr
  }
  
  WindowSizeUI = 20
  StepSizeUI = 5

  
  Haplo<-Haplo[,-positions]
  ObjectHEBV <- GetHEBVmat(ObjectGeno = Haplo,
                           ObjectBeta = Bhat ,
                           ObjectMap = map,
                           WindowSize = WindowSizeUI,
                           StepSize = StepSizeUI)
  
  d1<-dim(ObjectHEBV$HEBV)
  HEBV<-c()
  for(i in seq(1,d1[1],by=2)){
    ind.hebv<-rbind(ObjectHEBV$HEBV[i,],ObjectHEBV$HEBV[i+1,])
    temp.hebv<-c()
    for (j in seq(d1[2])){
      temp.hebv<-c(temp.hebv,max(ind.hebv[,j]))
    }
    HEBV<-rbind(HEBV, temp.hebv)
  }
  
  line.names <- row.names(ObjectHEBV$HEBV) %>%
    str_replace(pattern = "\\.[0-9]$", replacement = "") %>%
    unique()
  
  rownames(HEBV)<-line.names
  colnames(HEBV)<- colnames(ObjectHEBV$HEBV)
  
  ObjectHEBV<-list(HEBV=HEBV,POSITION=ObjectHEBV$POSITION)
  ########
  ########
  
  SelDonor <- c()
  
  HSelDonor <- ComputeH(ObjectHEBVmat = ObjectHEBV$HEBV,# elites only
                         ListPopE = PopE,
                         ListPopD = NULL)
  
  invisible(lapply(1:(p.sel/2),function(iter){
    tmp <- ComputeH(ObjectHEBVmat = ObjectHEBV$HEBV,
                    ListPopE = c(SelDonor,PopE),# PopE + selected donors
                    ListPopD = PopD[!PopD%in%c(SelDonor)] # yet unselected donors
    )
    #Increment with selected donor
    SelDonor <<- c(SelDonor,tmp[order(tmp$H,decreasing =TRUE),]$LINE[1])
    # Increment with H criterion value of elites + selected donors
    HSelDonor <<- rbind(HSelDonor,tmp[order(tmp$H,decreasing =TRUE),][1,])
  }))
  
  # Assuming we have selected the donors:
  SelDonors4UC <- HSelDonor$LINE[seq(2,p.sel/2+1)]

      
      UCtmp <-getGaSolutionsFrontier(
        Markers=Parent.genos[PopE,],
        Markers2=Parent.genos[SelDonors4UC,],
        K=A.sc[c(PopE,SelDonors4UC),c(PopE,SelDonors4UC)],
        markereffects=Bhat,
        markermap=NULL,
        nmates=p.sel/2,
        npopGA=100,
        nitGA=100,
        mc.cores=1,
        mutprob=0.999,
        noself=TRUE,
        method=3,
        type=2L,
        generation=1L,
        plotiters=F)
      
      n.iters=length(UCtmp[[2]])
      
  cros<-UCtmp[[2]][[n.iters]]
  colnames(cros)<-c("Parent1","Parent2")
  
  selection<-c(cros[,1], cros[,2] )


  parent.sel<-c()
  parent.sel$lines.sel<-unique(selection)
  
  parent.sel$value.sel<-GEBV[unique(selection),]
  
  return(list(parent.selections.i=parent.sel,crossing.block.i=cros))
  
}
