compute_taxa_milkType_association <- function(asfbf_fdrFile=fdrCI_lefse_phylum_level_6m_asfbf, tfbf_fdrFile = fdrCI_lefse_phylum_level_6m_tfbf, taxaFile = phylumTaxa)
{
  dataSet <-thisDataInstance_all
  dataSetInf1m  <- dataSet[grepl("Inf.Bl",dataSet$X.SampleID),]
  dataSetInf6m  <- dataSet[grepl("Inf.6m",dataSet$X.SampleID),]

  fdr_asfbf <- asfbf_fdrFile
  fdr_tfbf <- tfbf_fdrFile

  taxaTable <- taxaFile
  taxaTable <- as.data.frame(taxaTable)
  taxaTable$taxa <- gsub("\\| ",".",taxaTable$taxa)
  row.names(taxaTable) <- taxaTable$taxa
  taxaTable_6m<- taxaTable[,grepl("Inf.6m",colnames(taxaTable))]
  myT <- t(taxaTable_6m)
  myT2 <- (myT/(rowSums(myT)+1))*(sum(rowSums(myT))/nrow(myT))
  myTLogged <- log10((myT2 +1))
  myTLogged_6m<- myTLogged[,(colSums(myTLogged==0)/nrow(myTLogged)) <= 0.9]

  taxaTable_1m <- taxaTable[,grepl("Inf.Bl",colnames(taxaTable))]
  myT_1m <- t(taxaTable_1m)
  myT2_1m <- ( myT_1m/(rowSums(myT_1m)+1))*(sum(rowSums(myT_1m))/nrow(myT_1m))
  myTLogged_1m <- log10((myT2_1m +1))

  taxaMeta_6m <- merge(dataSetInf6m,myTLogged_6m,by.x="X.SampleID",by.y="row.names")
  taxaMeta_1m <- merge(dataSetInf1m,myTLogged_1m,by.x="X.SampleID",by.y="row.names")

  exposureNamesList <- character(0);allOutcomes <- character(0);
  meanBF_list <- numeric(0);meanTF_list <- numeric(0);meanASF_list <- numeric(0);
  sdBF_list <- numeric(0);sdTF_list <- numeric(0);sdASF_list <- numeric(0);
  diffTFBF_list <- numeric(0);diffASFBF_list <- numeric(0);
  lwrTFBF_list <- numeric(0);lwrASFBF_list <- numeric(0);
  uprTFBF_list <- numeric(0);uprASFBF_list <- numeric(0);
  pValTFBF_list <- numeric(0);pValASFBF_list <- numeric(0);
  tfbf_fdrci_list <- character(0);asfbf_fdrci_list <- character(0);

  for(i in (length(dataSetInf6m)+1):(length(taxaMeta_6m)))
  {
    allOutcomes[[length(allOutcomes)+1]] <- names(taxaMeta_6m)[i];

    modelForm=as.formula(paste("thisMicrobe","~milkConsumptionType+mother_age+prepreg_bmi_kgm2 +mode_of_delivery+mom_BMI+sex+ baby_age+inf_weight_kg+fruit_incj_Inf+thisMicrobe_1m"  ));

    thisDataInstance_6m <- na.omit(data.frame(dyad_id=taxaMeta_6m$dyad_id,
                                              thisMicrobe=taxaMeta_6m[,i],
                                              milkConsumptionType=taxaMeta_6m$milkConsumptionType,
                                              mother_age=taxaMeta_6m$mother_age,
                                              prepreg_bmi_kgm2=taxaMeta_6m$prepreg_bmi_kgm2,
                                              mode_of_delivery=taxaMeta_6m$mode_of_delivery,
                                              mom_BMI=taxaMeta_6m$mom_current_BMI,
                                              sex=taxaMeta_6m$sex,
                                              baby_age=taxaMeta_6m$baby_age,
                                              inf_weight_kg=taxaMeta_6m$inf_weight_kg_1m6m,
                                              fruit_incj_Inf=taxaMeta_6m$fruit_incj_Inf_1m6m))

    thisDataInstance_1m <- na.omit(data.frame(dyad_id=taxaMeta_1m$dyad_id,
                                              thisMicrobe_1m=taxaMeta_1m[,colnames(taxaMeta_1m) %in% names(taxaMeta_6m[i])]))
    thisDataInstance_1m6m <- merge(thisDataInstance_6m,thisDataInstance_1m,by="dyad_id")
    thisDataInstance_1m6m$milkConsumptionType <- factor(thisDataInstance_1m6m$milkConsumptionType)

    modelInfo <- lm(modelForm, data = thisDataInstance_1m6m)
    statModelAnova <- aov(modelInfo)
    tukey <- suppressWarnings( TukeyHSD(statModelAnova,"milkConsumptionType"))

    meanBF <- format(mean(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "BF",names(thisDataInstance_1m6m) %in% "thisMicrobe"]),digits=2)
    meanBF_list[[length(meanBF_list)+1]] <- meanBF;

    meanTF <- format(mean(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "TF",names(thisDataInstance_1m6m) %in% "thisMicrobe"]),digits=2)
    meanTF_list[[length(meanTF_list)+1]] <- meanTF;

    meanASF <- format(mean(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "ASF",names(thisDataInstance_1m6m) %in% "thisMicrobe"]),digits=2)
    meanASF_list[[length(meanASF_list)+1]] <- meanASF;

    sdBF <- format(sd(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "BF",names(thisDataInstance_1m6m) %in% "thisMicrobe"]),digits=2)
    sdBF_list[[length(sdBF_list)+1]] <- sdBF;

    sdTF <- format(sd(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "TF",names(thisDataInstance_1m6m) %in% "thisMicrobe"]),digits=2)
    sdTF_list[[length(sdTF_list)+1]] <- sdTF;

    sdASF <- format(sd(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "ASF",names(thisDataInstance_1m6m) %in% "thisMicrobe"]),digits=2)
    sdASF_list[[length(sdASF_list)+1]] <- sdASF;


    pValTFBF <- format.pval(tukey$milkConsumptionType[1,4],digits=2);
    pValTFBF_list[[length(pValTFBF_list)+1]] <- pValTFBF;
    # fdr_tfbf_value <- paste(format.pval(fdr_tfbf[fdr_tfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[1,4])*10^3)/10^3),][2],digits=2)," (",
    #                         format.pval(fdr_tfbf[fdr_tfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[1,4])*10^3)/10^3),][3],digits=2),",",
    #                         format.pval(fdr_tfbf[fdr_tfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[1,4])*10^3)/10^3),][4],digits=2),")",sep="")
    # tfbf_fdrci_list[[length(tfbf_fdrci_list)+1]] <- fdr_tfbf_value;

    pValASFBF <- format.pval(tukey$milkConsumptionType[2,4],digits=2);
    pValASFBF_list[[length(pValASFBF_list)+1]] <- pValASFBF;
    # fdr_asfbf_value <- paste(format.pval(fdr_asfbf[fdr_asfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[2,4])*10^3)/10^3),][2],digits=2)," (",
    #                          format.pval(fdr_asfbf[fdr_asfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[2,4])*10^3)/10^3),][3],digits=2),",",
    #                          format.pval(fdr_asfbf[fdr_asfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[2,4])*10^3)/10^3),][4],digits=2),")",sep="")
    # asfbf_fdrci_list[[length(asfbf_fdrci_list)+1]] <- fdr_asfbf_value;

    diffTFBF <- format(tukey$milkConsumptionType[1,1],digits=2);
    diffTFBF_list[[length(diffTFBF_list)+1]] <- diffTFBF;

    diffASFBF <- format(tukey$milkConsumptionType[2,1],digits=2);
    diffASFBF_list[[length(diffASFBF_list)+1]] <- diffASFBF;

    lwrTFBF <- format(tukey$milkConsumptionType[1,2],digits=2);
    lwrTFBF_list[[length(lwrTFBF_list)+1]] <- lwrTFBF;

    lwrASFBF <- format(tukey$milkConsumptionType[2,2],digits=2);
    lwrASFBF_list[[length(lwrASFBF_list)+1]] <- lwrASFBF;

    uprTFBF <- format(tukey$milkConsumptionType[1,3],digits=2);
    uprTFBF_list[[length(uprTFBF_list)+1]] <- uprTFBF;

    uprASFBF <- format(tukey$milkConsumptionType[2,3],digits=2);
    uprASFBF_list[[length(uprASFBF_list)+1]] <- uprASFBF;

  }

  bfList <- paste(meanBF_list," (",sdBF_list,")",sep="")
  tfList <- paste(meanTF_list," (",sdTF_list,")",sep="")
  asfList <- paste(meanASF_list," (",sdASF_list,")",sep="")

  difftfbfList <- paste(diffTFBF_list," (",lwrTFBF_list,", ",uprTFBF_list,")",sep="")
  diffasfbfList <- paste(diffASFBF_list," (",lwrASFBF_list,", ",uprASFBF_list,")",sep="")

  # allOutcomes <- cbind(allOutcomes,
  #                       bfList,tfList,asfList,
  #                       difftfbfList,pValTFBF_list,tfbf_fdrci_list,
  #                       diffasfbfList,pValASFBF_list,asfbf_fdrci_list);

  allOutcomes <- cbind(allOutcomes,
                        bfList,tfList,asfList,
                        difftfbfList,pValTFBF_list,
                        diffasfbfList,pValASFBF_list);

  allOutcomes <- data.frame(allOutcomes)
  return(allOutcomes)
}

compute_alphaDiversity_milkType_association <- function()
{
  dataSet <-thisDataInstance_all
  dataSetInf1m  <- dataSet[grepl("Inf.Bl",dataSet$X.SampleID),]
  dataSetInf6m  <- dataSet[grepl("Inf.6m",dataSet$X.SampleID),]

  fdr_asfbf <- alphaFDRCI_asfbf
  fdr_tfbf <- alphaFDRCI_tfbf

  taxaTable_6m<- alphaDiversityTable[grepl("Inf.6m",rownames(alphaDiversityTable)),]
  taxaTable_1m <- alphaDiversityTable[grepl("Inf.Bl",rownames(alphaDiversityTable)),]

  taxaMeta_6m <- merge(dataSetInf6m,taxaTable_6m,by.x="X.SampleID",by.y="row.names")
  taxaMeta_1m <- merge(dataSetInf1m,taxaTable_1m,by.x="X.SampleID",by.y="row.names")

  exposureNamesList <- character(0);allOutcomes <- character(0);
  meanBF_list <- numeric(0);meanTF_list <- numeric(0);meanASF_list <- numeric(0);
  sdBF_list <- numeric(0);sdTF_list <- numeric(0);sdASF_list <- numeric(0);
  diffTFBF_list <- numeric(0);diffASFBF_list <- numeric(0);
  lwrTFBF_list <- numeric(0);lwrASFBF_list <- numeric(0);
  uprTFBF_list <- numeric(0);uprASFBF_list <- numeric(0);
  pValTFBF_list <- numeric(0);pValASFBF_list <- numeric(0);
  tfbf_fdrci_list <- character(0);asfbf_fdrci_list <- character(0);

  for(i in (length(dataSetInf6m)+1):(length(taxaMeta_6m)))
  {

    allOutcomes[[length(allOutcomes)+1]] <- names(taxaMeta_6m)[i];

    modelForm=as.formula(paste("thisMicrobe","~milkConsumptionType+mother_age+prepreg_bmi_kgm2 +mode_of_delivery+mom_BMI+sex+ baby_age+inf_weight_kg+fruit_incj_Inf+thisMicrobe_1m"  ));

    thisDataInstance_6m <- na.omit(data.frame(dyad_id=taxaMeta_6m$dyad_id,
                                              thisMicrobe=taxaMeta_6m[,i],
                                              milkConsumptionType=taxaMeta_6m$milkConsumptionType,
                                              mother_age=taxaMeta_6m$mother_age,
                                              prepreg_bmi_kgm2=taxaMeta_6m$prepreg_bmi_kgm2,
                                              mode_of_delivery=taxaMeta_6m$mode_of_delivery,
                                              mom_BMI=taxaMeta_6m$mom_current_BMI,
                                              sex=taxaMeta_6m$sex,
                                              baby_age=taxaMeta_6m$baby_age,
                                              inf_weight_kg=taxaMeta_6m$inf_weight_kg_1m6m,
                                              fruit_incj_Inf=taxaMeta_6m$fruit_incj_Inf_1m6m))

    thisDataInstance_1m <- na.omit(data.frame(dyad_id=taxaMeta_1m$dyad_id,
                                              thisMicrobe_1m=taxaMeta_1m[,colnames(taxaMeta_1m) %in% names(taxaMeta_6m[i])]))
    thisDataInstance_1m6m <- merge(thisDataInstance_6m,thisDataInstance_1m,by="dyad_id")
    thisDataInstance_1m6m$milkConsumptionType <- factor(thisDataInstance_1m6m$milkConsumptionType)

    modelInfo <- lm(modelForm, data = thisDataInstance_1m6m)
    statModelAnova <- aov(modelInfo)
    tukey <- suppressWarnings( TukeyHSD(statModelAnova,"milkConsumptionType"))

    meanBF <- format(mean(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "BF",names(thisDataInstance_1m6m) %in% "thisMicrobe"]),digits=2)
    meanBF_list[[length(meanBF_list)+1]] <- meanBF;

    meanTF <- format(mean(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "TF",names(thisDataInstance_1m6m) %in% "thisMicrobe"]),digits=2)
    meanTF_list[[length(meanTF_list)+1]] <- meanTF;

    meanASF <- format(mean(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "ASF",names(thisDataInstance_1m6m) %in% "thisMicrobe"]),digits=2)
    meanASF_list[[length(meanASF_list)+1]] <- meanASF;

    sdBF <- format(sd(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "BF",names(thisDataInstance_1m6m) %in% "thisMicrobe"]),digits=2)
    sdBF_list[[length(sdBF_list)+1]] <- sdBF;

    sdTF <- format(sd(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "TF",names(thisDataInstance_1m6m) %in% "thisMicrobe"]),digits=2)
    sdTF_list[[length(sdTF_list)+1]] <- sdTF;

    sdASF <- format(sd(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "ASF",names(thisDataInstance_1m6m) %in% "thisMicrobe"]),digits=2)
    sdASF_list[[length(sdASF_list)+1]] <- sdASF;


    pValTFBF <- format.pval(tukey$milkConsumptionType[1,4],digits=2);
    pValTFBF_list[[length(pValTFBF_list)+1]] <- pValTFBF;
    # fdr_tfbf_value <- paste(format.pval(fdr_tfbf[fdr_tfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[1,4])*10^3)/10^3),][2],digits=2)," (",
    #                         format.pval(fdr_tfbf[fdr_tfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[1,4])*10^3)/10^3),][3],digits=2),",",
    #                         format.pval(fdr_tfbf[fdr_tfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[1,4])*10^3)/10^3),][4],digits=2),")",sep="")
    # tfbf_fdrci_list[[length(tfbf_fdrci_list)+1]] <- fdr_tfbf_value;

    pValASFBF <- format.pval(tukey$milkConsumptionType[2,4],digits=2);
    pValASFBF_list[[length(pValASFBF_list)+1]] <- pValASFBF;
    # fdr_asfbf_value <- paste(format.pval(fdr_asfbf[fdr_asfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[2,4])*10^3)/10^3),][2],digits=2)," (",
    #                          format.pval(fdr_asfbf[fdr_asfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[2,4])*10^3)/10^3),][3],digits=2),",",
    #                          format.pval(fdr_asfbf[fdr_asfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[2,4])*10^3)/10^3),][4],digits=2),")",sep="")
    # asfbf_fdrci_list[[length(asfbf_fdrci_list)+1]] <- fdr_asfbf_value;

    diffTFBF <- format(tukey$milkConsumptionType[1,1],digits=2);
    diffTFBF_list[[length(diffTFBF_list)+1]] <- diffTFBF;

    diffASFBF <- format(tukey$milkConsumptionType[2,1],digits=2);
    diffASFBF_list[[length(diffASFBF_list)+1]] <- diffASFBF;


    lwrTFBF <- format(tukey$milkConsumptionType[1,2],digits=2);
    lwrTFBF_list[[length(lwrTFBF_list)+1]] <- lwrTFBF;

    lwrASFBF <- format(tukey$milkConsumptionType[2,2],digits=2);
    lwrASFBF_list[[length(lwrASFBF_list)+1]] <- lwrASFBF;

    uprTFBF <- format(tukey$milkConsumptionType[1,3],digits=2);
    uprTFBF_list[[length(uprTFBF_list)+1]] <- uprTFBF;

    uprASFBF <- format(tukey$milkConsumptionType[2,3],digits=2);
    uprASFBF_list[[length(uprASFBF_list)+1]] <- uprASFBF;

  }

  bfList <- paste(meanBF_list," (",sdBF_list,")",sep="")
  tfList <- paste(meanTF_list," (",sdTF_list,")",sep="")
  asfList <- paste(meanASF_list," (",sdASF_list,")",sep="")

  difftfbfList <- paste(diffTFBF_list," (",lwrTFBF_list,", ",uprTFBF_list,")",sep="")
  diffasfbfList <- paste(diffASFBF_list," (",lwrASFBF_list,", ",uprASFBF_list,")",sep="")

  # allOutcomes <- cbind(allOutcomes,
  #                       bfList,tfList,asfList,
  #                       difftfbfList,pValTFBF_list,tfbf_fdrci_list,
  #                       diffasfbfList,pValASFBF_list,asfbf_fdrci_list);

  allOutcomes <- cbind(allOutcomes,
                        bfList,tfList,asfList,
                        difftfbfList,pValTFBF_list,tfbf_fdrci_list,
                        diffasfbfList,pValASFBF_list);

  allOutcomes <- data.frame(allOutcomes)
  return(allOutcomes)
}

compute_betaDiversity_milkType_association <- function()
{
  dataSet <-thisDataInstance_all
  dataSetInf1m  <- dataSet[grepl("Inf.Bl",dataSet$X.SampleID),]
  dataSetInf6m  <- dataSet[grepl("Inf.6m",dataSet$X.SampleID),]

  fdr_asfbf <- betaFDRCI_asfbf
  fdr_tfbf <- betaFDRCI_tfbf

  loadings6m <- eigen6m

  taxaMeta_1m <- merge(dataSetInf1m,ordinationData1m,by.x="X.SampleID",by.y="row.names")
  taxaMeta_6m <- merge(dataSetInf6m,ordinationData6m,by.x="X.SampleID",by.y="row.names")

  exposureNamesList <- character(0);allOutcomes <- character(0);
  meanBF_list <- numeric(0);meanTF_list <- numeric(0);meanASF_list <- numeric(0);
  sdBF_list <- numeric(0); sdTF_list <- numeric(0);sdASF_list <- numeric(0);
  diffTFBF_list <- numeric(0);diffASFBF_list <- numeric(0);
  lwrTFBF_list <- numeric(0);lwrASFBF_list <- numeric(0);
  uprTFBF_list <- numeric(0);uprASFBF_list <- numeric(0);
  pValTFBF_list <- numeric(0);pValASFBF_list <- numeric(0);
  tfbf_fdrci_list <- character(0);asfbf_fdrci_list <- character(0);

  for(i in (length(dataSetInf6m)+1):(length(taxaMeta_6m)))
  {
    allOutcomes[[length(allOutcomes)+1]] <- names(taxaMeta_6m)[i];

    modelForm=as.formula(paste("thisAxis","~milkConsumptionType+mother_age+prepreg_bmi_kgm2 +mode_of_delivery+mom_BMI+sex+ baby_age+inf_weight_kg+fruit_incj_Inf+thisAxis_1m"  ));

    thisDataInstance_6m <- na.omit(data.frame(dyad_id=taxaMeta_6m$dyad_id,
                                              thisAxis=taxaMeta_6m[,i],
                                              milkConsumptionType=taxaMeta_6m$milkConsumptionType,
                                              mother_age=taxaMeta_6m$mother_age,
                                              prepreg_bmi_kgm2=taxaMeta_6m$prepreg_bmi_kgm2,
                                              mode_of_delivery=taxaMeta_6m$mode_of_delivery,
                                              mom_BMI=taxaMeta_6m$mom_current_BMI,
                                              sex=taxaMeta_6m$sex,
                                              baby_age=taxaMeta_6m$baby_age,
                                              inf_weight_kg=taxaMeta_6m$inf_weight_kg_1m6m,
                                              fruit_incj_Inf=taxaMeta_6m$fruit_incj_Inf_1m6m))

    thisDataInstance_1m <- na.omit(data.frame(dyad_id=taxaMeta_1m$dyad_id,
                                              thisAxis_1m=taxaMeta_1m[,colnames(taxaMeta_1m) %in% names(taxaMeta_6m[i])]))
    thisDataInstance_1m6m <- merge(thisDataInstance_6m,thisDataInstance_1m,by="dyad_id")
    thisDataInstance_1m6m$milkConsumptionType <- factor(thisDataInstance_1m6m$milkConsumptionType)

    modelInfo <- lm(modelForm, data = thisDataInstance_1m6m)
    statModelAnova <- aov(modelInfo)
    tukey <- suppressWarnings( TukeyHSD(statModelAnova,"milkConsumptionType"))

    meanBF <- format(mean(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "BF",names(thisDataInstance_1m6m) %in% "thisAxis"]),digits=2)
    meanBF_list[[length(meanBF_list)+1]] <- meanBF;

    meanTF <- format(mean(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "TF",names(thisDataInstance_1m6m) %in% "thisAxis"]),digits=2)
    meanTF_list[[length(meanTF_list)+1]] <- meanTF;

    meanASF <- format(mean(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "ASF",names(thisDataInstance_1m6m) %in% "thisAxis"]),digits=2)
    meanASF_list[[length(meanASF_list)+1]] <- meanASF;

    sdBF <- format(sd(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "BF",names(thisDataInstance_1m6m) %in% "thisAxis"]),digits=2)
    sdBF_list[[length(sdBF_list)+1]] <- sdBF;

    sdTF <- format(sd(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "TF",names(thisDataInstance_1m6m) %in% "thisAxis"]),digits=2)
    sdTF_list[[length(sdTF_list)+1]] <- sdTF;

    sdASF <- format(sd(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "ASF",names(thisDataInstance_1m6m) %in% "thisAxis"]),digits=2)
    sdASF_list[[length(sdASF_list)+1]] <- sdASF;


    pValTFBF <- format.pval(tukey$milkConsumptionType[1,4],digits=2);
    pValTFBF_list[[length(pValTFBF_list)+1]] <- pValTFBF;
    # fdr_tfbf_value <- paste(format.pval(fdr_tfbf[fdr_tfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[1,4])*10^3)/10^3),][2],digits=2)," (",
    #                         format.pval(fdr_tfbf[fdr_tfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[1,4])*10^3)/10^3),][3],digits=2),",",
    #                         format.pval(fdr_tfbf[fdr_tfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[1,4])*10^3)/10^3),][4],digits=2),")",sep="")
    # tfbf_fdrci_list[[length(tfbf_fdrci_list)+1]] <- fdr_tfbf_value;

    pValASFBF <- format.pval(tukey$milkConsumptionType[2,4],digits=2);
    pValASFBF_list[[length(pValASFBF_list)+1]] <- pValASFBF;
    # fdr_asfbf_value <- paste(format.pval(fdr_asfbf[fdr_asfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[2,4])*10^3)/10^3),][2],digits=2)," (",
    #                          format.pval(fdr_asfbf[fdr_asfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[2,4])*10^3)/10^3),][3],digits=2),",",
    #                          format.pval(fdr_asfbf[fdr_asfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[2,4])*10^3)/10^3),][4],digits=2),")",sep="")
    # asfbf_fdrci_list[[length(asfbf_fdrci_list)+1]] <- fdr_asfbf_value;

    diffTFBF <- format(tukey$milkConsumptionType[1,1],digits=2);
    diffTFBF_list[[length(diffTFBF_list)+1]] <- diffTFBF;

    diffASFBF <- format(tukey$milkConsumptionType[2,1],digits=2);
    diffASFBF_list[[length(diffASFBF_list)+1]] <- diffASFBF;


    lwrTFBF <- format(tukey$milkConsumptionType[1,2],digits=2);
    lwrTFBF_list[[length(lwrTFBF_list)+1]] <- lwrTFBF;

    lwrASFBF <- format(tukey$milkConsumptionType[2,2],digits=2);
    lwrASFBF_list[[length(lwrASFBF_list)+1]] <- lwrASFBF;

    uprTFBF <- format(tukey$milkConsumptionType[1,3],digits=2);
    uprTFBF_list[[length(uprTFBF_list)+1]] <- uprTFBF;

    uprASFBF <- format(tukey$milkConsumptionType[2,3],digits=2);
    uprASFBF_list[[length(uprASFBF_list)+1]] <- uprASFBF;

  }

  bfList <- paste(meanBF_list," (",sdBF_list,")",sep="")
  tfList <- paste(meanTF_list," (",sdTF_list,")",sep="")
  asfList <- paste(meanASF_list," (",sdASF_list,")",sep="")

  difftfbfList <- paste(diffTFBF_list," (",lwrTFBF_list,", ",uprTFBF_list,")",sep="")
  diffasfbfList <- paste(diffASFBF_list," (",lwrASFBF_list,", ",uprASFBF_list,")",sep="")

  # allOutcomes <- cbind(allOutcomes,
  #                       bfList,tfList,asfList,
  #                       difftfbfList,pValTFBF_list,tfbf_fdrci_list,
  #                       diffasfbfList,pValASFBF_list,asfbf_fdrci_list);
  allOutcomes <- cbind(allOutcomes,
                        bfList,tfList,asfList,
                        difftfbfList,pValTFBF_list,
                        diffasfbfList,pValASFBF_list);

  allOutcomes <- data.frame(allOutcomes)
  return(allOutcomes)
}

compute_keggModule_milkType_association <- function()
{
  dataSet <-thisDataInstance_all
  dataSetInf1m  <- dataSet[grepl("Inf.Bl",dataSet$X.SampleID),]
  dataSetInf6m  <- dataSet[grepl("Inf.6m",dataSet$X.SampleID),]

  fdr_asfbf <- keggModuleFDRCI_asfbf
  fdr_tfbf <- keggModuleFDRCI_tfbf

  moduleTable <- humann_modulec_picrust
  moduleTable <- moduleTable[-c(1,2),]
  colnames(moduleTable) <- substr(colnames(moduleTable),1,nchar(colnames(moduleTable))-10)
  moduleTableBaby <- moduleTable[,grepl("Inf.6m",colnames(moduleTable))]
  myT <- t(moduleTableBaby)
  myTLogged_6m <- log10((myT/(rowSums(myT)+1))*(sum(rowSums(myT))/nrow(myT)) +1)
  myTLogged_6m <- myTLogged_6m[,(colSums(myTLogged_6m==0)/nrow(myTLogged_6m)) <= 0.9]

  moduleTableBaby_1m <- moduleTable[,grepl("Inf.Bl",colnames(moduleTable))]
  myT_1m <- t(moduleTableBaby_1m)
  myTLogged_1m<- log10((myT_1m/(rowSums(myT_1m)+1))*(sum(rowSums(myT_1m))/nrow(myT_1m)) +1)

  moduleMeta_6m <- merge(dataSetInf6m,myTLogged_6m,by.x="X.SampleID",by.y="row.names")
  moduleMeta_1m <- merge(dataSetInf1m,myTLogged_1m,by.x="X.SampleID",by.y="row.names")

  exposureNamesList <- character(0);allOutcomes <- character(0);
  meanBF_list <- numeric(0);meanTF_list <- numeric(0);meanASF_list <- numeric(0);
  sdBF_list <- numeric(0);sdTF_list <- numeric(0);sdASF_list <- numeric(0);
  diffTFBF_list <- numeric(0);diffASFBF_list <- numeric(0);
  lwrTFBF_list <- numeric(0);lwrASFBF_list <- numeric(0);
  uprTFBF_list <- numeric(0);uprASFBF_list <- numeric(0);
  pValTFBF_list <- numeric(0);pValASFBF_list <- numeric(0);
  tfbf_fdrci_list <- character(0);asfbf_fdrci_list <- character(0);

  for(i in (length(dataSetInf6m)+1):(length(moduleMeta_6m)))
  {
    allOutcomes[[length(allOutcomes)+1]] <- names(moduleMeta_6m)[i];

    modelForm=as.formula(paste("thisModule","~milkConsumptionType+mother_age+prepreg_bmi_kgm2 +mode_of_delivery+mom_BMI+sex+ baby_age+inf_weight_kg+fruit_incj_Inf+thisModule_1m"  ));

    thisDataInstance_6m <- na.omit(data.frame(dyad_id=moduleMeta_6m$dyad_id,
                                              thisModule=moduleMeta_6m[,i],
                                              milkConsumptionType=moduleMeta_6m$milkConsumptionType,
                                              mother_age=moduleMeta_6m$mother_age,
                                              prepreg_bmi_kgm2=moduleMeta_6m$prepreg_bmi_kgm2,
                                              mode_of_delivery=moduleMeta_6m$mode_of_delivery,
                                              mom_BMI=moduleMeta_6m$mom_current_BMI,
                                              sex=moduleMeta_6m$sex,
                                              baby_age=moduleMeta_6m$baby_age,
                                              inf_weight_kg=moduleMeta_6m$inf_weight_kg_1m6m,
                                              fruit_incj_Inf=moduleMeta_6m$fruit_incj_Inf_1m6m))

    thisDataInstance_1m <- na.omit(data.frame(dyad_id=moduleMeta_1m$dyad_id,
                                              thisModule_1m=moduleMeta_1m[,i]))
    thisDataInstance_1m6m <- merge(thisDataInstance_6m,thisDataInstance_1m,by="dyad_id")
    thisDataInstance_1m6m$milkConsumptionType <- factor(thisDataInstance_1m6m$milkConsumptionType)


    modelInfo <- lm(modelForm, data = thisDataInstance_1m6m)
    statModelAnova <- aov(modelInfo)
    tukey <- suppressWarnings( TukeyHSD(statModelAnova,"milkConsumptionType"))

    meanBF <- format(mean(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "BF",names(thisDataInstance_1m6m) %in% "thisModule"]),digits=2)
    meanBF_list[[length(meanBF_list)+1]] <- meanBF;

    meanTF <- format(mean(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "TF",names(thisDataInstance_1m6m) %in% "thisModule"]),digits=2)
    meanTF_list[[length(meanTF_list)+1]] <- meanTF;

    meanASF <- format(mean(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "ASF",names(thisDataInstance_1m6m) %in% "thisModule"]),digits=2)
    meanASF_list[[length(meanASF_list)+1]] <- meanASF;

    sdBF <- format(sd(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "BF",names(thisDataInstance_1m6m) %in% "thisModule"]),digits=2)
    sdBF_list[[length(sdBF_list)+1]] <- sdBF;

    sdTF <- format(sd(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "TF",names(thisDataInstance_1m6m) %in% "thisModule"]),digits=2)
    sdTF_list[[length(sdTF_list)+1]] <- sdTF;

    sdASF <- format(sd(thisDataInstance_1m6m[thisDataInstance_1m6m$milkConsumptionType %in% "ASF",names(thisDataInstance_1m6m) %in% "thisModule"]),digits=2)
    sdASF_list[[length(sdASF_list)+1]] <- sdASF;


    pValTFBF <- format.pval(tukey$milkConsumptionType[1,4],digits=2);
    pValTFBF_list[[length(pValTFBF_list)+1]] <- pValTFBF;
    # fdr_tfbf_value <- paste(format.pval(fdr_tfbf[fdr_tfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[1,4])*10^3)/10^3),][2],digits=2)," (",
    #                         format.pval(fdr_tfbf[fdr_tfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[1,4])*10^3)/10^3),][3],digits=2),",",
    #                         format.pval(fdr_tfbf[fdr_tfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[1,4])*10^3)/10^3),][4],digits=2),")",sep="")
    # tfbf_fdrci_list[[length(tfbf_fdrci_list)+1]] <- fdr_tfbf_value;

    pValASFBF <- format.pval(tukey$milkConsumptionType[2,4],digits=2);
    pValASFBF_list[[length(pValASFBF_list)+1]] <- pValASFBF;
    # fdr_asfbf_value <- paste(format.pval(fdr_asfbf[fdr_asfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[2,4])*10^3)/10^3),][2],digits=2)," (",
    #                          format.pval(fdr_asfbf[fdr_asfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[2,4])*10^3)/10^3),][3],digits=2),",",
    #                          format.pval(fdr_asfbf[fdr_asfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[2,4])*10^3)/10^3),][4],digits=2),")",sep="")
    # asfbf_fdrci_list[[length(asfbf_fdrci_list)+1]] <- fdr_asfbf_value;

    diffTFBF <- format(tukey$milkConsumptionType[1,1],digits=2);
    diffTFBF_list[[length(diffTFBF_list)+1]] <- diffTFBF;

    diffASFBF <- format(tukey$milkConsumptionType[2,1],digits=2);
    diffASFBF_list[[length(diffASFBF_list)+1]] <- diffASFBF;

    lwrTFBF <- format(tukey$milkConsumptionType[1,2],digits=2);
    lwrTFBF_list[[length(lwrTFBF_list)+1]] <- lwrTFBF;

    lwrASFBF <- format(tukey$milkConsumptionType[2,2],digits=2);
    lwrASFBF_list[[length(lwrASFBF_list)+1]] <- lwrASFBF;

    uprTFBF <- format(tukey$milkConsumptionType[1,3],digits=2);
    uprTFBF_list[[length(uprTFBF_list)+1]] <- uprTFBF;

    uprASFBF <- format(tukey$milkConsumptionType[2,3],digits=2);
    uprASFBF_list[[length(uprASFBF_list)+1]] <- uprASFBF;
  }

  bfList <- paste(meanBF_list," (",sdBF_list,")",sep="")
  tfList <- paste(meanTF_list," (",sdTF_list,")",sep="")
  asfList <- paste(meanASF_list," (",sdASF_list,")",sep="")

  difftfbfList <- paste(diffTFBF_list," (",lwrTFBF_list,", ",uprTFBF_list,")",sep="")
  diffasfbfList <- paste(diffASFBF_list," (",lwrASFBF_list,", ",uprASFBF_list,")",sep="")

  # allOutcomes <- cbind(allOutcomes,
  #                       bfList,tfList,asfList,
  #                       difftfbfList,pValTFBF_list,tfbf_fdrci_list,
  #                       diffasfbfList,pValASFBF_list,asfbf_fdrci_list);
  allOutcomes <- cbind(allOutcomes,
                        bfList,tfList,asfList,
                        difftfbfList,pValTFBF_list,
                        diffasfbfList,pValASFBF_list);

  allOutcomes <- data.frame(allOutcomes)
  return(allOutcomes)
}

compute_possibleCovariates_milkType <- function()
{
  dataSet <-thisDataInstance_all
  dataSetInf6m  <- dataSet[grepl("Inf.6m",dataSet$X.SampleID),]

  fdr_sugar <- potentialCovariates

  potentialCovariates <- c("sex","baby_age","inf_weight_kg_1m6m","inf_length_cm_6m","zbmi_6m","zwfl_6m","zwei_6m","zlen_6m","skinf_tricep_mm_6m",
                           "skinf_subscap_mm_6m","skinf_supra_mm_6m","skinf_midthigh_mm_6m","circ_umb_cm_6m","mother_age","mom_current_BMI","prepreg_bmi_kgm2",
                           "mode_of_delivery","begun_solid_food","m_ener_Inf_6m","fruit_incj_Inf_1m6m","PAS","NEG","ORC")

  exposureNamesList <- character(0);allOutcomes <- character(0);

  meanBF_list <- numeric(0);meanTF_list <- numeric(0);meanASF_list <- numeric(0);
  sdBF_list <- numeric(0);sdTF_list <- numeric(0);sdASF_list <- numeric(0);
  sugar_fdrci_list <- numeric(0); pValSugarList <- numeric(0);

  i <- 1
  for( variable in potentialCovariates)
  {
    allOutcomes[[length(allOutcomes)+1]] <- variable;
    categorical <- ifelse(variable %in% c("sex","mode_of_delivery","begun_solid_food"),TRUE,FALSE)

    thisDataInstance <- na.omit(data.frame(variable=dataSetInf6m[,names(dataSetInf6m) %in% variable],
                                           milkConsumptionType=dataSetInf6m$milkConsumptionType))

    modelForm=as.formula(paste("variable~milkConsumptionType"  ));

    modelInfo <- lm(modelForm, data = thisDataInstance)
    statModelAnova <- aov(modelInfo)
    tukey <- suppressWarnings( TukeyHSD(statModelAnova,"milkConsumptionType"))

    meanBF <- format(mean(thisDataInstance[thisDataInstance$milkConsumptionType %in% "BF",names(thisDataInstance) %in% "variable"]),digits=2)
    meanBF_list[[length(meanBF_list)+1]] <- meanBF;

    meanTF <- format(mean(thisDataInstance[thisDataInstance$milkConsumptionType %in% "TF",names(thisDataInstance) %in% "variable"]),digits=2)
    meanTF_list[[length(meanTF_list)+1]] <- meanTF;

    meanASF <- format(mean(thisDataInstance[thisDataInstance$milkConsumptionType %in% "ASF",names(thisDataInstance) %in% "variable"]),digits=2)
    meanASF_list[[length(meanASF_list)+1]] <- meanASF;

    sdBF <- format(sd(thisDataInstance[thisDataInstance$milkConsumptionType %in% "BF",names(thisDataInstance) %in% "variable"]),digits=2)
    sdBF_list[[length(sdBF_list)+1]] <- sdBF;

    sdTF <- format(sd(thisDataInstance[thisDataInstance$milkConsumptionType %in% "TF",names(thisDataInstance) %in% "variable"]),digits=2)
    sdTF_list[[length(sdTF_list)+1]] <- sdTF;

    sdASF <- format(sd(thisDataInstance[thisDataInstance$milkConsumptionType %in% "ASF",names(thisDataInstance) %in% "variable"]),digits=2)
    sdASF_list[[length(sdASF_list)+1]] <- sdASF;

    if(categorical %in% TRUE)
    {
      pval=suppressWarnings(chisq.test(as.table(rbind(thisDataInstance$milkConsumptionType,thisDataInstance[,names(thisDataInstance) %in% "variable"])))$p.value)
    }
    if(categorical %in% FALSE)
    {
      modelInfo <- lm(modelForm, data = thisDataInstance)
      pval <- anova(modelInfo)$"Pr(>F)"[1];
    }

    pValSugar <- format.pval(pval,digits=3);
    pValSugarList[[length(pValSugarList)+1]] <- pValSugar;

    # fdr_sugar_value <- paste(format.pval(fdr_sugar[fdr_sugar$threshold %in% (floor(-log10(anova(modelInfo)$"Pr(>F)"[1])*10^1)/10^1),][2],digits=2)," (",
    #                          format.pval(fdr_sugar[fdr_sugar$threshold %in% (floor(-log10(anova(modelInfo)$"Pr(>F)"[1])*10^1)/10^1),][3],digits=2),",",
    #                          format.pval(fdr_sugar[fdr_sugar$threshold %in% (floor(-log10(anova(modelInfo)$"Pr(>F)"[1])*10^1)/10^1),][4],digits=2),")",sep="")
    #
    # sugar_fdrci_list[[length(sugar_fdrci_list)+1]] <- fdr_sugar_value;
  }
  bfList <- paste(meanBF_list," (",sdBF_list,")",sep="")
  tfList <- paste(meanTF_list," (",sdTF_list,")",sep="")
  asfList <- paste(meanASF_list," (",sdASF_list,")",sep="")

  # allOutcomes <- cbind(allOutcomes,
  #                       bfList,tfList,asfList,pValSugarList,sugar_fdrci_list);
  allOutcomes <- cbind(allOutcomes,
                        bfList,tfList,asfList,pValSugarList);

  allOutcomes <- data.frame(allOutcomes)
  return(allOutcomes)
}



compute_bayley_milkType_association <- function()
{
  dataSet <-thisDataInstance_all
  dataSetInf6m  <- dataSet[grepl("Inf.6m",dataSet$X.SampleID),]

  bayleysVariables <- c("bsid_cog_ss","bsid_lang_rc_ss",
                        "bsid_lang_ec_ss","bsid_lang_ss","bsid_mot_fm_ss","bsid_mot_gm_ss",
                        "bsid_mot_ss","bsid_se_ss","bsid_ab_com_ss","bsid_ab_cu_ss",
                        "bsid_ab_fa_ss","bsid_ab_hl_ss","bsid_ab_hs_ss","bsid_ab_ls_ss","bsid_ab_sc_ss","bsid_ab_sd_ss","bsid_ab_soc_ss",
                        "bsid_ab_mo_ss","bsid_ab_ss")
  bayleysNames <- c("Cognitive Scaled Score","Receptive Communication (RC) Total Scaled Score",
                    "Expressive Communication (EC) Total Scaled Score","Language Sum Scaled Score","Fine Motor (FM) Scaled Score","Gross Motor (GM) Scaled Score",
                    "Motor Sum Scaled Score","Social-Emotional Scaled Score","Communication (Com) Scaled Score","Community Use (CU) Scaled Score",
                    "Functional Pre-Academics (FA) Scaled Score","Home Living (HL) Scaled Score","Health and Safety (HS) Scaled Score","Leisure (LS) Scaled Score",
                    "Self-Care (SC) Scaled Score","Self-Direction (SD) Scaled Score","Social (Soc) Scaled Score",
                    "Motor (MO) Scaled Score","Adaptive Behavior Sum Scaled Score")

  fdr_asfbf <- bayleysFDRCI_asfbf
  fdr_tfbf <- bayleysFDRCI_tfbf

  exposureNamesList <- character(0);allOutcomes <- character(0);

  meanBF_list <- numeric(0);meanTF_list <- numeric(0);meanASF_list <- numeric(0);
  sdBF_list <- numeric(0);sdTF_list <- numeric(0);sdASF_list <- numeric(0);
  diffTFBF_list <- numeric(0);diffASFBF_list <- numeric(0);
  lwrTFBF_list <- numeric(0);lwrASFBF_list <- numeric(0);
  uprTFBF_list <- numeric(0);uprASFBF_list <- numeric(0);
  pValTFBF_list <- numeric(0);pValASFBF_list <- numeric(0);
  tfbf_fdrci_list <- character(0);asfbf_fdrci_list <- character(0);

  i <- 1
  for( variable in bayleysVariables)
  {
    allOutcomes[[length(allOutcomes)+1]] <- variable;
    thisDataInstance <- na.omit(data.frame(dyad_id=dataSetInf6m$dyad_id,
                                           thisMeasure=as.numeric(as.character(dataSetInf6m[,names(dataSetInf6m) %in% variable])),
                                           milkConsumptionType=dataSetInf6m$milkConsumptionType,
                                           prepreg_bmi_kgm2=dataSetInf6m$prepreg_bmi_kgm2,
                                           sex=factor(dataSetInf6m$sex),
                                           fruit_incj_Inf=dataSetInf6m$fruit_incj_Inf_24,
                                           inf_weight_kg=dataSetInf6m$inf_weight_kg))

    modelForm=as.formula(paste("thisMeasure","~milkConsumptionType+prepreg_bmi_kgm2 +sex+inf_weight_kg+fruit_incj_Inf"  ));


    modelInfo <- lm(modelForm, data = thisDataInstance)
    statModelAnova <- aov(modelInfo)
    tukey <- suppressWarnings( TukeyHSD(statModelAnova,"milkConsumptionType"))

    meanBF <- format(mean(thisDataInstance[thisDataInstance$milkConsumptionType %in% "BF",names(thisDataInstance) %in% "thisMeasure"]),digits=2)
    meanBF_list[[length(meanBF_list)+1]] <- meanBF;

    meanTF <- format(mean(thisDataInstance[thisDataInstance$milkConsumptionType %in% "TF",names(thisDataInstance) %in% "thisMeasure"]),digits=2)
    meanTF_list[[length(meanTF_list)+1]] <- meanTF;

    meanASF <- format(mean(thisDataInstance[thisDataInstance$milkConsumptionType %in% "ASF",names(thisDataInstance) %in% "thisMeasure"]),digits=2)
    meanASF_list[[length(meanASF_list)+1]] <- meanASF;

    sdBF <- format(sd(thisDataInstance[thisDataInstance$milkConsumptionType %in% "BF",names(thisDataInstance) %in% "thisMeasure"]),digits=2)
    sdBF_list[[length(sdBF_list)+1]] <- sdBF;

    sdTF <- format(sd(thisDataInstance[thisDataInstance$milkConsumptionType %in% "TF",names(thisDataInstance) %in% "thisMeasure"]),digits=2)
    sdTF_list[[length(sdTF_list)+1]] <- sdTF;

    sdASF <- format(sd(thisDataInstance[thisDataInstance$milkConsumptionType %in% "ASF",names(thisDataInstance) %in% "thisMeasure"]),digits=2)
    sdASF_list[[length(sdASF_list)+1]] <- sdASF;


    pValTFBF <- format.pval(tukey$milkConsumptionType[1,4],digits=2);
    pValTFBF_list[[length(pValTFBF_list)+1]] <- pValTFBF;
    # fdr_tfbf_value <- paste(format.pval(fdr_tfbf[fdr_tfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[1,4])*10^3)/10^3),][2],digits=2)," (",
    #                         format.pval(fdr_tfbf[fdr_tfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[1,4])*10^3)/10^3),][3],digits=2),",",
    #                         format.pval(fdr_tfbf[fdr_tfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[1,4])*10^3)/10^3),][4],digits=2),")",sep="")
    # tfbf_fdrci_list[[length(tfbf_fdrci_list)+1]] <- fdr_tfbf_value;

    pValASFBF <- format.pval(tukey$milkConsumptionType[2,4],digits=2);
    pValASFBF_list[[length(pValASFBF_list)+1]] <- pValASFBF;
    # fdr_asfbf_value <- paste(format.pval(fdr_asfbf[fdr_asfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[2,4])*10^3)/10^3),][2],digits=2)," (",
    #                          format.pval(fdr_asfbf[fdr_asfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[2,4])*10^3)/10^3),][3],digits=2),",",
    #                          format.pval(fdr_asfbf[fdr_asfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[2,4])*10^3)/10^3),][4],digits=2),")",sep="")
    # asfbf_fdrci_list[[length(asfbf_fdrci_list)+1]] <- fdr_asfbf_value;

    diffTFBF <- format(tukey$milkConsumptionType[1,1],digits=2);
    diffTFBF_list[[length(diffTFBF_list)+1]] <- diffTFBF;

    diffASFBF <- format(tukey$milkConsumptionType[2,1],digits=2);
    diffASFBF_list[[length(diffASFBF_list)+1]] <- diffASFBF;

    lwrTFBF <- format(tukey$milkConsumptionType[1,2],digits=2);
    lwrTFBF_list[[length(lwrTFBF_list)+1]] <- lwrTFBF;

    lwrASFBF <- format(tukey$milkConsumptionType[2,2],digits=2);
    lwrASFBF_list[[length(lwrASFBF_list)+1]] <- lwrASFBF;

    uprTFBF <- format(tukey$milkConsumptionType[1,3],digits=2);
    uprTFBF_list[[length(uprTFBF_list)+1]] <- uprTFBF;

    uprASFBF <- format(tukey$milkConsumptionType[2,3],digits=2);
    uprASFBF_list[[length(uprASFBF_list)+1]] <- uprASFBF;
    i <- i+1;
  }
  bfList <- paste(meanBF_list," (",sdBF_list,")",sep="")
  tfList <- paste(meanTF_list," (",sdTF_list,")",sep="")
  asfList <- paste(meanASF_list," (",sdASF_list,")",sep="")

  difftfbfList <- paste(diffTFBF_list," (",lwrTFBF_list,", ",uprTFBF_list,")",sep="")
  diffasfbfList <- paste(diffASFBF_list," (",lwrASFBF_list,", ",uprASFBF_list,")",sep="")

  # allOutcomes <- cbind(allOutcomes,
  #                       bfList,tfList,asfList,
  #                       difftfbfList,pValTFBF_list,tfbf_fdrci_list,
  #                       diffasfbfList,pValASFBF_list,asfbf_fdrci_list);

  allOutcomes <- cbind(allOutcomes,
                        bfList,tfList,asfList,
                        difftfbfList,pValTFBF_list,
                        diffasfbfList,pValASFBF_list);

  allOutcomes <- data.frame(allOutcomes)
  return(allOutcomes)
}


compute_somaticGrowth_milkType_association <- function( )
{
  dataSet <-thisDataInstance_all
  dataSetInf6m  <- dataSet[grepl("Inf.6m",dataSet$X.SampleID),]

  somaticVariables <- c("inf_weight_kg","zbmi","zwfl","zwei","zlen","skinf_tricep_mm","skinf_subscap_mm",
                        "skinf_supra_mm","skinf_midthigh_mm","circ_umb_cm")
  somaticNames <- c("Infant Weight","BMI Z-score","Weight-for-Length Z-score","Weight Z-score","Length Z-score","Tricep Skinfold (mm)","Subscapular Skinfold (mm)",
                    "Suprailiac Skinfold (mm)","Midthigh Skinfold (mm)","Abdominal Circumference (cm)")

  fdr_asfbf <- anthroFDRCI_asfbf
  fdr_tfbf <- anthroFDRCI_tfbf


  exposureNamesList <- character(0);allOutcomes <- character(0);

  meanBF_list <- numeric(0);meanTF_list <- numeric(0);meanASF_list <- numeric(0);
  sdBF_list <- numeric(0);sdTF_list <- numeric(0);sdASF_list <- numeric(0);
  diffTFBF_list <- numeric(0);diffASFBF_list <- numeric(0);
  lwrTFBF_list <- numeric(0);lwrASFBF_list <- numeric(0);
  uprTFBF_list <- numeric(0);uprASFBF_list <- numeric(0);
  pValTFBF_list <- numeric(0);pValASFBF_list <- numeric(0);
  tfbf_fdrci_list <- character(0);asfbf_fdrci_list <- character(0);

  i <- 1
  for( variable in somaticVariables)
  {
    allOutcomes[[length(allOutcomes)+1]] <- variable;
    thisDataInstance <- na.omit(data.frame(dyad_id=dataSetInf6m$dyad_id,
                                           thisMeasure=as.numeric(as.character(dataSetInf6m[,names(dataSetInf6m) %in% variable])),
                                           milkConsumptionType=dataSetInf6m$milkConsumptionType,
                                           prepreg_bmi_kgm2=dataSetInf6m$prepreg_bmi_kgm2,
                                           sex=factor(dataSetInf6m$sex),
                                           fruit_incj_Inf=dataSetInf6m$fruit_incj_Inf_24))


    modelForm=as.formula(paste("thisMeasure","~milkConsumptionType+prepreg_bmi_kgm2 +sex+fruit_incj_Inf"  ));


    modelInfo <- lm(modelForm, data = thisDataInstance)
    statModelAnova <- aov(modelInfo)
    tukey <- suppressWarnings( TukeyHSD(statModelAnova,"milkConsumptionType"))

    meanBF <- format(mean(thisDataInstance[thisDataInstance$milkConsumptionType %in% "BF",names(thisDataInstance) %in% "thisMeasure"]),digits=2)
    meanBF_list[[length(meanBF_list)+1]] <- meanBF;

    meanTF <- format(mean(thisDataInstance[thisDataInstance$milkConsumptionType %in% "TF",names(thisDataInstance) %in% "thisMeasure"]),digits=2)
    meanTF_list[[length(meanTF_list)+1]] <- meanTF;

    meanASF <- format(mean(thisDataInstance[thisDataInstance$milkConsumptionType %in% "ASF",names(thisDataInstance) %in% "thisMeasure"]),digits=2)
    meanASF_list[[length(meanASF_list)+1]] <- meanASF;

    sdBF <- format(sd(thisDataInstance[thisDataInstance$milkConsumptionType %in% "BF",names(thisDataInstance) %in% "thisMeasure"]),digits=2)
    sdBF_list[[length(sdBF_list)+1]] <- sdBF;

    sdTF <- format(sd(thisDataInstance[thisDataInstance$milkConsumptionType %in% "TF",names(thisDataInstance) %in% "thisMeasure"]),digits=2)
    sdTF_list[[length(sdTF_list)+1]] <- sdTF;

    sdASF <- format(sd(thisDataInstance[thisDataInstance$milkConsumptionType %in% "ASF",names(thisDataInstance) %in% "thisMeasure"]),digits=2)
    sdASF_list[[length(sdASF_list)+1]] <- sdASF;


    pValTFBF <- format.pval(tukey$milkConsumptionType[1,4],digits=2);
    pValTFBF_list[[length(pValTFBF_list)+1]] <- pValTFBF;
    # fdr_tfbf_value <- paste(format.pval(fdr_tfbf[fdr_tfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[1,4])*10^3)/10^3),][2],digits=2)," (",
    #                         format.pval(fdr_tfbf[fdr_tfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[1,4])*10^3)/10^3),][3],digits=2),",",
    #                         format.pval(fdr_tfbf[fdr_tfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[1,4])*10^3)/10^3),][4],digits=2),")",sep="")
    # tfbf_fdrci_list[[length(tfbf_fdrci_list)+1]] <- fdr_tfbf_value;

    pValASFBF <- format.pval(tukey$milkConsumptionType[2,4],digits=2);
    pValASFBF_list[[length(pValASFBF_list)+1]] <- pValASFBF;
    # fdr_asfbf_value <- paste(format.pval(fdr_asfbf[fdr_asfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[2,4])*10^3)/10^3),][2],digits=2)," (",
    #                          format.pval(fdr_asfbf[fdr_asfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[2,4])*10^3)/10^3),][3],digits=2),",",
    #                          format.pval(fdr_asfbf[fdr_asfbf$threshold %in% (floor(-log10(tukey$milkConsumptionType[2,4])*10^3)/10^3),][4],digits=2),")",sep="")
    # asfbf_fdrci_list[[length(asfbf_fdrci_list)+1]] <- fdr_asfbf_value;

    diffTFBF <- format(tukey$milkConsumptionType[1,1],digits=2);
    diffTFBF_list[[length(diffTFBF_list)+1]] <- diffTFBF;

    diffASFBF <- format(tukey$milkConsumptionType[2,1],digits=2);
    diffASFBF_list[[length(diffASFBF_list)+1]] <- diffASFBF;


    lwrTFBF <- format(tukey$milkConsumptionType[1,2],digits=2);
    lwrTFBF_list[[length(lwrTFBF_list)+1]] <- lwrTFBF;

    lwrASFBF <- format(tukey$milkConsumptionType[2,2],digits=2);
    lwrASFBF_list[[length(lwrASFBF_list)+1]] <- lwrASFBF;


    uprTFBF <- format(tukey$milkConsumptionType[1,3],digits=2);
    uprTFBF_list[[length(uprTFBF_list)+1]] <- uprTFBF;

    uprASFBF <- format(tukey$milkConsumptionType[2,3],digits=2);
    uprASFBF_list[[length(uprASFBF_list)+1]] <- uprASFBF;

    i <- i+1;
  }
  bfList <- paste(meanBF_list," (",sdBF_list,")",sep="")
  tfList <- paste(meanTF_list," (",sdTF_list,")",sep="")
  asfList <- paste(meanASF_list," (",sdASF_list,")",sep="")

  difftfbfList <- paste(diffTFBF_list," (",lwrTFBF_list,", ",uprTFBF_list,")",sep="")
  diffasfbfList <- paste(diffASFBF_list," (",lwrASFBF_list,", ",uprASFBF_list,")",sep="")

  # allOutcomes <- cbind(allOutcomes,
  #                       bfList,tfList,asfList,
  #                       difftfbfList,pValTFBF_list,tfbf_fdrci_list,
  #                       diffasfbfList,pValASFBF_list,asfbf_fdrci_list);

  allOutcomes <- cbind(allOutcomes,
                        bfList,tfList,asfList,
                        difftfbfList,pValTFBF_list,
                        diffasfbfList,pValASFBF_list);

  allOutcomes <- data.frame(allOutcomes)
  return(allOutcomes)
}
