plot_phylum_relativeAbundances_barplot <- function( )
{
  dataSet <-thisDataInstance_all
  dataSetInf6m  <- dataSet[grepl("Inf.6m",dataSet$X.SampleID),]
  dataSetInf6m <- dataSetInf6m[,names(dataSetInf6m) %in% c("X.SampleID","dyad_id","milkConsumptionType")]

  taxaTable <- phylumTaxa
  taxaTable <- as.data.frame(taxaTable)
  taxaTable$taxa <- gsub("\\| ",".",taxaTable$taxa)
  row.names(taxaTable) <- taxaTable$taxa
  taxaTable_6m<- taxaTable[,grepl("Inf.6m",colnames(taxaTable))]
  myT <- t(taxaTable_6m)
  myT <- myT[order(row.names(myT)),]

  otherMicrobes <- as.data.frame(rowSums(myT[,!colnames(myT) %in% c("k__Bacteria.p__Actinobacteria","k__Bacteria.p__Bacteroidetes","k__Bacteria.p__Firmicutes","k__Bacteria.p__Proteobacteria","k__Bacteria.p__Verrucomicrobia")]))
  names(otherMicrobes)[1] <- "Others"
  myT <- myT[,colnames(myT) %in% c("k__Bacteria.p__Actinobacteria","k__Bacteria.p__Bacteroidetes","k__Bacteria.p__Firmicutes","k__Bacteria.p__Proteobacteria","k__Bacteria.p__Verrucomicrobia")]
  myT <- merge(myT,otherMicrobes,by="row.names")
  row.names(myT) <- myT[,1]
  myT <- myT[,-1]
  myT <- myT[c("k__Bacteria.p__Actinobacteria","k__Bacteria.p__Bacteroidetes","k__Bacteria.p__Firmicutes","k__Bacteria.p__Proteobacteria","k__Bacteria.p__Verrucomicrobia","Others")]
  colnames(myT)[1:5] <- substr(colnames(myT)[1:5],13,nchar(colnames(myT)[1:5]))
  myT2 <- myT/rowSums(myT)
  myT3 <- myT2[,(colSums(myT2==0)/nrow(myT2)) <= 0.99]

  taxaMeta <- merge(dataSetInf6m,myT3,by.x="X.SampleID",by.y="row.names")

  #taxaMeta2 <-  taxaMeta[order(taxaMeta$frucLevel,taxaMeta$CAGE),];
  taxaMeta2 <-  taxaMeta[order(taxaMeta$milkConsumptionType,taxaMeta$`p__Actinobacteria`,taxaMeta$`p__Bacteroidetes`,taxaMeta$`p__Firmicutes`),];
  taxaMeta3 <- taxaMeta2[,!colnames(taxaMeta2) %in% c("X.SampleID", "dyad_id")]

  phylum_colors <- c("#8ca63a",
                     "#9964cb",
                     "#5fa271",
                     "#c95a80",
                     "#798bc6",
                     "#c76f3e")

  summaries <- taxaMeta3 %>%
    group_by(milkConsumptionType) %>%
    summarize_if(is.numeric, mean)
  summariesnames <- summaries$milkConsumptionType
  summaries <- summaries[,-1]
  summaries <- as.matrix(summaries)
  row.names(summaries) <- summariesnames
  dat <- melt(summaries,id.vars="milkConsumptionType")
  names(dat)[names(dat) %in% "Var1"] <- "milkConsumptionType"

  p <- ggplot(dat, aes(x=milkConsumptionType, y=value, fill = Var2))
  print(p +geom_bar(position="stack", stat="identity",
                    colour="black",# Use black outlines,
                    size=.5) +      # Thinner lines
          scale_fill_manual(values=phylum_colors) +
          #guides(fill=FALSE)+
          ggtitle("Phylum-level")+
          ylab("Percentage") + xlab("Milk ConsumptionType") +
          theme_classic(base_size = 20)+
          theme(axis.line=element_line(size=1),
                axis.ticks.y=element_line(size=1),
                axis.ticks.x=element_blank(),
                axis.text=element_text(face="bold",size=16),
                text=element_text(face="bold",size=20),
                legend.position="left",
                legend.title=element_blank()
          )+
          theme(axis.line.x = element_line(color="black", size = 2),
                axis.line.y = element_line(color="black", size = 2)
          )
  )
}

plot_family_relativeAbundances_barplot <- function( )
{
  dataSet <-thisDataInstance_all
  dataSetInf6m  <- dataSet[grepl("Inf.6m",dataSet$X.SampleID),]
  dataSetInf6m <- dataSetInf6m[,names(dataSetInf6m) %in% c("X.SampleID","dyad_id","milkConsumptionType")]

  taxaTable <- familyTaxa
  taxaTable <- as.data.frame(taxaTable)
  taxaTable$taxa <- gsub("\\| ",".",taxaTable$taxa)
  row.names(taxaTable) <- taxaTable$taxa
  taxaTable_6m<- taxaTable[,grepl("Inf.6m",colnames(taxaTable))]
  myT <- t(taxaTable_6m)
  myT <- myT[order(row.names(myT)),]

  mainMicrobes <- names(colSums(myT)[colSums(myT)>30000])
  microbeName <- character(0)
  for(i in 1:length(mainMicrobes))
    microbeName[[length(microbeName)+1]] <- str_split(mainMicrobes,"__")[[i]][6]
  otherMicrobes <- as.data.frame(rowSums(myT[,!colnames(myT) %in% c(mainMicrobes)]))
  names(otherMicrobes)[1] <- "Others"
  myT <- myT[,colnames(myT) %in% c(mainMicrobes)]
  myT <- merge(myT,otherMicrobes,by="row.names")
  row.names(myT) <- myT[,1]
  myT <- myT[,-1]
  myT <- myT[c(mainMicrobes,"Others")]
  colnames(myT)[1:length(mainMicrobes)] <- microbeName

  myT2 <- myT/rowSums(myT)
  myT3 <- myT2[,(colSums(myT2==0)/nrow(myT2)) <= 0.99]

  taxaMeta <- merge(dataSetInf6m,myT3,by.x="X.SampleID",by.y="row.names")

  #taxaMeta2 <-  taxaMeta[order(taxaMeta$frucLevel,taxaMeta$CAGE),];
  taxaMeta2 <-  taxaMeta[order(taxaMeta$milkConsumptionType,taxaMeta$`Bacteroidaceae`,taxaMeta$`Enterobacteriaceae`,taxaMeta$`Bifidobacteriaceae`,taxaMeta$`Lachnospiraceae`),];
  taxaMeta3 <- taxaMeta2[,!colnames(taxaMeta2) %in% c("X.SampleID", "dyad_id")]

  family_colors <- c("#aa60da",
                     "#70b836",
                     "#da45af",
                     "#57bf78",
                     "#db5680",
                     "#39c2ca",
                     "#d95d48",
                     "#6780d9",
                     "#aaa13c",
                     "#c07dc2",
                     "#5f924c",
                     "#d0883f")

  summaries <- taxaMeta3 %>%
    group_by(milkConsumptionType) %>%
    summarize_if(is.numeric, mean)
  summariesnames <- summaries$milkConsumptionType
  summaries <- summaries[,-1]
  summaries <- as.matrix(summaries)
  row.names(summaries) <- summariesnames
  dat <- melt(summaries,id.vars="milkConsumptionType")
  names(dat)[names(dat) %in% "Var1"] <- "milkConsumptionType"

  p <- ggplot(dat, aes(x=milkConsumptionType, y=value, fill = Var2))
  print(p +geom_bar(position="stack", stat="identity",
                    colour="black",# Use black outlines,
                    size=.5) +      # Thinner lines
          scale_fill_manual(values=family_colors) +
          ggtitle("Family-level")+
          ylab("Percentage") + xlab("Milk ConsumptionType") +
          theme_classic(base_size = 20)+
          theme(axis.line=element_line(size=1),
                axis.ticks.y=element_line(size=1),
                axis.ticks.x=element_blank(),
                axis.text=element_text(face="bold",size=16),
                text=element_text(face="bold",size=20),
                legend.position="left",
                legend.title=element_blank()
          )+
          theme(axis.line.x = element_line(color="black", size = 2),
                axis.line.y = element_line(color="black", size = 2)
          )
  )
}

plot_taxa_milkType_association <- function(taxaData = phylumTaxa,variable = "k__Bacteria.p__Actinobacteria",variableName = "Logged Abundance of Actinobacteria")
{
  dataSet <-thisDataInstance_all
  dataSetInf1m  <- dataSet[grepl("Inf.Bl",dataSet$X.SampleID),]
  dataSetInf6m  <- dataSet[grepl("Inf.6m",dataSet$X.SampleID),]

  taxaTable <- taxaData
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

  modelForm=as.formula(paste("thisIndex","~milkConsumptionType+mother_age+prepreg_bmi_kgm2 +mode_of_delivery+mom_BMI+sex+ baby_age+inf_weight_kg+fruit_incj_Inf+thisIndex_1m"  ));

  thisDataInstance_6m <- na.omit(data.frame(dyad_id=taxaMeta_6m$dyad_id,
                                            thisIndex=as.numeric(as.character(taxaMeta_6m[,names(taxaMeta_6m) %in% variable])),
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
                                            thisIndex_1m=taxaMeta_1m[,colnames(taxaMeta_1m) %in% variable]))
  thisDataInstance_1m6m <- merge(thisDataInstance_6m,thisDataInstance_1m,by="dyad_id")

  modelInfo <- lm(modelForm, data = thisDataInstance_1m6m)
  statModelAnova <- aov(modelInfo)
  tukey <- suppressWarnings( TukeyHSD(statModelAnova,"milkConsumptionType"))

  p <- ggplot(thisDataInstance_1m6m,aes(x=factor(milkConsumptionType),y=thisIndex))
  print(p + geom_boxplot(outlier.size=-1) +
          geom_jitter(colour="grey46",size=5, position=position_jitter(0.2),aes(shape=milkConsumptionType)) +
          guides(shape=FALSE)+
          ggtitle("")+
          geom_signif(comparisons=list(c("TF","BF")),
                      annotations=ifelse(tukey$milkConsumptionType[1,4]<0.001,"****",
                                         ifelse(tukey$milkConsumptionType[1,4]<0.01,"***",
                                                ifelse(tukey$milkConsumptionType[1,4]<0.05,"**",
                                                       ifelse(tukey$milkConsumptionType[1,4]<0.1,"*","NS")))),
                      y_position = max(thisDataInstance_1m6m$thisIndex)+0.1, tip_length = 0.1, vjust=0.3) +
          geom_signif(comparisons=list(c("ASF","BF")),
                      annotations=ifelse(tukey$milkConsumptionType[2,4]<0.001,"****",
                                         ifelse(tukey$milkConsumptionType[2,4]<0.01,"***",
                                                ifelse(tukey$milkConsumptionType[2,4]<0.05,"**",
                                                       ifelse(tukey$milkConsumptionType[2,4]<0.1,"*","NS")))),
                      y_position = max(thisDataInstance_1m6m$thisIndex)+0.6, tip_length = 0.1, vjust=0.3) +

          xlab("Milk Consumption Type") +
          ylab(variableName)+
          theme_classic(base_size = 20)+
          theme(axis.line=element_line(size=1),
                axis.ticks.y=element_line(size=1),
                axis.ticks.x=element_blank(),
                axis.text=element_text(face="bold",size=20),
                text=element_text(face="bold",size=16),
                legend.position="bottom",
                legend.title=element_blank(),
                plot.title = element_text(hjust = 0.5)
          )+
          theme(axis.line.x = element_line(color="black", size = 2),
                axis.line.y = element_line(color="black", size = 2)
          )
  )
}

plot_alphaDiversity_milkType_association <- function(variable = "shannonDiversity",variableName = "Shannon Diversity")
{
  dataSet <-thisDataInstance_all
  dataSetInf1m  <- dataSet[grepl("Inf.Bl",dataSet$X.SampleID),]
  dataSetInf6m  <- dataSet[grepl("Inf.6m",dataSet$X.SampleID),]

  taxaTable_6m<- alphaDiversityTable[grepl("Inf.6m",rownames(alphaDiversityTable)),]
  taxaTable_1m <- alphaDiversityTable[grepl("Inf.Bl",rownames(alphaDiversityTable)),]

  taxaMeta_6m <- merge(dataSetInf6m,taxaTable_6m,by.x="X.SampleID",by.y="row.names")
  taxaMeta_1m <- merge(dataSetInf1m,taxaTable_1m,by.x="X.SampleID",by.y="row.names")

  modelForm=as.formula(paste("thisIndex","~milkConsumptionType+mother_age+prepreg_bmi_kgm2 +mode_of_delivery+mom_BMI+sex+ baby_age+inf_weight_kg+fruit_incj_Inf+thisIndex_1m"  ));

  thisDataInstance_6m <- na.omit(data.frame(dyad_id=taxaMeta_6m$dyad_id,
                                            thisIndex=as.numeric(as.character(taxaMeta_6m[,names(taxaMeta_6m) %in% variable])),
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
                                            thisIndex_1m=taxaMeta_1m[,colnames(taxaMeta_1m) %in% variable]))
  thisDataInstance_1m6m <- merge(thisDataInstance_6m,thisDataInstance_1m,by="dyad_id")

  modelInfo <- lm(modelForm, data = thisDataInstance_1m6m)
  statModelAnova <- aov(modelInfo)
  tukey <- suppressWarnings( TukeyHSD(statModelAnova,"milkConsumptionType"))

  p <- ggplot(thisDataInstance_1m6m,aes(x=factor(milkConsumptionType),y=thisIndex))
  print(p + geom_boxplot(outlier.size=-1) +
          geom_jitter(colour="grey46",size=5, position=position_jitter(0.2),aes(shape=milkConsumptionType)) +
          guides(shape=FALSE)+
          ggtitle("")+
          geom_signif(comparisons=list(c("TF","BF")),
                      annotations=ifelse(tukey$milkConsumptionType[1,4]<0.001,"****",
                                         ifelse(tukey$milkConsumptionType[1,4]<0.01,"***",
                                                ifelse(tukey$milkConsumptionType[1,4]<0.05,"**",
                                                       ifelse(tukey$milkConsumptionType[1,4]<0.1,"*","NS")))),
                      y_position = max(thisDataInstance_1m6m$thisIndex)+0.1, tip_length = 0.1, vjust=0.3) +
          geom_signif(comparisons=list(c("ASF","BF")),
                      annotations=ifelse(tukey$milkConsumptionType[2,4]<0.001,"****",
                                         ifelse(tukey$milkConsumptionType[2,4]<0.01,"***",
                                                ifelse(tukey$milkConsumptionType[2,4]<0.05,"**",
                                                       ifelse(tukey$milkConsumptionType[2,4]<0.1,"*","NS")))),
                      y_position = max(thisDataInstance_1m6m$thisIndex)+0.6, tip_length = 0.1, vjust=0.3) +

          xlab("Milk Consumption Type") +
          ylab(variableName)+
          theme_classic(base_size = 20)+
          theme(axis.line=element_line(size=1),
                axis.ticks.y=element_line(size=1),
                axis.ticks.x=element_blank(),
                axis.text=element_text(face="bold",size=20),
                text=element_text(face="bold",size=16),
                legend.position="bottom",
                legend.title=element_blank(),
                plot.title = element_text(hjust = 0.5)
          )+
          theme(axis.line.x = element_line(color="black", size = 2),
                axis.line.y = element_line(color="black", size = 2)
          )
  )
}

plot_betaDiversity_milkType_association <- function(axis = 1)
{
  dataSet <-thisDataInstance_all
  dataSetInf1m  <- dataSet[grepl("Inf.Bl",dataSet$X.SampleID),]
  dataSetInf6m  <- dataSet[grepl("Inf.6m",dataSet$X.SampleID),]

  loadings6m <- eigen6m


  taxaMeta_1m <- merge(dataSetInf1m,ordinationData1m,by.x="X.SampleID",by.y="row.names")
  taxaMeta_6m <- merge(dataSetInf6m,ordinationData6m,by.x="X.SampleID",by.y="row.names")

  modelForm=as.formula(paste("thisIndex","~milkConsumptionType+mother_age+prepreg_bmi_kgm2 +mode_of_delivery+mom_BMI+sex+ baby_age+inf_weight_kg+fruit_incj_Inf+thisIndex_1m"  ));

  thisDataInstance_6m <- na.omit(data.frame(dyad_id=taxaMeta_6m$dyad_id,
                                            thisIndex=taxaMeta_6m[,names(taxaMeta_6m) %in% paste("Axis",axis,sep = "")],
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
                                            thisIndex_1m=taxaMeta_1m[,colnames(taxaMeta_1m) %in% paste("Axis",axis,sep = "")]))
  thisDataInstance_1m6m <- merge(thisDataInstance_6m,thisDataInstance_1m,by="dyad_id")

  modelInfo <- lm(modelForm, data = thisDataInstance_1m6m)
  statModelAnova <- aov(modelInfo)
  tukey <- suppressWarnings( TukeyHSD(statModelAnova,"milkConsumptionType"))

  p <- ggplot(thisDataInstance_1m6m,aes(x=factor(milkConsumptionType),y=thisIndex))
  print(p + geom_boxplot(outlier.size=-1) +
          geom_jitter(colour="grey46",size=5, position=position_jitter(0.2),aes(shape=milkConsumptionType)) +
          guides(shape=FALSE)+
          ggtitle("")+
          geom_signif(comparisons=list(c("TF","BF")),
                      annotations=ifelse(tukey$milkConsumptionType[1,4]<0.001,"****",
                                         ifelse(tukey$milkConsumptionType[1,4]<0.01,"***",
                                                ifelse(tukey$milkConsumptionType[1,4]<0.05,"**",
                                                       ifelse(tukey$milkConsumptionType[1,4]<0.1,"*","NS")))),
                      y_position = max(thisDataInstance_1m6m$thisIndex)+0.1, tip_length = 0.1, vjust=0.3) +
          geom_signif(comparisons=list(c("ASF","BF")),
                      annotations=ifelse(tukey$milkConsumptionType[2,4]<0.001,"****",
                                         ifelse(tukey$milkConsumptionType[2,4]<0.01,"***",
                                                ifelse(tukey$milkConsumptionType[2,4]<0.05,"**",
                                                       ifelse(tukey$milkConsumptionType[2,4]<0.1,"*","NS")))),
                      y_position = max(thisDataInstance_1m6m$thisIndex)+0.6, tip_length = 0.1, vjust=0.3) +

          xlab("Milk Consumption Type") +
          ylab(paste("Axis",axis," ",(round(loadings6m$ProportionExplained[axis],2))*100,"%",sep=""))+
          theme_classic(base_size = 20)+
          theme(axis.line=element_line(size=1),
                axis.ticks.y=element_line(size=1),
                axis.ticks.x=element_blank(),
                axis.text=element_text(face="bold",size=20),
                text=element_text(face="bold",size=16),
                legend.position="bottom",
                legend.title=element_blank(),
                plot.title = element_text(hjust = 0.5)
          )+
          theme(axis.line.x = element_line(color="black", size = 2),
                axis.line.y = element_line(color="black", size = 2)
          )
  )
}

# plot_deicode_pca_plots <- function(dataFile="data/metaData.RData", betaDiversityFile_1m = "data/betaDiversity_1m.RData",betaDiversityFile_6m = "data/betaDiversity_6m.RData",loadings = "data/loadings_6m.RData",axis = 1)
# {
#   load(dataFile)
#   dataSet <-thisDataInstance_all
#   dataSetInf1m  <- dataSet[grepl("Inf.Bl",dataSet$X.SampleID),]
#   dataSetInf6m  <- dataSet[grepl("Inf.6m",dataSet$X.SampleID),]
#
#   load(betaDiversityFile_6m)
#   load(betaDiversityFile_1m)
#   load(loadings);
#   loadings6m <- eigen6m
#
#   taxaMeta_1m <- merge(dataSetInf1m,ordinationData1m,by.x="X.SampleID",by.y="row.names")
#   taxaMeta_6m <- merge(dataSetInf6m,ordinationData6m,by.x="X.SampleID",by.y="row.names")
#
#   modelForm1=as.formula(paste("Axis1","~milkConsumptionType+mother_age+prepreg_bmi_kgm2 +mode_of_delivery+mom_BMI+sex+ baby_age+inf_weight_kg+fruit_incj_Inf+Axis1_1m"  ));
#   modelForm2=as.formula(paste("Axis2","~milkConsumptionType+mother_age+prepreg_bmi_kgm2 +mode_of_delivery+mom_BMI+sex+ baby_age+inf_weight_kg+fruit_incj_Inf+Axis2_1m"  ));
#
#   thisDataInstance_6m <- na.omit(data.frame(dyad_id=taxaMeta_6m$dyad_id,
#                                             Axis1=taxaMeta_6m$Axis1,
#                                             Axis2=taxaMeta_6m$Axis2,
#                                             milkConsumptionType=taxaMeta_6m$milkConsumptionType,
#                                             mother_age=taxaMeta_6m$mother_age,
#                                             prepreg_bmi_kgm2=taxaMeta_6m$prepreg_bmi_kgm2,
#                                             mode_of_delivery=taxaMeta_6m$mode_of_delivery,
#                                             mom_BMI=taxaMeta_6m$mom_current_BMI,
#                                             sex=taxaMeta_6m$sex,
#                                             baby_age=taxaMeta_6m$baby_age,
#                                             inf_weight_kg=taxaMeta_6m$inf_weight_kg_1m6m,
#                                             fruit_incj_Inf=taxaMeta_6m$fruit_incj_Inf_1m6m))
#
#   thisDataInstance_1m <- na.omit(data.frame(dyad_id=taxaMeta_1m$dyad_id,
#                                             Axis1_1m=taxaMeta_1m$Axis1,
#                                             Axis2_1m=taxaMeta_1m$Axis2))
#   thisDataInstance_1m6m <- merge(thisDataInstance_6m,thisDataInstance_1m,by="dyad_id")
#
#   reducedModMDS1_6m <- lm(Axis1~mother_age+prepreg_bmi_kgm2 +mode_of_delivery+mom_BMI+sex+ baby_age+inf_weight_kg+Axis1_1m,data=thisDataInstance_1m6m);
#   fullModMDS1_6m <- lm(Axis1~factor(milkConsumptionType)+mother_age+prepreg_bmi_kgm2 +mode_of_delivery+mom_BMI+sex+ baby_age+inf_weight_kg+Axis1_1m,data=thisDataInstance_1m6m);
#   pVal1 <- as.numeric(format(anova(reducedModMDS1_1m,fullModMDS1_1m)$`Pr(>F)`[2],digits=3));
#
#   reducedModMDS2_6m <- lm(Axis2~mother_age+prepreg_bmi_kgm2 +mode_of_delivery+mom_BMI+sex+ baby_age+inf_weight_kg+Axis2_1m,data=thisDataInstance_1m6m);
#   fullModMDS2_6m <- lm(Axis2~factor(milkConsumptionType)+mother_age+prepreg_bmi_kgm2 +mode_of_delivery+mom_BMI+sex+ baby_age+inf_weight_kg+Axis2_1m,data=thisDataInstance_1m6m);
#   pVal2 <- as.numeric(format(anova(reducedModMDS2_1m,fullModMDS2_1m)$`Pr(>F)`[2],digits=3));
#
#   comp1<-as.character(paste("Axis1"," ",(round(loadings6m$ProportionExplained[1],2))*100,"%, p-val = ",format.pval(pVal1,2),sep=""));
#   comp2<-as.character(paste("Axis2"," ",(round(loadings6m$ProportionExplained[2],2))*100,"%, p-val = ",format.pval(pVal2,2),sep=""));
#
#   thisDataInstance <- na.omit(data.frame(Axis1=ordinationMeta1m$Axis1,Axis2=ordinationMeta1m$Axis2,milkConsumptionType=ordinationMeta1m$milkConsumptionType))
#   p <- ggplot(thisDataInstance_1m6m)
#   print(p + geom_point(aes(Axis1,Axis2,colour = milkConsumptionType,shape=milkConsumptionType),size = 10) +
#           #scale_color_manual(values=c("#999999", "#E69F00", "darkorchid1"))+
#           scale_color_manual(values=c("black", "grey46", "grey76"))+
#           #guides(colour=FALSE)+
#           #guides(shape=FALSE)+
#           xlab(comp1) + ylab(comp2) +
#           theme_classic(base_size = 30)+
#           theme(text=element_text(face="bold",size=20),
#                 axis.line=element_line(size=1),
#                 axis.ticks=element_line(size=1),
#                 axis.text=element_text(face="bold",size=16),
#                 axis.title=element_text(face="bold",size=30),
#                 plot.margin=unit(c(9.0, 5.5, 6.5, 5.5), "points"),
#                 legend.position = "bottom"
#           )+
#           theme(axis.line.x = element_line(color="black", size = 2),
#                 axis.line.y = element_line(color="black", size = 2),
#                 legend.title = element_blank(),
#                 legend.text=element_text(size=20)
#           )
#   )
# }


plot_keggModule_milkType_association <- function(variable = "M00221", variableName = "M00221: Putative simple sugar transport system")
{
  dataSet <- thisDataInstance_all
  dataSetInf1m  <- dataSet[grepl("Inf.Bl",dataSet$X.SampleID),]
  dataSetInf6m  <- dataSet[grepl("Inf.6m",dataSet$X.SampleID),]

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

  modelForm=as.formula(paste("thisModule","~milkConsumptionType+mother_age+prepreg_bmi_kgm2 +mode_of_delivery+mom_BMI+sex+ baby_age+inf_weight_kg+fruit_incj_Inf+thisModule_1m"  ));

  thisDataInstance_6m <- na.omit(data.frame(dyad_id=moduleMeta_6m$dyad_id,
                                            thisModule=as.numeric(as.character(moduleMeta_6m[,names(moduleMeta_6m) %in% variable])),
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
                                            thisModule_1m=moduleMeta_1m[,colnames(moduleMeta_1m) %in% variable]))
  thisDataInstance_1m6m <- merge(thisDataInstance_6m,thisDataInstance_1m,by="dyad_id")

  modelInfo <- lm(modelForm, data = thisDataInstance_1m6m)
  statModelAnova <- aov(modelInfo)
  tukey <- suppressWarnings( TukeyHSD(statModelAnova,"milkConsumptionType"))

  p <- ggplot(thisDataInstance_1m6m,aes(x=factor(milkConsumptionType),y=thisModule))
  print(p + geom_boxplot(outlier.size=-1) +
          geom_jitter(colour="grey46",size=5, position=position_jitter(0.2),aes(shape=milkConsumptionType)) +
          guides(shape=FALSE)+
          ggtitle("")+
          geom_signif(comparisons=list(c("TF","BF")),
                      annotations=ifelse(tukey$milkConsumptionType[1,4]<0.001,"****",
                                         ifelse(tukey$milkConsumptionType[1,4]<0.01,"***",
                                                ifelse(tukey$milkConsumptionType[1,4]<0.05,"**",
                                                       ifelse(tukey$milkConsumptionType[1,4]<0.1,"*","NS")))),
                      y_position = max(thisDataInstance_1m6m$thisModule)+0.1, tip_length = 0.1, vjust=0.3) +
          geom_signif(comparisons=list(c("ASF","BF")),
                      annotations=ifelse(tukey$milkConsumptionType[2,4]<0.001,"****",
                                         ifelse(tukey$milkConsumptionType[2,4]<0.01,"***",
                                                ifelse(tukey$milkConsumptionType[2,4]<0.05,"**",
                                                       ifelse(tukey$milkConsumptionType[2,4]<0.1,"*","NS")))),
                      y_position = max(thisDataInstance_1m6m$thisModule)+0.6, tip_length = 0.1, vjust=0.3) +

          xlab("Milk Consumption Type") +
          ylab(variableName)+
          theme_classic(base_size = 20)+
          theme(axis.line=element_line(size=1),
                axis.ticks.y=element_line(size=1),
                axis.ticks.x=element_blank(),
                axis.text=element_text(face="bold",size=20),
                text=element_text(face="bold",size=16),
                legend.position="bottom",
                legend.title=element_blank(),
                plot.title = element_text(hjust = 0.5)
          )+
          theme(axis.line.x = element_line(color="black", size = 2),
                axis.line.y = element_line(color="black", size = 2)
          )
  )
}



plot_bayley_milkType_association <- function(variable = "bsid_cog_ss",variableName = "Cognitive Scaled Score")
{
  dataSet <-thisDataInstance_all
  dataSetInf6m  <- dataSet[grepl("Inf.6m",dataSet$X.SampleID),]

  thisDataInstance <- na.omit(data.frame(dyad_id=dataSetInf6m$dyad_id,
                                         thisMeasure=as.numeric(as.character(dataSetInf6m[,names(dataSetInf6m) %in% variable])),
                                         milkConsumptionType=dataSetInf6m$milkConsumptionType,
                                         prepreg_bmi_kgm2=dataSetInf6m$prepreg_bmi_kgm2,
                                         sex=factor(dataSetInf6m$sex),
                                         fruit_incj_Inf=dataSetInf6m$fruit_incj_Inf_24,
                                         inf_weight_kg=dataSetInf6m$inf_weight_kg))

  modelForm=as.formula(paste("thisMeasure","~milkConsumptionType+prepreg_bmi_kgm2 +sex+inf_weight_kg+fruit_incj_Inf"));

  modelInfo <- lm(modelForm, data = thisDataInstance)
  statModelAnova <- aov(modelInfo)
  tukey <- suppressWarnings( TukeyHSD(statModelAnova,"milkConsumptionType"))

  p <- ggplot(thisDataInstance,aes(x=factor(milkConsumptionType),y=thisMeasure))
  print(p + geom_boxplot(outlier.size=-1) +
          geom_jitter(colour="grey46",size=5, position=position_jitter(0.2),aes(shape=milkConsumptionType)) +
          guides(shape=FALSE)+
          ggtitle("")+
          geom_signif(comparisons=list(c("TF","BF")),
                      annotations=ifelse(tukey$milkConsumptionType[1,4]<0.001,"****",
                                         ifelse(tukey$milkConsumptionType[1,4]<0.01,"***",
                                                ifelse(tukey$milkConsumptionType[1,4]<0.05,"**",
                                                       ifelse(tukey$milkConsumptionType[1,4]<0.1,"*","NS")))),
                      y_position = max(thisDataInstance$thisMeasure)+0.1, tip_length = 0.1, vjust=0.3) +
          geom_signif(comparisons=list(c("ASF","BF")),
                      annotations=ifelse(tukey$milkConsumptionType[2,4]<0.001,"****",
                                         ifelse(tukey$milkConsumptionType[2,4]<0.01,"***",
                                                ifelse(tukey$milkConsumptionType[2,4]<0.05,"**",
                                                       ifelse(tukey$milkConsumptionType[2,4]<0.1,"*","NS")))),
                      y_position = max(thisDataInstance$thisMeasure)+0.6, tip_length = 0.1, vjust=0.3) +

          xlab("6m Milk Consumption Type") +
          ylab(paste("24m",variableName))+
          theme_classic(base_size = 20)+
          theme(axis.line=element_line(size=1),
                axis.ticks.y=element_line(size=1),
                axis.ticks.x=element_blank(),
                axis.text=element_text(face="bold",size=20),
                text=element_text(face="bold",size=16),
                legend.position="bottom",
                legend.title=element_blank(),
                plot.title = element_text(hjust = 0.5)
          )+
          theme(axis.line.x = element_line(color="black", size = 2),
                axis.line.y = element_line(color="black", size = 2)
          )
  )
}

plot_somaticGrowth_milkType_association <- function(variable = "inf_weight_kg",variableName = "Infant Weight")
{
  dataSet <-thisDataInstance_all
  dataSetInf6m  <- dataSet[grepl("Inf.6m",dataSet$X.SampleID),]

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

  p <- ggplot(thisDataInstance,aes(x=factor(milkConsumptionType),y=thisMeasure))
  print(p + geom_boxplot(outlier.size=-1) +
          geom_jitter(colour="grey46",size=5, position=position_jitter(0.2),aes(shape=milkConsumptionType)) +
          guides(shape=FALSE)+
          ggtitle("")+
          geom_signif(comparisons=list(c("TF","BF")),
                      annotations=ifelse(tukey$milkConsumptionType[1,4]<0.001,"****",
                                         ifelse(tukey$milkConsumptionType[1,4]<0.01,"***",
                                                ifelse(tukey$milkConsumptionType[1,4]<0.05,"**",
                                                       ifelse(tukey$milkConsumptionType[1,4]<0.1,"*","NS")))),
                      y_position = max(thisDataInstance$thisMeasure)+0.1, tip_length = 0.1, vjust=0.3) +
          geom_signif(comparisons=list(c("ASF","BF")),
                      annotations=ifelse(tukey$milkConsumptionType[2,4]<0.001,"****",
                                         ifelse(tukey$milkConsumptionType[2,4]<0.01,"***",
                                                ifelse(tukey$milkConsumptionType[2,4]<0.05,"**",
                                                       ifelse(tukey$milkConsumptionType[2,4]<0.1,"*","NS")))),
                      y_position = max(thisDataInstance$thisMeasure)+0.6, tip_length = 0.1, vjust=0.3) +

          xlab("6m Milk Consumption Type") +
          ylab(paste(variableName," at 24 m",sep=""))+
          theme_classic(base_size = 20)+
          theme(axis.line=element_line(size=1),
                axis.ticks.y=element_line(size=1),
                axis.ticks.x=element_blank(),
                axis.text=element_text(face="bold",size=20),
                text=element_text(face="bold",size=16),
                legend.position="bottom",
                legend.title=element_blank(),
                plot.title = element_text(hjust = 0.5)
          )+
          theme(axis.line.x = element_line(color="black", size = 2),
                axis.line.y = element_line(color="black", size = 2)
          )
  )
}

