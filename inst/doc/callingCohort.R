## ----createCountTables, eval=F------------------------------------------------
#  countTables <- baseCountsFromBamList( list.files("bam", full.names = T) )

## ----downloadCountTables------------------------------------------------------
countTables <- readRDS(url("http://steinmetzlab.embl.de/mutaseq/nuc.count.per.position.RDS"))
countTables <- countTables[grep('Q59',names(countTables),invert=T)]
print(head(countTables[[1]]))

## ----parsePatient-------------------------------------------------------------
patient <- gsub("_.+","",names(countTables))
print(table(patient))

## ----pressure-----------------------------------------------------------------
library(mitoClone)
result <- mutationCallsFromCohort(countTables, patient)
print(colnames(result$P342@M))
print(colnames(result$HRK@M))

## ----get_association, fig.width=6,fig.height=4--------------------------------
#retrieve the binary mutation status for nuclear and mitochondrial genomic mutation from read count tables (corresponds to result$P342 , but additionally contains nuclear mutation calls)
mutations <- mutationCallsFromMatrix(as.matrix(M_P1), as.matrix(N_P1))


getasso <- function(m1,m2) {
      testm <- matrix( data = c(
        sum(m1 =="1" &m2 == "1"),
        sum(m1 =="1" &m2 == "0"),
        sum(m1 == "0" & m2 == "1"),
        sum(m1 == "0" & m2 == "0")
      ), ncol = 2  )
      fisher.test(testm)
}


#run association tests between binary mitochondrial mutation status and binary SRSF2 mutation status
mitochondrial <- colnames(mutations@ternary)[grep("^X\\d+",colnames(mutations@ternary))]
tests <- apply(mutations@ternary[,mitochondrial], 2, getasso, m2 = mutations@ternary[,"SRSF2"])

p.val <- sapply(tests, "[[", "p.value")
est <- sapply(tests, "[[", "estimate")
plf <- data.frame(name = mitochondrial,
                  pval = -log10(p.val),
                  or = est)

library(ggplot2)
qplot(x = name , y = pval, color = paste0(or <1, pval > 1), data=na.omit(plf)) + coord_flip() + theme_bw()+ theme(panel.grid = element_blank(), axis.title.x = element_text(color="black"), axis.text = element_text(color="black"), axis.text.y = element_text(size = 6)) + xlab("") + ylab("Association w/ SRSF2 mutation (-log10 p)")+ scale_color_manual(values = c("TRUETRUE" = "blue","FALSETRUE" = "red", "FALSEFALSE" = "grey", "TRUEFALSE" = "grey"), guide=F)



