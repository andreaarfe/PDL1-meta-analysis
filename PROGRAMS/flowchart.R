
# Computes the number of trials excluded at each step in the flowchart diagram
d <- read.csv("./LITERATURE SEARCH/IncludedTrials.csv")

# Number of abstracts published before November 31, 2018
n0 <- nrow(d)
n0

# Criteria 1: manuscript unaccessible by the authors 
list_of_reasons_c1 <- c("CANNOT RECOVER PDF")
d$C1 <- as.numeric(d$Reason %in% list_of_reasons_c1)
n1 <- sum(d[,"C1"])


# Criteria 2: not randomized phase II or III of PD-1/L1 inhibitors in adults (≥18 years) 
# with solid tumors

list_of_reasons_c2 <- c("Review",
                     "Retrospective",
                     "Design",
                     "Not Randomized",
                     "Case Study",
                     "Commentary",
                     "Non-randomized (post-hoc analysis)",
                     "Opinion",
                     "Economic",
                     "Non-randomized",
                     "Pooled Checkmate 57/17",
                     "FDA Approval Summary",
                     "First In-Human",
                     "Subgrop Analysis",
                     "subgrop Analysis",
                     "Protocol",
                     "Liquid tumor",
                     "Immunodisease positive subpop analysis",
                     "Not PD-L1 focused")
d$C2 <- as.numeric(d$Reason %in% list_of_reasons_c2)
n2 <- sum(d[d$C1==0,"C2"])


# Criteria 3: did not compare a PD-1/L1 inhibitor in monotherapy to a 
# standard cytotoxic agent

list_of_reasons_c3 <- c("Not PD-L1 specific therapy",
                        "Not monotherapy")
d$C3 <- as.numeric(d$Reason %in% list_of_reasons_c3)
n3 <- sum(d[d$C1==0 & d$C2==0,"C3"])


# Criteria 4: did not report Kaplan-Meier curves for overall survival 
# stratified by PD-L1 expression

list_of_reasons_c4 <- c("No OS",
                        "No OS data stratified by PD-L1",
                        "No PD-L1 KM",
                        "Uses H-score",
                        "Ongoing/No results"
)
d$C4 <- as.numeric(d$Reason %in% list_of_reasons_c4)
n4 <- sum(d[d$C1==0 & d$C2==0 & d$C3==0, "C4"])


# Criteria 5: reported data available from another publication

list_of_reasons_c5 <- c("Duplicate (data already in 29884413)")
d$C5 <- as.numeric(d$Reason %in% list_of_reasons_c5)
n5 <- sum(d[d$C1==0 & d$C2==0 & d$C3==0 & d$C4==0, "C5"])


# Final number of included trials
nfinal <- n0-n1-n2-n3-n4-n5
nfinal

# List Pubmed IDs of inclded trials
PUBMEDID_selected <- d[d$C1==0 & d$C2==0 & d$C3==0 & d$C4==0 & d$C5==0, "PMID"]

# Collects the results
flowchart <- data.frame(N=c(n0,n1,n2,n3,n4,n5,nfinal),
                        Reason=c(
  "# Number of identified abstracts published before November 31, 2018",
  "# Criteria 1: manuscript unaccessible by the authors",
  "# Criteria 2: not randomized phase II or III studies of PD-1/L1 inihibitors in adults (≥18 years) with solid tumors",
  "# Criteria 3: did not compare a PD-1/L1 inhibitor in monotherapy to a standard cytotoxic agent)",
  "# Criteria 4: did not report Kaplan-Meier curves for overall survival stratified by PD-L1 expression",
  "# Criteria 5: reported data available from another publication",
  "# Final number of included trials"))

sink(file="./RESULTS/flowchart_data.txt")
print("Flowchart data")
print(flowchart)
print("PUBMED ID of included trials")
print(PUBMEDID_selected)
sink()
