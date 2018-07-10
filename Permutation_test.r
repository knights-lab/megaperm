
## New permutation edition
outcome = map$Supplement
outcome[map$Supplement=="MCT" & map$StudyDayNo <= 8] = "preMCT"
outcome[map$Supplement=="EVOO" & map$StudyDayNo <= 8] = "preEVOO"
outcome = factor(outcome)
comp = c("MCT","EVOO","preMCT","preEVOO")
names(comp) = c("preMCT","preEVOO","MCT","EVOO")

nstrains = nrow(strains)
niter = 999999

#st.t = t(strains)
# Get baseline
#diffs = numeric(nstrains)
#for (i in which(outcome=="MCT"))
#  for (j in which(outcome=="EVOO"))
#    diffs = diffs + strains[,i] - strains[,j]
whichMCT = which(outcome=="MCT"); ml = length(whichMCT)
whichEVOO = which(outcome=="EVOO"); el = length(whichEVOO)
diffs = rowSums(strains[,whichMCT])*el - 
  rowSums(strains[,whichEVOO])*ml
absdiffs = abs(diffs*(1/(ml*el)))

sigdf = integer(nstrains)
pre.labels = as.character(outcome[map$StudyDayNo <= 8])
ppl = unique(map$UserName)
nppl = length(ppl)
testlab = character(length(outcome))
states = map[!duplicated(map$UserName),]$Supplement
pops = table(map$UserName)
pop_pre = tapply(outcome,map$UserName,function(x) sum(x=="preMCT" | x=="preEVOO"))
pop_pst = pops-pop_pre
pop_ix = integer(nppl+1)
pop_ix[1] = 1; pop_ix[nppl+1]=length(outcome)+1
for (i in 2:nppl) pop_ix[i] = pop_ix[i-1]+pops[i-1]
nppl_iter = 1:nppl

t0 = Sys.time()
for (i in 1:niter) {
  state.s = sample(states)
  state.pr = unname(comp[state.s])
  for (p in nppl_iter) 
    testlab[pop_ix[p]:(pop_ix[p+1]-1)] = 
      sample(c(rep(state.s[p],pop_pst[p]),rep(state.pr[p],pop_pre[p])))
  whichMCT = which(testlab=="MCT"); ml = length(whichMCT)
  whichEVOO = which(testlab=="EVOO"); el = length(whichEVOO)
  test = rowSums(strains[,whichMCT])*el - 
    rowSums(strains[,whichEVOO])*ml
  sigdf = sigdf + (abs(test*(1/(el*ml))) >= absdiffs)
}
print(Sys.time()-t0)
pvals = (sigdf+1) / (niter+1)
pvals.adj = p.adjust(pvals,method = "BH")
sort(pvals.adj)[1:10]