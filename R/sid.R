# count the ordered pairs (i, j) for which the parents of i in learned are not
# a valid adjustment set, in true, for the causal effect of i on j (the
# generalized adjustment criterion of Shpitser et al. / Perkovic et al.).
sid.dag.vs.dag = function(learned, true, debug = FALSE) {

  .Call(call_sid,
        learned = learned,
        true = true,
        debug = debug)

}#SID.DAG.VS.DAG
