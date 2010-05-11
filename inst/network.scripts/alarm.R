n = 20000

BOOL = c("TRUE", "FALSE")
LV2 = c("NORMAL", "HIGH")
LV3 = c("LOW", "NORMAL", "HIGH")
LV4 = c("ZERO", "LOW", "NORMAL", "HIGH")

HYPOVOLEMIA = sample(BOOL, n, prob = c(0.2, 0.8), replace = TRUE)
LVFAILURE = sample(BOOL, n, prob = c(0.05, 0.95), replace = TRUE)

HISTORY = LVFAILURE
HISTORY[HISTORY == "TRUE"] = sample(BOOL, length(which(HISTORY == "TRUE")), prob = c(0.9, 0.1), replace = TRUE)
HISTORY[HISTORY == "FALSE"] = sample(BOOL, length(which(HISTORY == "FALSE")), prob = c(0.01, 0.99), replace = TRUE)

LVEDVOLUME = apply(cbind(HYPOVOLEMIA, LVFAILURE), 1, paste, collapse = ":")
LVEDVOLUME[LVEDVOLUME == "TRUE:TRUE"] = sample(LV3, length(which(LVEDVOLUME == "TRUE:TRUE")), prob = c(0.95, 0.04, 0.01), replace = TRUE)
LVEDVOLUME[LVEDVOLUME == "FALSE:TRUE"] = sample(LV3, length(which(LVEDVOLUME == "FALSE:TRUE")), prob = c(0.98, 0.01, 0.01), replace = TRUE)
LVEDVOLUME[LVEDVOLUME == "TRUE:FALSE"] = sample(LV3, length(which(LVEDVOLUME == "TRUE:FALSE")), prob = c(0.01, 0.09, 0.9), replace = TRUE)
LVEDVOLUME[LVEDVOLUME == "FALSE:FALSE"] = sample(LV3, length(which(LVEDVOLUME == "FALSE:FALSE")), prob = c(0.05, 0.9, 0.05), replace = TRUE)

STROKEVOLUME = apply(cbind(HYPOVOLEMIA, LVFAILURE), 1, paste, collapse = ":")
STROKEVOLUME[STROKEVOLUME == "TRUE:TRUE"] = sample(LV3, length(which(STROKEVOLUME == "TRUE:TRUE")), prob = c(0.98, 0.01, 0.01), replace = TRUE)
STROKEVOLUME[STROKEVOLUME == "FALSE:TRUE"] = sample(LV3, length(which(STROKEVOLUME == "FALSE:TRUE")), prob = c(0.95, 0.04, 0.01), replace = TRUE)
STROKEVOLUME[STROKEVOLUME == "TRUE:FALSE"] = sample(LV3, length(which(STROKEVOLUME == "TRUE:FALSE")), prob = c(0.5, 0.49, 0.01), replace = TRUE)
STROKEVOLUME[STROKEVOLUME == "FALSE:FALSE"] = sample(LV3, length(which(STROKEVOLUME == "FALSE:FALSE")), prob = c(0.05, 0.9, 0.05), replace = TRUE)

CVP = LVEDVOLUME
CVP[CVP == "LOW"] = sample(LV3, length(which(CVP == "LOW")), prob = c(0.95, 0.04, 0.01), replace = TRUE)
CVP[CVP == "NORMAL"] = sample(LV3, length(which(CVP == "NORMAL")), prob = c(0.04, 0.95, 0.01), replace = TRUE)
CVP[CVP == "HIGH"] = sample(LV3, length(which(CVP == "HIGH")), prob = c(0.01, 0.29, 0.7), replace = TRUE)

PCWP = LVEDVOLUME
PCWP[PCWP == "LOW"] = sample(LV3, length(which(PCWP == "LOW")), prob = c(0.95, 0.04, 0.01), replace = TRUE)
PCWP[PCWP == "NORMAL"] = sample(LV3, length(which(PCWP == "NORMAL")), prob = c(0.04, 0.95, 0.01), replace = TRUE)
PCWP[PCWP == "HIGH"] = sample(LV3, length(which(PCWP == "HIGH")), prob = c(0.01, 0.04, 0.95), replace = TRUE)

ERRLOWOUTPUT = sample(BOOL, n, prob = c(0.05, 0.95), replace = TRUE)

ERRCAUTER = sample(BOOL, n, prob = c(0.1, 0.9), replace = TRUE)
INSUFFANESTH = sample(BOOL, n, prob = c(0.1, 0.9), replace = TRUE)
ANAPHYLAXIS = sample(BOOL, n, prob = c(0.01, 0.99), replace = TRUE)

TPR = ANAPHYLAXIS
TPR[TPR == "TRUE"] = sample(LV3, length(which(TPR == "TRUE")), prob = c(0.98, 0.01, 0.01), replace = TRUE)
TPR[TPR == "FALSE"] = sample(LV3, length(which(TPR == "FALSE")), prob = c(0.3, 0.4, 0.3), replace = TRUE)

KINKEDTUBE = sample(BOOL, n, prob = c(0.04, 0.96), replace = TRUE)
FIO2 = sample(c("LOW", "NORMAL"), n, prob = c(0.05, 0.95), replace = TRUE)
PULMEMBOLUS = sample(BOOL, n, prob = c(0.01, 0.99), replace = TRUE)

PAP = PULMEMBOLUS
PAP[PAP == "TRUE"] = sample(LV3, length(which(PAP == "TRUE")), prob = c(0.01, 0.19, 0.8), replace = TRUE)
PAP[PAP == "FALSE"] = sample(LV3, length(which(PAP == "FALSE")), prob = c(0.05, 0.9, 0.05), replace = TRUE)

INTUBATION = sample(c("NORMAL", "ESOPHAGEAL", "ONESIDED"), n, prob = c(0.92, 0.03, 0.05), replace = TRUE)

SHUNT = apply(cbind(INTUBATION, PULMEMBOLUS), 1, paste, collapse = ":")
SHUNT[SHUNT == "NORMAL:TRUE"] = sample(c("NORMAL", "HIGH"), length(which(SHUNT == "NORMAL:TRUE")), prob = c(0.1, 0.9), replace = TRUE)
SHUNT[SHUNT == "ESOPHAGEAL:TRUE"] = sample(c("NORMAL", "HIGH"), length(which(SHUNT == "ESOPHAGEAL:TRUE")), prob = c(0.1, 0.9), replace = TRUE)
SHUNT[SHUNT == "ONESIDED:TRUE"] = sample(c("NORMAL", "HIGH"), length(which(SHUNT == "ONESIDED:TRUE")), prob = c(0.01, 0.99), replace = TRUE)
SHUNT[SHUNT == "NORMAL:FALSE"] = sample(c("NORMAL", "HIGH"), length(which(SHUNT == "NORMAL:FALSE")), prob = c(0.95, 0.05), replace = TRUE)
SHUNT[SHUNT == "ESOPHAGEAL:FALSE"] = sample(c("NORMAL", "HIGH"), length(which(SHUNT == "ESOPHAGEAL:FALSE")), prob = c(0.95, 0.05), replace = TRUE)
SHUNT[SHUNT == "ONESIDED:FALSE"] = sample(c("NORMAL", "HIGH"), length(which(SHUNT == "ONESIDED:FALSE")), prob = c(0.05, 0.95), replace = TRUE)

DISCONNECT = sample(BOOL, n, prob = c(0.1, 0.9), replace = TRUE)
MINVOLSET = sample(LV3, n, prob = c(0.05, 0.9, 0.05), replace = TRUE)

VENTMACH = MINVOLSET
VENTMACH[VENTMACH == "LOW"] = sample(LV4, length(which(VENTMACH == "LOW")), prob = c(0.05, 0.93, 0.01, 0.01), replace = TRUE)
VENTMACH[VENTMACH == "NORMAL"] = sample(LV4, length(which(VENTMACH == "NORMAL")), prob = c(0.05, 0.01, 0.93, 0.01), replace = TRUE)
VENTMACH[VENTMACH == "HIGH"] = sample(LV4, length(which(VENTMACH == "HIGH")), prob = c(0.05, 0.01, 0.01, 0.93), replace = TRUE)

VENTTUBE = apply(cbind(DISCONNECT, VENTMACH), 1, paste, collapse = ":")
VENTTUBE[VENTTUBE == "TRUE:ZERO"] = sample(LV4, length(which(VENTTUBE == "TRUE:ZERO")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTTUBE[VENTTUBE == "FALSE:ZERO"] = sample(LV4, length(which(VENTTUBE == "FALSE:ZERO")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTTUBE[VENTTUBE == "TRUE:LOW"] = sample(LV4, length(which(VENTTUBE == "TRUE:LOW")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTTUBE[VENTTUBE == "FALSE:LOW"] = sample(LV4, length(which(VENTTUBE == "FALSE:LOW")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTTUBE[VENTTUBE == "TRUE:NORMAL"] = sample(LV4, length(which(VENTTUBE == "TRUE:NORMAL")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTTUBE[VENTTUBE == "FALSE:NORMAL"] = sample(LV4, length(which(VENTTUBE == "FALSE:NORMAL")), prob = c(0.01, 0.97, 0.01, 0.01), replace = TRUE)
VENTTUBE[VENTTUBE == "TRUE:HIGH"] = sample(LV4, length(which(VENTTUBE == "TRUE:HIGH")), prob = c(0.01, 0.01, 0.97, 0.01), replace = TRUE)
VENTTUBE[VENTTUBE == "FALSE:HIGH"] = sample(LV4, length(which(VENTTUBE == "FALSE:HIGH")), prob = c(0.01, 0.01, 0.01, 0.97), replace = TRUE)

PRESS = apply(cbind(INTUBATION, KINKEDTUBE, VENTTUBE), 1, paste, collapse = ":")
PRESS[PRESS == "NORMAL:TRUE:ZERO"] = sample(LV4, length(which(PRESS == "NORMAL:TRUE:ZERO")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
PRESS[PRESS == "ESOPHAGEAL:TRUE:ZERO"] = sample(LV4, length(which(PRESS == "ESOPHAGEAL:TRUE:ZERO")), prob = c(0.01, 0.30, 0.49, 0.20), replace = TRUE)
PRESS[PRESS == "ONESIDED:TRUE:ZERO"] = sample(LV4, length(which(PRESS == "ONESIDED:TRUE:ZERO")), prob = c(0.01, 0.01, 0.08, 0.90), replace = TRUE)
PRESS[PRESS == "NORMAL:FALSE:ZERO"] = sample(LV4, length(which(PRESS == "NORMAL:FALSE:ZERO")), prob = c(0.01, 0.01, 0.01, 0.97), replace = TRUE)
PRESS[PRESS == "ESOPHAGEAL:FALSE:ZERO"] = sample(LV4, length(which(PRESS == "ESOPHAGEAL:FALSE:ZERO")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
PRESS[PRESS == "ONESIDED:FALSE:ZERO"] = sample(LV4, length(which(PRESS == "ONESIDED:FALSE:ZERO")), prob = c(0.10, 0.84, 0.05, 0.01), replace = TRUE)
PRESS[PRESS == "NORMAL:TRUE:LOW"] = sample(LV4, length(which(PRESS == "NORMAL:TRUE:LOW")), prob = c(0.05, 0.25, 0.25, 0.45), replace = TRUE)
PRESS[PRESS == "ESOPHAGEAL:TRUE:LOW"] = sample(LV4, length(which(PRESS == "ESOPHAGEAL:TRUE:LOW")), prob = c(0.01, 0.15, 0.25, 0.59), replace = TRUE)
PRESS[PRESS == "ONESIDED:TRUE:LOW"] = sample(LV4, length(which(PRESS == "ONESIDED:TRUE:LOW")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
PRESS[PRESS == "NORMAL:FALSE:LOW"] = sample(LV4, length(which(PRESS == "NORMAL:FALSE:LOW")), prob = c(0.01, 0.29, 0.30, 0.40), replace = TRUE)
PRESS[PRESS == "ESOPHAGEAL:FALSE:LOW"] =  sample(LV4, length(which(PRESS == "ESOPHAGEAL:FALSE:LOW")), prob = c(0.01, 0.01, 0.08, 0.90), replace = TRUE)
PRESS[PRESS == "ONESIDED:FALSE:LOW"] =  sample(LV4, length(which(PRESS == "ONESIDED:FALSE:LOW")), prob = c(0.01, 0.01, 0.01, 0.97), replace = TRUE)
PRESS[PRESS == "NORMAL:TRUE:NORMAL"] =  sample(LV4, length(which(PRESS == "NORMAL:TRUE:NORMAL")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
PRESS[PRESS == "ESOPHAGEAL:TRUE:NORMAL"] =  sample(LV4, length(which(PRESS == "ESOPHAGEAL:TRUE:NORMAL")), prob = c(0.01, 0.97, 0.01, 0.01), replace = TRUE)
PRESS[PRESS == "ONESIDED:TRUE:NORMAL"] =  sample(LV4, length(which(PRESS == "ONESIDED:TRUE:NORMAL")), prob = c(0.01, 0.01, 0.97, 0.01), replace = TRUE)
PRESS[PRESS == "NORMAL:FALSE:NORMAL"] =  sample(LV4, length(which(PRESS == "NORMAL:FALSE:NORMAL")), prob = c(0.01, 0.01, 0.01, 0.97), replace = TRUE)
PRESS[PRESS == "ESOPHAGEAL:FALSE:NORMAL"] =  sample(LV4, length(which(PRESS == "ESOPHAGEAL:FALSE:NORMAL")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
PRESS[PRESS == "ONESIDED:FALSE:NORMAL"] =  sample(LV4, length(which(PRESS == "ONESIDED:FALSE:NORMAL")), prob = c(0.40, 0.58, 0.01, 0.01), replace = TRUE)
PRESS[PRESS == "NORMAL:TRUE:HIGH"] =  sample(LV4, length(which(PRESS == "NORMAL:TRUE:HIGH")), prob = c(0.20, 0.75, 0.04, 0.01), replace = TRUE)
PRESS[PRESS == "ESOPHAGEAL:TRUE:HIGH"] =  sample(LV4, length(which(PRESS == "ESOPHAGEAL:TRUE:HIGH")), prob = c(0.20, 0.70, 0.09, 0.01), replace = TRUE)
PRESS[PRESS == "ONESIDED:TRUE:HIGH"] =  sample(LV4, length(which(PRESS == "ONESIDED:TRUE:HIGH")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
PRESS[PRESS == "NORMAL:FALSE:HIGH"] =  sample(LV4, length(which(PRESS == "NORMAL:FALSE:HIGH")), prob = c(0.01, 0.90, 0.08, 0.01), replace = TRUE)
PRESS[PRESS == "ESOPHAGEAL:FALSE:HIGH"] =  sample(LV4, length(which(PRESS == "ESOPHAGEAL:FALSE:HIGH")), prob = c(0.01, 0.01, 0.38, 0.60), replace = TRUE)
PRESS[PRESS == "ONESIDED:FALSE:HIGH"] =  sample(LV4, length(which(PRESS == "ONESIDED:FALSE:HIGH")), prob = c(0.01, 0.01, 0.01, 0.97), replace = TRUE)

VENTLUNG = apply(cbind(INTUBATION, KINKEDTUBE, VENTTUBE), 1, paste, collapse = ":")
VENTLUNG[VENTLUNG == "NORMAL:TRUE:ZERO"] = sample(LV4, length(which(VENTLUNG == "NORMAL:TRUE:ZERO")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "ESOPHAGEAL:TRUE:ZERO"] = sample(LV4, length(which(VENTLUNG == "ESOPHAGEAL:TRUE:ZERO")), prob = c(0.95, 0.03, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "ONESIDED:TRUE:ZERO"] = sample(LV4, length(which(VENTLUNG == "ONESIDED:TRUE:ZERO")), prob = c(0.4, 0.58, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "NORMAL:FALSE:ZERO"] = sample(LV4, length(which(VENTLUNG == "NORMAL:FALSE:ZERO")), prob = c(0.3, 0.68, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "ESOPHAGEAL:FALSE:ZERO"] = sample(LV4, length(which(VENTLUNG == "ESOPHAGEAL:FALSE:ZERO")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "ONESIDED:FALSE:ZERO"] = sample(LV4, length(which(VENTLUNG == "ONESIDED:FALSE:ZERO")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "NORMAL:TRUE:LOW"] = sample(LV4, length(which(VENTLUNG == "NORMAL:TRUE:LOW")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "ESOPHAGEAL:TRUE:LOW"] = sample(LV4, length(which(VENTLUNG == "ESOPHAGEAL:TRUE:LOW")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "ONESIDED:TRUE:LOW"] = sample(LV4, length(which(VENTLUNG == "ONESIDED:TRUE:LOW")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "NORMAL:FALSE:LOW"] = sample(LV4, length(which(VENTLUNG == "NORMAL:FALSE:LOW")), prob = c(0.95, 0.03, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "ESOPHAGEAL:FALSE:LOW"] =  sample(LV4, length(which(VENTLUNG == "ESOPHAGEAL:FALSE:LOW")), prob = c(0.5, 0.48, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "ONESIDED:FALSE:LOW"] =  sample(LV4, length(which(VENTLUNG == "ONESIDED:FALSE:LOW")), prob = c(0.3, 0.68, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "NORMAL:TRUE:NORMAL"] =  sample(LV4, length(which(VENTLUNG == "NORMAL:TRUE:NORMAL")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "ESOPHAGEAL:TRUE:NORMAL"] =  sample(LV4, length(which(VENTLUNG == "ESOPHAGEAL:TRUE:NORMAL")), prob = c(0.01, 0.97, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "ONESIDED:TRUE:NORMAL"] =  sample(LV4, length(which(VENTLUNG == "ONESIDED:TRUE:NORMAL")), prob = c(0.01, 0.01, 0.97, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "NORMAL:FALSE:NORMAL"] =  sample(LV4, length(which(VENTLUNG == "NORMAL:FALSE:NORMAL")), prob = c(0.01, 0.01, 0.01, 0.97), replace = TRUE)
VENTLUNG[VENTLUNG == "ESOPHAGEAL:FALSE:NORMAL"] =  sample(LV4, length(which(VENTLUNG == "ESOPHAGEAL:FALSE:NORMAL")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "ONESIDED:FALSE:NORMAL"] =  sample(LV4, length(which(VENTLUNG == "ONESIDED:FALSE:NORMAL")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "NORMAL:TRUE:HIGH"] =  sample(LV4, length(which(VENTLUNG == "NORMAL:TRUE:HIGH")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "ESOPHAGEAL:TRUE:HIGH"] =  sample(LV4, length(which(VENTLUNG == "ESOPHAGEAL:TRUE:HIGH")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "ONESIDED:TRUE:HIGH"] =  sample(LV4, length(which(VENTLUNG == "ONESIDED:TRUE:HIGH")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "NORMAL:FALSE:HIGH"] =  sample(LV4, length(which(VENTLUNG == "NORMAL:FALSE:HIGH")), prob = c(0.01, 0.97, 0.01, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "ESOPHAGEAL:FALSE:HIGH"] =  sample(LV4, length(which(VENTLUNG == "ESOPHAGEAL:FALSE:HIGH")), prob = c(0.01, 0.01, 0.97, 0.01), replace = TRUE)
VENTLUNG[VENTLUNG == "ONESIDED:FALSE:HIGH"] =  sample(LV4, length(which(VENTLUNG == "ONESIDED:FALSE:HIGH")), prob = c(0.01, 0.01, 0.01, 0.97), replace = TRUE)

MINVOL = apply(cbind(INTUBATION, VENTLUNG), 1, paste, collapse = ":")
MINVOL[MINVOL == "NORMAL:ZERO"] = sample(LV4, length(which(MINVOL == "NORMAL:ZERO")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
MINVOL[MINVOL == "ESOPHAGEAL:ZERO"] = sample(LV4, length(which(MINVOL == "ESOPHAGEAL:ZERO")), prob = c(0.01, 0.97, 0.01, 0.01), replace = TRUE)
MINVOL[MINVOL == "ONESIDED:ZERO"] = sample(LV4, length(which(MINVOL == "ONESIDED:ZERO")), prob = c(0.01, 0.01, 0.97, 0.01), replace = TRUE)
MINVOL[MINVOL == "NORMAL:LOW"] = sample(LV4, length(which(MINVOL == "NORMAL:LOW")), prob = c(0.01, 0.01, 0.01, 0.97), replace = TRUE)
MINVOL[MINVOL == "ESOPHAGEAL:LOW"] = sample(LV4, length(which(MINVOL == "ESOPHAGEAL:LOW")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
MINVOL[MINVOL == "ONESIDED:LOW"] = sample(LV4, length(which(MINVOL == "ONESIDED:LOW")), prob = c(0.60, 0.38, 0.01, 0.01), replace = TRUE)
MINVOL[MINVOL == "NORMAL:NORMAL"] = sample(LV4, length(which(MINVOL == "NORMAL:NORMAL")), prob = c(0.50, 0.48, 0.01, 0.01), replace = TRUE)
MINVOL[MINVOL == "ESOPHAGEAL:NORMAL"] = sample(LV4, length(which(MINVOL == "ESOPHAGEAL:NORMAL")), prob = c(0.50, 0.48, 0.01, 0.01), replace = TRUE)
MINVOL[MINVOL == "ONESIDED:NORMAL"] = sample(LV4, length(which(MINVOL == "ONESIDED:NORMAL")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
MINVOL[MINVOL == "NORMAL:HIGH"] = sample(LV4, length(which(MINVOL == "NORMAL:HIGH")), prob = c(0.01, 0.97, 0.01, 0.01), replace = TRUE)
MINVOL[MINVOL == "ESOPHAGEAL:HIGH"] = sample(LV4, length(which(MINVOL == "ESOPHAGEAL:HIGH")), prob = c(0.01, 0.01, 0.97, 0.01), replace = TRUE)
MINVOL[MINVOL == "ONESIDED:HIGH"] = sample(LV4, length(which(MINVOL == "ONESIDED:HIGH")), prob = c(0.01, 0.01, 0.01, 0.97), replace = TRUE)

VENTALV = apply(cbind(INTUBATION, VENTLUNG), 1, paste, collapse = ":")
VENTALV[VENTALV == "NORMAL:ZERO"] = sample(LV4, length(which(VENTALV == "NORMAL:ZERO")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTALV[VENTALV == "ESOPHAGEAL:ZERO"] = sample(LV4, length(which(VENTALV == "ESOPHAGEAL:ZERO")), prob = c(0.01, 0.97, 0.01, 0.01), replace = TRUE)
VENTALV[VENTALV == "ONESIDED:ZERO"] = sample(LV4, length(which(VENTALV == "ONESIDED:ZERO")), prob = c(0.01, 0.01, 0.97, 0.01), replace = TRUE)
VENTALV[VENTALV == "NORMAL:LOW"] = sample(LV4, length(which(VENTALV == "NORMAL:LOW")), prob = c(0.01, 0.01, 0.01, 0.97), replace = TRUE)
VENTALV[VENTALV == "ESOPHAGEAL:LOW"] = sample(LV4, length(which(VENTALV == "ESOPHAGEAL:LOW")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTALV[VENTALV == "ONESIDED:LOW"] = sample(LV4, length(which(VENTALV == "ONESIDED:LOW")), prob = c(0.01, 0.97, 0.01, 0.01), replace = TRUE)
VENTALV[VENTALV == "NORMAL:NORMAL"] = sample(LV4, length(which(VENTALV == "NORMAL:NORMAL")), prob = c(0.01, 0.01, 0.97, 0.01), replace = TRUE)
VENTALV[VENTALV == "ESOPHAGEAL:NORMAL"] = sample(LV4, length(which(VENTALV == "ESOPHAGEAL:NORMAL")), prob = c(0.01, 0.01, 0.01, 0.97), replace = TRUE)
VENTALV[VENTALV == "ONESIDED:NORMAL"] = sample(LV4, length(which(VENTALV == "ONESIDED:NORMAL")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
VENTALV[VENTALV == "NORMAL:HIGH"] = sample(LV4, length(which(VENTALV == "NORMAL:HIGH")), prob = c(0.03, 0.95, 0.01, 0.01), replace = TRUE)
VENTALV[VENTALV == "ESOPHAGEAL:HIGH"] = sample(LV4, length(which(VENTALV == "ESOPHAGEAL:HIGH")), prob = c(0.01, 0.94, 0.04, 0.01), replace = TRUE)
VENTALV[VENTALV == "ONESIDED:HIGH"] = sample(LV4, length(which(VENTALV == "ONESIDED:HIGH")), prob = c(0.01, 0.88, 0.10, 0.01), replace = TRUE)

ARTCO2 = VENTALV
ARTCO2[ARTCO2 == "ZERO"] = sample(LV3, length(which(ARTCO2 == "ZERO")), prob = c(0.01, 0.01, 0.98), replace = TRUE)
ARTCO2[ARTCO2 == "LOW"] = sample(LV3, length(which(ARTCO2 == "LOW")), prob = c(0.01, 0.01, 0.98), replace = TRUE)
ARTCO2[ARTCO2 == "NORMAL"] = sample(LV3, length(which(ARTCO2 == "NORMAL")), prob = c(0.04, 0.92, 0.04), replace = TRUE)
ARTCO2[ARTCO2 == "HIGH"] = sample(LV3, length(which(ARTCO2 == "HIGH")), prob = c(0.90, 0.09, 0.01), replace = TRUE)

EXPCO2 =  apply(cbind(ARTCO2, VENTLUNG), 1, paste, collapse = ":")
EXPCO2[EXPCO2 == "LOW:ZERO"] = sample(LV4, length(which(EXPCO2 == "LOW:ZERO")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
EXPCO2[EXPCO2 == "NORMAL:ZERO"] = sample(LV4, length(which(EXPCO2 == "NORMAL:ZERO")), prob = c(0.01, 0.97, 0.01, 0.01), replace = TRUE)
EXPCO2[EXPCO2 == "HIGH:ZERO"] = sample(LV4, length(which(EXPCO2 == "HIGH:ZERO")), prob = c(0.01, 0.97, 0.01, 0.01), replace = TRUE)
EXPCO2[EXPCO2 == "LOW:LOW"] = sample(LV4, length(which(EXPCO2 == "LOW:LOW")), prob = c(0.01, 0.97, 0.01, 0.01), replace = TRUE)
EXPCO2[EXPCO2 == "NORMAL:LOW"] = sample(LV4, length(which(EXPCO2 == "NORMAL:LOW")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
EXPCO2[EXPCO2 == "HIGH:LOW"] = sample(LV4, length(which(EXPCO2 == "HIGH:LOW")), prob = c(0.01, 0.01, 0.97, 0.01), replace = TRUE)
EXPCO2[EXPCO2 == "LOW:NORMAL"] = sample(LV4, length(which(EXPCO2 == "LOW:NORMAL")), prob = c(0.01, 0.01, 0.97, 0.01), replace = TRUE)
EXPCO2[EXPCO2 == "NORMAL:NORMAL"] = sample(LV4, length(which(EXPCO2 == "NORMAL:NORMAL")), prob = c(0.01, 0.01, 0.97, 0.01), replace = TRUE)
EXPCO2[EXPCO2 == "HIGH:NORMAL"] = sample(LV4, length(which(EXPCO2 == "HIGH:NORMAL")), prob = c(0.97, 0.01, 0.01, 0.01), replace = TRUE)
EXPCO2[EXPCO2 == "LOW:HIGH"] = sample(LV4, length(which(EXPCO2 == "LOW:HIGH")), prob = c(0.01, 0.01, 0.01, 0.97), replace = TRUE)
EXPCO2[EXPCO2 == "NORMAL:HIGH"] = sample(LV4, length(which(EXPCO2 == "NORMAL:HIGH")), prob = c(0.01, 0.01, 0.01, 0.97), replace = TRUE)
EXPCO2[EXPCO2 == "HIGH:HIGH"] = sample(LV4, length(which(EXPCO2 == "HIGH:HIGH")), prob = c(0.01, 0.01, 0.01, 0.97), replace = TRUE)

PVSAT =  apply(cbind(FIO2, VENTALV), 1, paste, collapse = ":")
PVSAT[PVSAT == "LOW:ZERO"] = sample(LV3, length(which(PVSAT == "LOW:ZERO")), prob = c(1, 0, 0), replace = TRUE)
PVSAT[PVSAT == "NORMAL:ZERO"] = sample(LV3, length(which(PVSAT == "NORMAL:ZERO")), prob = c(0.99, 0.01, 0), replace = TRUE)
PVSAT[PVSAT == "LOW:LOW"] = sample(LV3, length(which(PVSAT == "LOW:LOW")), prob = c(0.95, 0.04, 0.01), replace = TRUE)
PVSAT[PVSAT == "NORMAL:LOW"] = sample(LV3, length(which(PVSAT == "NORMAL:LOW")), prob = c(0.95, 0.04, 0.01), replace = TRUE)
PVSAT[PVSAT == "LOW:NORMAL"] = sample(LV3, length(which(PVSAT == "LOW:NORMAL")), prob = c(1, 0, 0), replace = TRUE)
PVSAT[PVSAT == "NORMAL:NORMAL"] = sample(LV3, length(which(PVSAT == "NORMAL:NORMAL")), prob = c(0.95, 0.04, 0.01), replace = TRUE)
PVSAT[PVSAT == "LOW:HIGH"] = sample(LV3, length(which(PVSAT == "LOW:HIGH")), prob = c(0.01, 0.95, 0.04), replace = TRUE)
PVSAT[PVSAT == "NORMAL:HIGH"] = sample(LV3, length(which(PVSAT == "NORMAL:HIGH")), prob = c(0.01, 0.01, 0.98), replace = TRUE)

SAO2 = apply(cbind(PVSAT, SHUNT), 1, paste, collapse = ":")
SAO2[SAO2 == "LOW:NORMAL"] = sample(LV3, length(which(SAO2 == "LOW:NORMAL")), prob = c(0.98, 0.01, 0.01), replace = TRUE)
SAO2[SAO2 == "NORMAL:NORMAL"] = sample(LV3, length(which(SAO2 == "NORMAL:NORMAL")), prob = c(0.01, 0.98, 0.01), replace = TRUE)
SAO2[SAO2 == "HIGH:NORMAL"] = sample(LV3, length(which(SAO2 == "HIGH:NORMAL")), prob = c(0.01, 0.01, 0.98), replace = TRUE)
SAO2[SAO2 == "LOW:HIGH"] = sample(LV3, length(which(SAO2 == "LOW:HIGH")), prob = c(0.98, 0.01, 0.01), replace = TRUE)
SAO2[SAO2 == "NORMAL:HIGH"] = sample(LV3, length(which(SAO2 == "NORMAL:HIGH")), prob = c(0.98, 0.01, 0.01), replace = TRUE)
SAO2[SAO2 == "HIGH:HIGH"] = sample(LV3, length(which(SAO2 == "HIGH:HIGH")), prob = c(0.69, 0.3, 0.01), replace = TRUE)

CATECHOL = apply(cbind(ARTCO2, INSUFFANESTH, SAO2, TPR), 1, paste, collapse = ":")
CATECHOL[CATECHOL == "LOW:TRUE:LOW:LOW"] = sample(LV2, length(which(CATECHOL == "LOW:TRUE:LOW:LOW")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "NORMAL:TRUE:LOW:LOW"] = sample(LV2, length(which(CATECHOL == "NORMAL:TRUE:LOW:LOW")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "HIGH:TRUE:LOW:LOW"] = sample(LV2, length(which(CATECHOL == "HIGH:TRUE:LOW:LOW")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "LOW:FALSE:LOW:LOW"] = sample(LV2, length(which(CATECHOL == "LOW:FALSE:LOW:LOW")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "NORMAL:FALSE:LOW:LOW"] = sample(LV2, length(which(CATECHOL == "NORMAL:FALSE:LOW:LOW")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "HIGH:FALSE:LOW:LOW"] = sample(LV2, length(which(CATECHOL == "HIGH:FALSE:LOW:LOW")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "LOW:TRUE:NORMAL:LOW"] = sample(LV2, length(which(CATECHOL == "LOW:TRUE:NORMAL:LOW")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "NORMAL:TRUE:NORMAL:LOW"] = sample(LV2, length(which(CATECHOL == "NORMAL:TRUE:NORMAL:LOW")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "HIGH:TRUE:NORMAL:LOW"] = sample(LV2, length(which(CATECHOL == "HIGH:TRUE:NORMAL:LOW")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "LOW:FALSE:NORMAL:LOW"] = sample(LV2, length(which(CATECHOL == "LOW:FALSE:NORMAL:LOW")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "NORMAL:FALSE:NORMAL:LOW"] = sample(LV2, length(which(CATECHOL == "NORMAL:FALSE:NORMAL:LOW")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "HIGH:FALSE:NORMAL:LOW"] = sample(LV2, length(which(CATECHOL == "HIGH:FALSE:NORMAL:LOW")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "LOW:TRUE:HIGH:LOW"] = sample(LV2, length(which(CATECHOL == "LOW:TRUE:HIGH:LOW")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "NORMAL:TRUE:HIGH:LOW"] = sample(LV2, length(which(CATECHOL == "NORMAL:TRUE:HIGH:LOW")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "HIGH:TRUE:HIGH:LOW"] = sample(LV2, length(which(CATECHOL == "HIGH:TRUE:HIGH:LOW")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "LOW:FALSE:HIGH:LOW"] = sample(LV2, length(which(CATECHOL == "LOW:FALSE:HIGH:LOW")), prob = c(0.05, 0.95), replace = TRUE)
CATECHOL[CATECHOL == "NORMAL:FALSE:HIGH:LOW"] = sample(LV2, length(which(CATECHOL == "NORMAL:FALSE:HIGH:LOW")), prob = c(0.05, 0.95), replace = TRUE)
CATECHOL[CATECHOL == "HIGH:FALSE:HIGH:LOW"] = sample(LV2, length(which(CATECHOL == "HIGH:FALSE:HIGH:LOW")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "LOW:TRUE:LOW:NORMAL"] = sample(LV2, length(which(CATECHOL == "LOW:TRUE:LOW:NORMAL")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "NORMAL:TRUE:LOW:NORMAL"] = sample(LV2, length(which(CATECHOL == "NORMAL:TRUE:LOW:NORMAL")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "HIGH:TRUE:LOW:NORMAL"] = sample(LV2, length(which(CATECHOL == "HIGH:TRUE:LOW:NORMAL")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "LOW:FALSE:LOW:NORMAL"] = sample(LV2, length(which(CATECHOL == "LOW:FALSE:LOW:NORMAL")), prob = c(0.05, 0.95), replace = TRUE)
CATECHOL[CATECHOL == "NORMAL:FALSE:LOW:NORMAL"] = sample(LV2, length(which(CATECHOL == "NORMAL:FALSE:LOW:NORMAL")), prob = c(0.05, 0.95), replace = TRUE)
CATECHOL[CATECHOL == "HIGH:FALSE:LOW:NORMAL"] = sample(LV2, length(which(CATECHOL == "HIGH:FALSE:LOW:NORMAL")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "LOW:TRUE:NORMAL:NORMAL"] = sample(LV2, length(which(CATECHOL == "LOW:TRUE:NORMAL:NORMAL")), prob = c(0.05, 0.95), replace = TRUE)
CATECHOL[CATECHOL == "NORMAL:TRUE:NORMAL:NORMAL"] = sample(LV2, length(which(CATECHOL == "NORMAL:TRUE:NORMAL:NORMAL")), prob = c(0.05, 0.95), replace = TRUE)
CATECHOL[CATECHOL == "HIGH:TRUE:NORMAL:NORMAL"] = sample(LV2, length(which(CATECHOL == "HIGH:TRUE:NORMAL:NORMAL")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "LOW:FALSE:NORMAL:NORMAL"] = sample(LV2, length(which(CATECHOL == "LOW:FALSE:NORMAL:NORMAL")), prob = c(0.05, 0.95), replace = TRUE)
CATECHOL[CATECHOL == "NORMAL:FALSE:NORMAL:NORMAL"] = sample(LV2, length(which(CATECHOL == "NORMAL:FALSE:NORMAL:NORMAL")), prob = c(0.05, 0.95), replace = TRUE)
CATECHOL[CATECHOL == "HIGH:FALSE:NORMAL:NORMAL"] = sample(LV2, length(which(CATECHOL == "HIGH:FALSE:NORMAL:NORMAL")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "LOW:TRUE:HIGH:NORMAL"] = sample(LV2, length(which(CATECHOL == "LOW:TRUE:HIGH:NORMAL")), prob = c(0.05, 0.95), replace = TRUE)
CATECHOL[CATECHOL == "NORMAL:TRUE:HIGH:NORMAL"] = sample(LV2, length(which(CATECHOL == "NORMAL:TRUE:HIGH:NORMAL")), prob = c(0.05, 0.95), replace = TRUE)
CATECHOL[CATECHOL == "HIGH:TRUE:HIGH:NORMAL"] = sample(LV2, length(which(CATECHOL == "HIGH:TRUE:HIGH:NORMAL")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "LOW:FALSE:HIGH:NORMAL"] = sample(LV2, length(which(CATECHOL == "LOW:FALSE:HIGH:NORMAL")), prob = c(0.05, 0.95), replace = TRUE)
CATECHOL[CATECHOL == "NORMAL:FALSE:HIGH:NORMAL"] = sample(LV2, length(which(CATECHOL == "NORMAL:FALSE:HIGH:NORMAL")), prob = c(0.05, 0.95), replace = TRUE)
CATECHOL[CATECHOL == "HIGH:FALSE:HIGH:NORMAL"] = sample(LV2, length(which(CATECHOL == "HIGH:FALSE:HIGH:NORMAL")), prob = c(0.01, 0.99), replace = TRUE)
CATECHOL[CATECHOL == "LOW:TRUE:LOW:HIGH"] = sample(LV2, length(which(CATECHOL == "LOW:TRUE:LOW:HIGH")), prob = c(0.70, 0.30), replace = TRUE)
CATECHOL[CATECHOL == "NORMAL:TRUE:LOW:HIGH"] = sample(LV2, length(which(CATECHOL == "NORMAL:TRUE:LOW:HIGH")), prob = c(0.70, 0.30), replace = TRUE)
CATECHOL[CATECHOL == "HIGH:TRUE:LOW:HIGH"] = sample(LV2, length(which(CATECHOL == "HIGH:TRUE:LOW:HIGH")), prob = c(0.10, 0.90), replace = TRUE)
CATECHOL[CATECHOL == "LOW:FALSE:LOW:HIGH"] = sample(LV2, length(which(CATECHOL == "LOW:FALSE:LOW:HIGH")), prob = c(0.70, 0.30), replace = TRUE)
CATECHOL[CATECHOL == "NORMAL:FALSE:LOW:HIGH"] = sample(LV2, length(which(CATECHOL == "NORMAL:FALSE:LOW:HIGH")), prob = c(0.70, 0.30), replace = TRUE)
CATECHOL[CATECHOL == "HIGH:FALSE:LOW:HIGH"] = sample(LV2, length(which(CATECHOL == "HIGH:FALSE:LOW:HIGH")), prob = c(0.10, 0.90), replace = TRUE)
CATECHOL[CATECHOL == "LOW:TRUE:NORMAL:HIGH"] = sample(LV2, length(which(CATECHOL == "LOW:TRUE:NORMAL:HIGH")), prob = c(0.70, 0.30), replace = TRUE)
CATECHOL[CATECHOL == "NORMAL:TRUE:NORMAL:HIGH"] = sample(LV2, length(which(CATECHOL == "NORMAL:TRUE:NORMAL:HIGH")), prob = c(0.70, 0.30), replace = TRUE)
CATECHOL[CATECHOL == "HIGH:TRUE:NORMAL:HIGH"] = sample(LV2, length(which(CATECHOL == "HIGH:TRUE:NORMAL:HIGH")), prob = c(0.10, 0.90), replace = TRUE)
CATECHOL[CATECHOL == "LOW:FALSE:NORMAL:HIGH"] = sample(LV2, length(which(CATECHOL == "LOW:FALSE:NORMAL:HIGH")), prob = c(0.95, 0.05), replace = TRUE)
CATECHOL[CATECHOL == "NORMAL:FALSE:NORMAL:HIGH"] = sample(LV2, length(which(CATECHOL == "NORMAL:FALSE:NORMAL:HIGH")), prob = c(0.99, 0.01), replace = TRUE)
CATECHOL[CATECHOL == "HIGH:FALSE:NORMAL:HIGH"] = sample(LV2, length(which(CATECHOL == "HIGH:FALSE:NORMAL:HIGH")), prob = c(0.30, 0.70), replace = TRUE)
CATECHOL[CATECHOL == "LOW:TRUE:HIGH:HIGH"] = sample(LV2, length(which(CATECHOL == "LOW:TRUE:HIGH:HIGH")), prob = c(0.95, 0.05), replace = TRUE)
CATECHOL[CATECHOL == "NORMAL:TRUE:HIGH:HIGH"] = sample(LV2, length(which(CATECHOL == "NORMAL:TRUE:HIGH:HIGH")), prob = c(0.99, 0.01), replace = TRUE)
CATECHOL[CATECHOL == "HIGH:TRUE:HIGH:HIGH"] = sample(LV2, length(which(CATECHOL == "HIGH:TRUE:HIGH:HIGH")), prob = c(0.30, 0.70), replace = TRUE)
CATECHOL[CATECHOL == "LOW:FALSE:HIGH:HIGH"] = sample(LV2, length(which(CATECHOL == "LOW:FALSE:HIGH:HIGH")), prob = c(0.95, 0.05), replace = TRUE)
CATECHOL[CATECHOL == "NORMAL:FALSE:HIGH:HIGH"] = sample(LV2, length(which(CATECHOL == "NORMAL:FALSE:HIGH:HIGH")), prob = c(0.99, 0.01), replace = TRUE)
CATECHOL[CATECHOL == "HIGH:FALSE:HIGH:HIGH"] = sample(LV2, length(which(CATECHOL == "HIGH:FALSE:HIGH:HIGH")), prob = c(0.30, 0.70), replace = TRUE)

HR = CATECHOL
HR[HR == "NORMAL"] = sample(LV3, length(which(HR == "NORMAL")), prob = c(0.05, 0.90, 0.05), replace = TRUE)
HR[HR == "HIGH"] = sample(LV3, length(which(HR == "HIGH")), prob = c(0.01, 0.09, 0.90), replace = TRUE)

HRBP = apply(cbind(ERRLOWOUTPUT, HR), 1, paste, collapse = ":")
HRBP[HRBP == "TRUE:LOW"] = sample(LV3, length(which(HRBP == "TRUE:LOW")), prob = c(0.98, 0.01, 0.01), replace = TRUE)
HRBP[HRBP == "FALSE:LOW"] = sample(LV3, length(which(HRBP == "FALSE:LOW")), prob = c(0.40, 0.59, 0.01), replace = TRUE)
HRBP[HRBP == "TRUE:NORMAL"] = sample(LV3, length(which(HRBP == "TRUE:NORMAL")), prob = c(0.30, 0.40, 0.30), replace = TRUE)
HRBP[HRBP == "FALSE:NORMAL"] = sample(LV3, length(which(HRBP == "FALSE:NORMAL")), prob = c(0.98, 0.01, 0.01), replace = TRUE)
HRBP[HRBP == "TRUE:HIGH"] = sample(LV3, length(which(HRBP == "TRUE:HIGH")), prob = c(0.01, 0.98, 0.01), replace = TRUE)
HRBP[HRBP == "FALSE:HIGH"] = sample(LV3, length(which(HRBP == "FALSE:HIGH")), prob = c(0.01, 0.01, 0.98), replace = TRUE)

HREKG = apply(cbind(ERRCAUTER, HR), 1, paste, collapse = ":")
HREKG[HREKG == "TRUE:LOW"] = sample(LV3, length(which(HREKG == "TRUE:LOW")), prob = c(0.33, 0.33, 0.33), replace = TRUE)
HREKG[HREKG == "FALSE:LOW"] = sample(LV3, length(which(HREKG == "FALSE:LOW")), prob = c(0.33, 0.33, 0.33), replace = TRUE)
HREKG[HREKG == "TRUE:NORMAL"] = sample(LV3, length(which(HREKG == "TRUE:NORMAL")), prob = c(0.33, 0.33, 0.33), replace = TRUE)
HREKG[HREKG == "FALSE:NORMAL"] = sample(LV3, length(which(HREKG == "FALSE:NORMAL")), prob = c(0.98, 0.01, 0.01), replace = TRUE)
HREKG[HREKG == "TRUE:HIGH"] = sample(LV3, length(which(HREKG == "TRUE:HIGH")), prob = c(0.01, 0.98, 0.01), replace = TRUE)
HREKG[HREKG == "FALSE:HIGH"] = sample(LV3, length(which(HREKG == "FALSE:HIGH")), prob = c(0.01, 0.01, 0.98), replace = TRUE)

HRSAT = apply(cbind(ERRCAUTER, HR), 1, paste, collapse = ":")
HRSAT[HRSAT == "TRUE:LOW"] = sample(LV3, length(which(HRSAT == "TRUE:LOW")), prob = c(0.33, 0.33, 0.33), replace = TRUE)
HRSAT[HRSAT == "FALSE:LOW"] = sample(LV3, length(which(HRSAT == "FALSE:LOW")), prob = c(0.33, 0.33, 0.33), replace = TRUE)
HRSAT[HRSAT == "TRUE:NORMAL"] = sample(LV3, length(which(HRSAT == "TRUE:NORMAL")), prob = c(0.33, 0.33, 0.33), replace = TRUE)
HRSAT[HRSAT == "FALSE:NORMAL"] = sample(LV3, length(which(HRSAT == "FALSE:NORMAL")), prob = c(0.98, 0.01, 0.01), replace = TRUE)
HRSAT[HRSAT == "TRUE:HIGH"] = sample(LV3, length(which(HRSAT == "TRUE:HIGH")), prob = c(0.01, 0.98, 0.01), replace = TRUE)
HRSAT[HRSAT == "FALSE:HIGH"] = sample(LV3, length(which(HRSAT == "FALSE:HIGH")), prob = c(0.01, 0.01, 0.98), replace = TRUE)

CO = apply(cbind(HR, STROKEVOLUME), 1, paste, collapse = ":")
CO[CO == "LOW:LOW"] = sample(LV3, length(which(CO == "LOW:LOW")), prob = c(0.98, 0.01, 0.01), replace = TRUE)
CO[CO == "NORMAL:LOW"] = sample(LV3, length(which(CO == "NORMAL:LOW")), prob = c(0.95, 0.04, 0.01), replace = TRUE)
CO[CO == "HIGH:LOW"] = sample(LV3, length(which(CO == "HIGH:LOW")), prob = c(0.80, 0.19, 0.01), replace = TRUE)
CO[CO == "LOW:NORMAL"] = sample(LV3, length(which(CO == "LOW:NORMAL")), prob = c(0.95, 0.04, 0.01), replace = TRUE)
CO[CO == "NORMAL:NORMAL"] = sample(LV3, length(which(CO == "NORMAL:NORMAL")), prob = c(0.04, 0.95, 0.01), replace = TRUE)
CO[CO == "HIGH:NORMAL"] = sample(LV3, length(which(CO == "HIGH:NORMAL")), prob = c(0.01, 0.04, 0.95), replace = TRUE)
CO[CO == "LOW:HIGH"] = sample(LV3, length(which(CO == "LOW:HIGH")), prob = c(0.30, 0.69, 0.01), replace = TRUE)
CO[CO == "NORMAL:HIGH"] = sample(LV3, length(which(CO == "NORMAL:HIGH")), prob = c(0.01, 0.30, 0.69), replace = TRUE)
CO[CO == "HIGH:HIGH"] = sample(LV3, length(which(CO == "HIGH:HIGH")), prob = c(0.01, 0.01, 0.98), replace = TRUE)

BP = apply(cbind(CO, TPR), 1, paste, collapse = ":")
BP[BP == "LOW:LOW"] = sample(LV3, length(which(BP == "LOW:LOW")), prob = c(0.98, 0.01, 0.01), replace = TRUE)
BP[BP == "NORMAL:LOW"] = sample(LV3, length(which(BP == "NORMAL:LOW")), prob = c(0.98, 0.01, 0.01), replace = TRUE)
BP[BP == "HIGH:LOW"] = sample(LV3, length(which(BP == "HIGH:LOW")), prob = c(0.90, 0.09, 0.01), replace = TRUE)
BP[BP == "LOW:NORMAL"] = sample(LV3, length(which(BP == "LOW:NORMAL")), prob = c(0.98, 0.01, 0.01), replace = TRUE)
BP[BP == "NORMAL:NORMAL"] = sample(LV3, length(which(BP == "NORMAL:NORMAL")), prob = c(0.10, 0.85, 0.05), replace = TRUE)
BP[BP == "HIGH:NORMAL"] = sample(LV3, length(which(BP == "HIGH:NORMAL")), prob = c(0.05, 0.20, 0.75), replace = TRUE)
BP[BP == "LOW:HIGH"] = sample(LV3, length(which(BP == "LOW:HIGH")), prob = c(0.30, 0.60, 0.10), replace = TRUE)
BP[BP == "NORMAL:HIGH"] = sample(LV3, length(which(BP == "NORMAL:HIGH")), prob = c(0.05, 0.40, 0.55), replace = TRUE)
BP[BP == "HIGH:HIGH"] = sample(LV3, length(which(BP == "HIGH:HIGH")), prob = c(0.01, 0.09, 0.90), replace = TRUE)

alarm = data.frame(
  CVP  = factor(CVP, levels = LV3),
  PCWP = factor(PCWP, levels = LV3),
  HIST = factor(HISTORY, levels = BOOL),
  TPR  = factor(TPR, levels = LV3),
  BP   = factor(BP, levels = LV3),
  CO   = factor(CO, levels = LV3),
  HRBP = factor(HRBP, levels = LV3),
  HREK = factor(HREKG, levels = LV3),
  HRSA = factor(HRSAT, levels = LV3),
  PAP  = factor(PAP, levels = LV3),
  SAO2 = factor(SAO2, levels = LV3),
  FIO2 = factor(FIO2, levels = c("LOW", "NORMAL")),
  PRSS = factor(PRESS, levels = LV4),
  ECO2 = factor(EXPCO2, levels = LV4),
  MINV = factor(MINVOL, levels = LV4),
  MVS  = factor(MINVOLSET, levels = LV3),
  HYP  = factor(HYPOVOLEMIA, levels = BOOL),
  LVF  = factor(LVFAILURE, levels = BOOL),
  APL  = factor(ANAPHYLAXIS, levels = BOOL),
  ANES = factor(INSUFFANESTH, levels = BOOL),
  PMB  = factor(PULMEMBOLUS, levels = BOOL),
  INT  = factor(INTUBATION, levels = c("NORMAL", "ESOPHAGEAL", "ONESIDED")),
  KINK = factor(KINKEDTUBE, levels = BOOL),
  DISC = factor(DISCONNECT, levels = BOOL),
  LVV  = factor(LVEDVOLUME, levels = LV3),
  STKV = factor(STROKEVOLUME, levels = LV3),
  CCHL = factor(CATECHOL, levels = LV2),
  ERLO = factor(ERRLOWOUTPUT, levels = BOOL),
  HR   = factor(HR, levels = LV3),
  ERCA = factor(ERRCAUTER, levels = BOOL),
  SHNT = factor(SHUNT, levels = c("NORMAL", "HIGH")),
  PVS  = factor(PVSAT, levels = LV3),
  ACO2 = factor(ARTCO2, levels = LV3),
  VALV = factor(VENTALV, levels = LV4),
  VLNG = factor(VENTLUNG, levels = LV4),
  VTUB = factor(VENTTUBE, levels = LV4),
  VMCH = factor(VENTMACH, levels = LV4)
)
