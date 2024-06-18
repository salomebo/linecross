## 	Script to make the sysdata.rda in the R/ folder

S=c(P1 = 1, 
    P2 = 0, 
    F1 = 1/2, F1_12 = 1/2, F1_21 = 1/2, 
    F2 = 1/2, F2_12 = 1/2, F2_11 = 1/2, F2_22 = 1/2, F2_21 = 1/2,
    F3 = 1/2,
    B1 = 3/4, B1_12 = 3/4, B1_11a= 3/4, B1_11b= 3/4, B1_21 = 3/4, 
    B2 = 1/4, B2_21 = 1/4, B2_22a= 1/4, B2_22b= 1/4, B2_12 = 1/4,
    P1xF2 = 3/4,
    P1xB1 = 7/8,
    P1xB2 = 5/8,
    P2xF2 = 1/4,
    P2xB1 = 3/8,
    P2xB2 = 1/8,
    F1xF2 = 1/2,
    F1xB1 = 5/8,
    F1xB2 = 3/8,
    F2xB1 = 5/8,
    F2xB2 = 3/8,
    B1xB1 = 3/4,
    B1xB2 = 1/2,
    B2xB2 = 1/4)

H=c(P1 = 0, 
    P2 = 0, 
    F1 = 1, F1_12 = 1, F1_21 =1, 
    F2 = 1/2, F2_12 = 1/2, F2_11 = 1/2, F2_22 = 1/2, F2_21 = 1/2,
    F3 = 1/2,
    B1 = 1/2, B1_12 = 1/2, B1_11a= 1/2, B1_11b= 1/2, B1_21 = 1/2, 
    B2 = 1/2, B2_21 = 1/2, B2_22a= 1/2, B2_22b= 1/2, B2_12 = 1/2,
    P1xF2 = 1/2,
    P1xB1 = 1/4,
    P1xB2 = 3/4,
    P2xF2 = 1/2,
    P2xB1 = 3/4,
    P2xB2 = 1/4,
    F1xF2 = 1/2,
    F1xB1 = 1/2,
    F1xB2 = 1/2,
    F2xB1 = 1/2,
    F2xB2 = 1/2,
    B1xB1 = 3/8,
    B1xB2 = 5/8,
    B2xB2 = 3/8)

propP1 = S
hybrid.index = H

Mlist=list()
Mlist$P1=list()
Mlist$P1$reference[[1]]=c("P1")
Mlist$P1$additive[[1]]= cbind(1-S)
Mlist$P1$additive[[2]]= c("Y2")
Mlist$P1$additive[[3]]= "(1-S)*par[1]"
Mlist$P1$dominance[[1]]=cbind(H)
Mlist$P1$dominance[[2]]=c("Yh")
Mlist$P1$dominance[[3]]= "H*par[1]"
Mlist$P1$add_dom[[1]]=cbind(c(H), c(1-S-H/2))
Mlist$P1$add_dom[[2]]=c("Yh", "Y2")
Mlist$P1$add_dom[[3]]= "H*par[1] + (1-S-H/2)*par[2]"

Mlist$P1$general[[1]] = cbind(c(H), c(1-S-H/2), c(1/2*H^2), c(H*(1-S-H/2)), c(1/2*(1-S-H/2)^2))
Mlist$P1$general[[2]] = c("Yh", "Y2", "Ehh", "Eh2", "E22")
Mlist$P1$general[[3]] = "H*par[1] + (1-S-H/2)*par[2] + H^2*par[3] + 2*H*(1-S-H/2)*par[4] + (1-S-H/2)^2*par[5]"
Mlist$P1$general_dom[[1]] = cbind(c(H), c(H^2), c(2*H*(1-S-H/2)), c((1-S-H/2)^2))
Mlist$P1$general_dom[[2]] = c("Yh", "Ehh", "Eh2", "E22")
Mlist$P1$general_dom[[3]] = "H*par[1]  + H^2*par[2] + 2*H*(1-S-H/2)*par[3] + (1-S-H/2)^2*par[4]"

Mlist$P1$generalWB[[1]] = cbind(c(H), c(1-S-H/2), c(1/2*(1-S-H/2)^2), c(H*(1-S)))
Mlist$P1$generalWB[[2]] = c("Yh", "Y2", "Ew", "Eb")
Mlist$P1$generalWB[[3]] = "H*par[1] + (1-S-H/2)*par[2] + (1-S-H/2)^2*par[3] + 2*H*(1-S)*par[4]"
Mlist$P1$generalWB_dom[[1]] = cbind(c(H), c((1-S-H/2)^2), c(2*H*(1-S)))
Mlist$P1$generalWB_dom[[2]] = c("Yh",  "Ew", "Eb")
Mlist$P1$generalWB_dom[[3]] = "H*par[1] +  (1-S-H/2)^2*par[2] + 2*H*(1-S)*par[3]"

Mlist$P1$generalW[[1]] = cbind(c(H), c(1-S-H/2), c(1/2*(1-S-H/2)^2))
Mlist$P1$generalW[[2]] = c("Yh", "Y2", "Ew")
Mlist$P1$generalW[[3]] = "H*par[1] + (1-S-H/2)*par[2] + (1-S-H/2)^2*par[3]"
Mlist$P1$generalW_dom[[1]] = cbind(c(H), c((1-S-H/2)^2))
Mlist$P1$generalW_dom[[2]] = c("Yh",  "Ew")
Mlist$P1$generalW_dom[[3]] = "H*par[1] +  (1-S-H/2)^2*par[2]"

Mlist$P1$generalB[[1]] = cbind(c(H), c(1-S-H/2), c(2*H*(1-S)))
Mlist$P1$generalB[[2]] = c("Yh", "Y2", "Eb")
Mlist$P1$generalB[[3]] = "H*par[1] + (1-S-H/2)*par[2] + 2*H*(1-S)*par[3]"
Mlist$P1$generalB_dom[[1]] = cbind(c(H), c(2*H*(1-S)))
Mlist$P1$generalB_dom[[2]] = c("Yh",  "Eb")
Mlist$P1$generalB_dom[[3]] = "H*par[1] + 2*H*(1-S)*par[2]"

Mlist$P1$classic[[1]] = cbind(c(2*S-2), c(2*H), c((2*S-2)^2), c((4*S*H-4*H)), c(4*H^2))
Mlist$P1$classic[[2]] = c("alpha", "delta", "Eaa", "Ead", "Edd")
Mlist$P1$classic[[3]] = "(2*S-2)*par[1] + 2*H*par[2] + ((2*S-2)^2)*par[3] + (4*S*H-4*H)*par[4] + (4*H^2)*par[5]"

Mlist$P1$multilinear[[1]] = "H*par[1] + (1-S-H/2)*par[2] + (par[3]/2)*(H*par[1] + (1-S-H/2)*par[2])^2"
Mlist$P1$multilinear[[2]] = c("Yh", "Y2", "eps")
Mlist$P1$multilinear[[3]] = "H*par[1] + (1-S-H/2)*par[2]"
Mlist$P1$multilinear[[4]] = "(1/2)*(H*par[1] + (1-S-H/2)*par[2])^2"
Mlist$P1$multilinear[[5]] = "(H*par[1] + (1-S-H/2)*par[2] + (par[3]/2)*(H*par[1] + (1-S-H/2)*par[2])^2)+par[4]"

Mlist$P1$multilinear_add[[1]]="(1-S)*par[1] + (par[2]/2)*((1-S)*par[1])^2"
Mlist$P1$multilinear_add[[2]]= c("Y2", "eps")
Mlist$P1$multilinear_add[[3]]="(1-S)*par[1]"
Mlist$P1$multilinear_add[[4]]="(1/2)*((1-S)*par[1])^2"
Mlist$P1$multilinear_add[[5]]="((1-S)*par[1] + (par[2]/2)*((1-S)*par[1])^2)+par[3]"

Mlist$P1$canalization[[1]]= "H*par[1] + (1-S-H/2)*par[2] + sign(H*par[1] + (1-S-H/2)*par[2])*(par[3]/2)*(H*par[1] + (1-S-H/2)*par[2])^2"
Mlist$P1$canalization[[2]]=c("Yh", "Y2", "eps")
Mlist$P1$canalization[[3]]="H*par[1] + (1-S-H/2)*par[2]"
Mlist$P1$canalization[[4]]="sign(H*par[1] + (1-S-H/2)*par[2])*(1/2)*(H*par[1] + (1-S-H/2)*par[2])^2"
Mlist$P1$canalization[[5]]="(H*par[1] + (1-S-H/2)*par[2] + sign(H*par[1] + (1-S-H/2)*par[2])*(par[3]/2)*(H*par[1] + (1-S-H/2)*par[2])^2)+par[4]"

Mlist$P1$dominance.estimate[[1]] = "parameters[1]/2 - parameters[2]/4"
Mlist$P1$dominance.estimate[[2]] = "parameters.se[1]"



Mlist$P2=list()
Mlist$P2$reference[[1]]=c("P2")
Mlist$P2$additive[[1]]= cbind(S)
Mlist$P2$additive[[2]]= c("Y1")
Mlist$P2$additive[[3]]= "S*par[1]"
Mlist$P2$dominance[[1]]=cbind(H)
Mlist$P2$dominance[[2]]=c("Yh")
Mlist$P2$dominance[[3]]= "H*par[1]"
Mlist$P2$add_dom[[1]]=cbind(c(H), c(S-H/2))
Mlist$P2$add_dom[[2]]=c("Yh", "Y1")
Mlist$P2$add_dom[[3]]= "H*par[1] + (S-H/2)*par[2]"

Mlist$P2$general[[1]] = cbind(c(H), c(S-H/2),c(1/2*H^2), c(H*(S-H/2)) , c(1/2*(S-H/2)^2))
Mlist$P2$general[[2]] = c("Yh", "Y1", "Ehh", "Eh1", "E11")
Mlist$P2$general[[3]] = "H*par[1] + (S-H/2)*par[2] + H^2*par[3] + 2*H*(S-H/2)*par[4] + (S-H/2)^2*par[5]"
Mlist$P2$general_dom[[1]] = cbind(c(H), c(H^2), c(2*H*(S-H/2)) , c((S-H/2)^2))
Mlist$P2$general_dom[[2]] = c("Yh", "Ehh", "Eh1", "E11")
Mlist$P2$general_dom[[3]] = "H*par[1] + H^2*par[2] + 2*H*(S-H/2)*par[3] + (S-H/2)^2*par[4]"

Mlist$P2$generalWB[[1]] = cbind(c(H), c(S-H/2), c((S-H/2)^2), c(2*H*S))
Mlist$P2$generalWB[[2]] = c("Yh", "Y1", "Ew", "Eb")
Mlist$P2$generalWB[[3]] = "H*par[1] + (S-H/2)*par[2] + (S-H/2)^2*par[3] + 2*H*S*par[4]"
Mlist$P2$generalWB_dom[[1]] = cbind(c(H), c((S-H/2)^2), c(2*H*S))
Mlist$P2$generalWB_dom[[2]] = c("Yh", "Ew", "Eb")
Mlist$P2$generalWB_dom[[3]] = "H*par[1] + (S-H/2)^2*par[2] + 2*H*S*par[3]"

Mlist$P2$generalW[[1]] = cbind(c(H), c(S-H/2), c((S-H/2)^2))
Mlist$P2$generalW[[2]] = c("Yh", "Y1", "Ew")
Mlist$P2$generalW[[3]] = "H*par[1] + (S-H/2)*par[2] + (S-H/2)^2*par[3]"
Mlist$P2$generalW_dom[[1]] = cbind(c(H), c((S-H/2)^2))
Mlist$P2$generalW_dom[[2]] = c("Yh", "Ew")
Mlist$P2$generalW_dom[[3]] = "H*par[1] + (S-H/2)^2*par[2]"

Mlist$P2$generalB[[1]] = cbind(c(H), c(S-H/2), c(2*H*S))
Mlist$P2$generalB[[2]] = c("Yh", "Y1", "Eb")
Mlist$P2$generalB[[3]] = "H*par[1] + (S-H/2)*par[2] + 2*H*S*par[3]"
Mlist$P2$generalB_dom[[1]] = cbind(c(H), c(2*H*S))
Mlist$P2$generalB_dom[[2]] = c("Yh", "Eb")
Mlist$P2$generalB_dom[[3]] = "H*par[1] + 2*H*S*par[2]"

Mlist$P2$classic[[1]] = cbind(c(2*S), c(2*H), c(4*S^2), c(4*S*H), c(4*H^2))
Mlist$P2$classic[[2]] = c("alpha", "delta", "Eaa", "Ead", "Edd")
Mlist$P2$classic[[3]] = "(2*S)*par[1] + 2*H*par[2] + 4*S^2*par[3] + 4*S*H*par[4] + 4*H^2*par[5]"

Mlist$P2$multilinear[[1]] = "H*par[1] + (S-H/2)*par[2] + (par[3]/2)*(H*par[1] + (S-H/2)*par[2])^2"
Mlist$P2$multilinear[[2]] = c("Yh", "Y1", "eps")
Mlist$P2$multilinear[[3]] = "H*par[1] + (S-H/2)*par[2]"
Mlist$P2$multilinear[[4]] = "(1/2)*(H*par[1] + (S-H/2)*par[2])^2"
Mlist$P2$multilinear[[5]] = "(H*par[1] + (S-H/2)*par[2] + (par[3]/2)*(H*par[1] + (S-H/2)*par[2])^2)+par[4]"

Mlist$P2$multilinear_add[[1]]="(S-H/2)*par[1] + (par[2]/2)*((S-H/2)*par[1])^2"
Mlist$P2$multilinear_add[[2]]= c("Y1", "eps")
Mlist$P2$multilinear_add[[3]]="(S-H/2)*par[1]"
Mlist$P2$multilinear_add[[4]]="(1/2)*((S-H/2)*par[1])^2"
Mlist$P2$multilinear_add[[5]]= "((S-H/2)*par[1] + (par[2]/2)*((S-H/2)*par[1])^2)+par[3]"

Mlist$P2$canalization[[1]]= "H*par[1] + (S-H/2)*par[2] + sign(H*par[1] + (S-H/2)*par[2])*(par[3]/2)*(H*par[1] + (S-H/2)*par[2])^2"
Mlist$P2$canalization[[2]]= c("Yh", "Y1", "eps")
Mlist$P2$canalization[[3]]="H*par[1] + (S-H/2)*par[2]"
Mlist$P2$canalization[[4]]="sign(H*par[1] + (S-H/2)*par[2])*(1/2)*(H*par[1] + (S-H/2)*par[2])^2"
Mlist$P2$canalization[[5]]= "(H*par[1] + (S-H/2)*par[2] + sign(H*par[1] + (S-H/2)*par[2])*(par[3]/2)*(H*par[1] + (S-H/2)*par[2])^2)+par[4]"

Mlist$P2$dominance.estimate[[1]] = "parameters[1]/2 - parameters[2]/4"
Mlist$P2$dominance.estimate[[2]] = "parameters.se[1]"



Mlist$F1=list()
Mlist$F1$reference[[1]]= c("F1", "F1_12", "F1_21")
Mlist$F1$additive[[1]]= cbind(2*S-1)
Mlist$F1$additive[[2]]= c("Y")
Mlist$F1$additive[[3]]= "(2*S-1)*par[1]"
Mlist$F1$dominance[[1]]= cbind(2*H-2) ## is it right ? - could it be cbind(c(S-H/2),c(1-S-H/2)) ?
Mlist$F1$dominance[[2]]= c("Y")
Mlist$F1$dominance[[3]]= "(2*H-2)*par[1]"
Mlist$F1$add_dom[[1]]= cbind(c(S-H/2), c(1-S-H/2))
Mlist$F1$add_dom[[2]]= c("Y1", "Y2")
Mlist$F1$add_dom[[3]]= "(S-H/2)*par[1] +  (1-S-H/2)*par[2]"

Mlist$F1$general[[1]] = cbind(c(S-H/2), c(1-S-H/2), c((S-H/2)^2), c(2*(S-H/2)*(1-S-H/2)), c((1-S-H/2)^2))
Mlist$F1$general[[2]] = c("Y1", "Y2", "E11", "E12", "E22")
Mlist$F1$general[[3]] = "(S-H/2)*par[1] + (1-S-H/2)*par[2] + (S-H/2)^2*par[3] + 2*(S-H/2)*(1-S-H/2)*par[4] + (1-S-H/2)^2*par[5]"
Mlist$F1$general_dom[[1]] = cbind(c(2*H-2), c((S-H/2)^2), c(2*(S-H/2)*(1-S-H/2)), c((1-S-H/2)^2))
Mlist$F1$general_dom[[2]] = c("d", "E11", "E12", "E22")
Mlist$F1$general_dom[[3]] = "(2*H-2)*par[1] +  (S-H/2)^2*par[2] + 2*(S-H/2)*(1-S-H/2)*par[3] + (1-S-H/2)^2*par[4]"

Mlist$F1$generalWB[[1]] = cbind(c(S-H/2), c(1-S-H/2), c((S-H/2)^2+(1-S-H/2)^2), c(2*(S-H/2)*(1-S-H/2)))
Mlist$F1$generalWB[[2]] = c("Y1", "Y2", "Ew", "Eb")
Mlist$F1$generalWB[[3]] = "(S-H/2)*par[1] + (1-S-H/2)*par[2] + ((S-H/2)^2 +(1-S-H/2)^2)*par[3] + 2*(S-H/2)*(1-S-H/2)*par[4]"
Mlist$F1$generalWB_dom[[1]] = cbind(c(2*H-2), c((S-H/2)^2+(1-S-H/2)^2), c(2*(S-H/2)*(1-S-H/2)))
Mlist$F1$generalWB_dom[[2]] = c("d", "Ew", "Eb")
Mlist$F1$generalWB_dom[[3]] = "(2*H-2)*par[1] + ((S-H/2)^2+(1-S-H/2)^2)*par[2] + 2*(S-H/2)*(1-S-H/2)*par[3]"

Mlist$F1$generalW[[1]] = cbind(c(S-H/2), c(1-S-H/2), c((S-H/2)^2+(1-S-H/2)^2))
Mlist$F1$generalW[[2]] = c("Y1", "Y2", "Ew")
Mlist$F1$generalW[[3]] = "(S-H/2)*par[1] + (1-S-H/2)*par[2] + ((S-H/2)^2+(1-S-H/2)^2)*par[3]"
Mlist$F1$generalW_dom[[1]] = cbind(c(2*H-2), c((S-H/2)^2+(1-S-H/2)^2))
Mlist$F1$generalW_dom[[2]] = c("d", "Ew")
Mlist$F1$generalW_dom[[3]] = "(2*H-2)*par[1] + ((S-H/2)^2+(1-S-H/2)^2)*par[2]"

Mlist$F1$generalB[[1]] = cbind(c(S-H/2), c(1-S-H/2), c(2*(S-H/2)*(1-S-H/2)))
Mlist$F1$generalB[[2]] = c("Y1", "Y2", "Eb")
Mlist$F1$generalB[[3]] = "(S-H/2)*par[1] + (1-S-H/2)*par[2] + 2*(S-H/2)*(1-S-H/2)*par[3]"
Mlist$F1$generalB_dom[[1]] = cbind(c(2*H-2), c(2*(S-H/2)*(1-S-H/2)))
Mlist$F1$generalB_dom[[2]] = c("d", "Eb")
Mlist$F1$generalB_dom[[3]] = "(2*H-2)*par[1] + 2*(S-H/2)*(1-S-H/2)*par[2]"


Mlist$F1$classic[[1]] = cbind(c(2*S-1), c(2*H-2), c((2*S-1)^2), c((2*S-1)*(2*H-2)), c((2*H-2)^2))
Mlist$F1$classic[[2]] = c("alpha", "delta", "Eaa", "Ead", "Edd")
Mlist$F1$classic[[3]] = "(2*S-1)*par[1] + (2*H-2)*par[2] + ((2*S-1)^2)*par[3] + (2*S-1)*(2*H-2)*par[4] + ((2*H-2)^2)*par[5]"

Mlist$F1$multilinear[[1]] = "(S-H/2)*par[1] + (1-S-H/2)*par[2] + (par[3]/2)*((S-H/2)*par[1] + (1-S-H/2)*par[2])^2"
Mlist$F1$multilinear[[2]] = c("Y1", "Y2", "eps")
Mlist$F1$multilinear[[3]] = "(S-H/2)*par[1] + (1-S-H/2)*par[2]"
Mlist$F1$multilinear[[4]] = "(1/2)*((S-H/2)*par[1] + (1-S-H/2)*par[2])^2"
Mlist$F1$multilinear[[5]] = "((S-H/2)*par[1] + (1-S-H/2)*par[2] + (par[3]/2)*((S-H/2)*par[1] + (1-S-H/2)*par[2])^2)+par[4]"

Mlist$F1$multilinear_add[[1]]= "(S-H/2)*par[1]+ (1-S)*par[2] + (par[3]/2)*((S-H/2)*par[1] + (1-S)*par[2])^2"
Mlist$F1$multilinear_add[[2]]=c("Y", "eps")
Mlist$F1$multilinear_add[[3]]= "(S-H/2)*par[1]+ (1-S)*par[2]"
Mlist$F1$multilinear_add[[4]]= "(1/2)*((S-H/2)*par[1] + (1-S)*par[2])^2"
Mlist$F1$multilinear_add[[5]]= "((S-H/2)*par[1]+ (1-S)*par[2] + (par[3]/2)*((S-H/2)*par[1] + (1-S)*par[2])^2)+par[4]"

Mlist$F1$canalization[[1]]= "(S-H/2)*par[1] + (1-S-H/2)*par[2] + sign((S-H/2)*par[1] + (1-S-H/2)*par[2])*(par[3]/2)*((S-H/2)*par[1] + (1-S-H/2)*par[2])^2"
Mlist$F1$canalization[[2]]= c("Y1", "Y2", "eps")
Mlist$F1$canalization[[3]]= "(S-H/2)*par[1] + (1-S-H/2)*par[2]"
Mlist$F1$canalization[[4]]= "sign((S-H/2)*par[1] + (1-S-H/2)*par[2])*(1/2)*((S-H/2)*par[1] + (1-S-H/2)*par[2])^2"
Mlist$F1$canalization[[5]]= "((S-H/2)*par[1] + (1-S-H/2)*par[2] + sign((S-H/2)*par[1] + (1-S-H/2)*par[2])*(par[3]/2)*((S-H/2)*par[1] + (1-S-H/2)*par[2])^2)+par[4]"

Mlist$F1$dominance.estimate[[1]] = "(-sum(parameters[1],parameters[2]))/4"
Mlist$F1$dominance.estimate[[2]] = "NA"


Mlist$F2=list()
Mlist$F2$reference[[1]] = c("F2" ,"F2_12", "F2_11", "F2_22", "F2_21")
Mlist$F2$additive[[1]]= cbind(2*S-1)
Mlist$F2$additive[[2]]= c("Y")
Mlist$F2$additive[[3]]= "(2*S-1)*par[1]"
Mlist$F2$dominance[[1]]= cbind(2*H-1)
Mlist$F2$dominance[[2]]= c("Y")
Mlist$F2$dominance[[3]]= "(2*H-1)*par[1]"
Mlist$F2$add_dom[[1]]= cbind(c(S-H), c(1-S-H))
Mlist$F2$add_dom[[2]]= c("Y1", "Y2")
Mlist$F2$add_dom[[3]]= "(S-H)*par[1] + (1-S-H)*par[2]"

Mlist$F2$general[[1]] = cbind(c(S-H),c(1-S-H),c((S-H)^2), c(2*(S-H)*(1-S-H)), c((1-S-H)^2))
Mlist$F2$general[[2]] = c("Y1", "Y2", "E11", "E12", "E22")
Mlist$F2$general[[3]] =  "(S-H)*par[1] + (1-S-H)*par[2] + (S-H)^2*par[3] + 2*(S-H)*(1-S-H)*par[4] + (1-S-H)^2*par[5]"
Mlist$F2$general_dom[[1]] = cbind(c(2*H-1),c((S-H)^2), c(2*(S-H)*(1-S-H)), c((1-S-H)^2))
Mlist$F2$general_dom[[2]] = c("Yh", "E11", "E12", "E22")
Mlist$F2$general_dom[[3]] =  "(2*H-1)*par[1] + (S-H)^2*par[2] + 2*(S-H)*(1-S-H)*par[3] + (1-S-H)^2*par[4]"

Mlist$F2$generalWB[[1]] = cbind(c(S-H),c(1-S-H),c((S-H)^2 + (1-S-H)^2), c(2*(S-H)*(1-S-H)))
Mlist$F2$generalWB[[2]] = c("Y1", "Y2", "Ew", "Eb")
Mlist$F2$generalWB[[3]] = "(S-H)*par[1] + (1-S-H)*par[2] + par[3]*((S-H)^2 + (1-S-H)^2) + 2*(S-H)*(1-S-H)*par[4]"
Mlist$F2$generalWB_dom[[1]] = cbind(c(2*H-1),c((S-H)^2 + (1-S-H)^2), c(2*(S-H)*(1-S-H)))
Mlist$F2$generalWB_dom[[2]] = c("Yh", "Ew", "Eb")
Mlist$F2$generalWB_dom[[3]] = "(2*H-1)*par[1] + ((S-H)^2 + (1-S-H)^2)*par[2] + 2*(S-H)*(1-S-H)*par[3]"

Mlist$F2$generalW[[1]] = cbind(c(S-H),c(1-S-H),c((S-H)^2 + (1-S-H)^2))
Mlist$F2$generalW[[2]] = c("Y1", "Y2", "Ew")
Mlist$F2$generalW[[3]] =  "(S-H)*par[1] + (1-S-H)*par[2] + ((S-H)^2 + (1-S-H)^2)*par[3]"
Mlist$F2$generalW_dom[[1]] = cbind(c(2*H-1),c((S-H)^2 + (1-S-H)^2))
Mlist$F2$generalW_dom[[2]] = c("Yh", "Ew")
Mlist$F2$generalW_dom[[3]] =  "(2*H-1)*par[1] + ((S-H)^2 + (1-S-H)^2)*par[2]"

Mlist$F2$generalB[[1]] = cbind(c(S-H),c(1-S-H), c(2*(S-H)*(1-S-H)))
Mlist$F2$generalB[[2]] = c("Y1", "Y2", "Eb")
Mlist$F2$generalB[[3]] = "(S-H)*par[1] + (1-S-H)*par[2] + 2*(S-H)*(1-S-H)*par[3]"
Mlist$F2$generalB_dom[[1]] = cbind(c(2*H-1), c(2*(S-H)*(1-S-H)))
Mlist$F2$generalB_dom[[2]] = c("Yh", "Eb")
Mlist$F2$generalB_dom[[3]] = "(2*H-1)*par[1] + 2*(S-H)*(1-S-H)*par[2]"

Mlist$F2$classic[[1]] = cbind(c(2*S-1), c(2*H-1), c((2*S-1)^2), c((2*S-1)*(2*H-1)), c((2*H-1)^2))
Mlist$F2$classic[[2]] = c("alpha", "delta", "Eaa", "Ead", "Edd")
Mlist$F2$classic[[3]] = "(2*S-1)*par[1] + (2*H-1)*par[2] + ((2*S-1)^2)*par[3] + (2*S-1)*(2*H-1)*par[4] + ((2*H-1)^2)*par[5]"

Mlist$F2$multilinear[[1]] = "(S-H)*par[1] + (1-S-H)*par[2] + (par[3]/2)*((S-H)*par[1] + (1-S-H)*par[2])^2"
Mlist$F2$multilinear[[2]] = c("Y1", "Y2", "eps")
Mlist$F2$multilinear[[3]] = "(S-H)*par[1] + (1-S-H)*par[2]"
Mlist$F2$multilinear[[4]] = "(1/2)*((S-H)*par[1] + (1-S-H)*par[2])^2"
Mlist$F2$multilinear[[5]] = "((S-H)*par[1] + (1-S-H)*par[2] + (par[3]/2)*((S-H)*par[1] + (1-S-H)*par[2])^2)+par[4]"

Mlist$F2$multilinear_add[[1]]= "(S-H/2)*par[1]+ (1-S)*par[2] + (par[3]/2)*((S-H/2)*par[1] + (1-S)*par[2])^2"
Mlist$F2$multilinear_add[[2]]= c("Y", "eps")
Mlist$F2$multilinear_add[[3]]= "(S-H/2)*par[1]+ (1-S)*par[2]"
Mlist$F2$multilinear_add[[4]]= "(1/2)*((S-H/2)*par[1] + (1-S)*par[2])^2"
Mlist$F2$multilinear_add[[5]]= "((S-H/2)*par[1]+ (1-S)*par[2] + (par[3]/2)*((S-H/2)*par[1] + (1-S)*par[2])^2)+par[4]"

Mlist$F2$canalization[[1]]= "(S-H)*par[1] + (1-S-H)*par[2] + sign((S-H)*par[1] + (1-S-H)*par[2])*(par[3]/2)*((S-H)*par[1] + (1-S-H)*par[2])^2"
Mlist$F2$canalization[[2]]= c("Y1", "Y2", "eps")
Mlist$F2$canalization[[3]]= "(S-H)*par[1] + (1-S-H)*par[2]"
Mlist$F2$canalization[[4]]= "sign((S-H)*par[1] + (1-S-H)*par[2])*(1/2)*((S-H)*par[1] + (1-S-H)*par[2])^2"
Mlist$F2$canalization[[5]]= "((S-H)*par[1] + (1-S-H)*par[2] + sign((S-H)*par[1] + (1-S-H)*par[2])*(par[3]/2)*((S-H)*par[1] + (1-S-H)*par[2])^2)+par[4]"

Mlist$F2$dominance.estimate[[1]] = "(-sum(parameters[1],parameters[2]))/2"
Mlist$F2$dominance.estimate[[2]] = "sqrt(sum(parameters.se[1]^2,parameters.se[2]^2))"

# saving these elements to sysdata file in the R-folder
save(Mlist, S, H, propP1, hybrid.index, file = "R/sysdata.rda")

