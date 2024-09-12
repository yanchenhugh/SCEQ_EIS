$onMulti

set
t time /1*273/
isim simul # /1*1500/
inode GHQ node # /1*7/
auxi 
;
alias(tt,t);
alias(inode,jnode);

scalars
alfa /0.33/ 
bet  /0.96/
delt /0.1/
gam  /2.0/
fai  /0.975/
rou_om /0.9/
STD_om /0.013/

DT      number of periods for optimization in SCEQ /200/
Tstar   number of periods of interest / 60 /
Nsim    number of simulations /1500/
s       starting period #
middle_index

;

middle_index = ceil(card(inode)/2);

Set
endonames /y, k, ctot, inve, meu/
endonames_Euler /k, ctot, meu/
;

Parameters
ghq(inode, auxi)
ghq_nodes(inode)
ghq_weights(inode)
;

** Now load gauss-hermite quadrature nodes
$call csv2gdx GaussHermite_7.csv id=ghq index=1 values=2..lastCol useHeader=y
$gdxIn GaussHermite_7.gdx
$load inode = Dim1
$load auxi = Dim2
$load ghq

ghq_nodes(inode) = ghq(inode, 'node') * sqrt(2);
ghq_weights(inode) = ghq(inode, 'weight') / sqrt(pi);

display ghq, ghq_nodes, ghq_weights;

** Now compute steadystate
Scalars
yss
kss
ctotss
invess
meuss
omss /1.0/
;

invess = delt * ( (1/bet-1+delt)/alfa )**(1/(alfa-1));
kss = invess/delt;
yss = omss * kss**alfa;
ctotss = yss - invess;
meuss = 0;

display yss, ctotss, invess;
