scalar starttime;
starttime = jnow;

option decimals = 6;

$include IRR_params_simp.gms
*$include IRR_params.gms

*******************************************************
* define Model

parameter
om(t) aggregate TFP shock
*ompath(isim,t) path for tfp exog component
omshocks(isim,t) random normalshocks for TFP
;

execseed = 1 + gmillisec(jnow);
omshocks(isim,t) = normal(0,1);

*ompath(isim,'1') = omss;
*loop(isim$(ord(isim)<=card(isim)),
*  loop(t$(ord(t) <= Tstar),
*    ompath(isim,t+1) = ompath(isim,t)**rou_om*exp(STD_om*omshocks(isim,t));
*  );
*);


Variables
obj, y(t), k(t), ctot(t), inve(t), meu(t)
;

Equations
Objective
Objective_FSB
Objective_abs

ProdTech(t)
EulerK(t)
EulerK_EIS(t)
LOMK(t)
Resource(t)
InvConstraint(t)
InvConstraint_OBC(t)
InvConstraint_OBC_EIS(t)
InvConstraint_INEQ(t)
InvConstraint_KKT(t)

TermC(t)
;

Objective.. 
obj =e= sum( t$(ord(t)>=DT+s-5 and ord(t)<=DT+s),
power(y(t)-yss,2) + power(inve(t)-invess,2) + power(meu(t) - meuss,2) + 
power(k(t)-kss,2) + power(ctot(t)-ctotss,2) );

Objective_FSB..
obj =e= 1;

Objective_abs..
obj =e= sum( t$(ord(t)>=DT+s-5 and ord(t)<=DT+s),
abs(y(t)-yss) + abs(inve(t)-invess) + abs(meu(t) - meuss) + 
abs(k(t)-kss) + abs(ctot(t)-ctotss) );

ProdTech(t)$(ord(t)>=s and ord(t)<DT+s)..
y(t) =e= om(t)*k(t-1)**alfa;

EulerK(t)$(ord(t)>=s and ord(t)<DT+s)..
ctot(t)**(-gam) - (1-bet*(1-delt))*meu(t) =e= bet*ctot(t+1)**(-gam)*( 1-delt + alfa*om(t+1)*k(t)**(alfa-1) );

EulerK_EIS(t)$(ord(t)>s and ord(t)<DT+s)..
ctot(t)**(-gam) - (1-bet*(1-delt))*meu(t) =e= bet*ctot(t+1)**(-gam)*( 1-delt + alfa*om(t+1)*k(t)**(alfa-1) );

LOMK(t)$(ord(t)>=s and ord(t)<DT+s)..
k(t) =e= (1-delt)*k(t-1) + inve(t);

Resource(t)$(ord(t)>=s and ord(t)<DT+s)..
y(t) =e= ctot(t) + inve(t);

InvConstraint(t)$(ord(t)>=s and ord(t)<DT+s)..
meu(t) =e= 1.e-12;

InvConstraint_OBC(t)$(ord(t)>=s and ord(t)<DT+s)..
0 =e= min(meu(t), inve(t)-fai*invess);

InvConstraint_OBC_EIS(t)$(ord(t)>s and ord(t)<DT+s)..
0 =e= min(meu(t), inve(t)-fai*invess);

InvConstraint_INEQ(t)$(ord(t)>=s and ord(t)<DT+s)..
inve(t) =g= fai*invess;

InvConstraint_KKT(t)$(ord(t)>=s and ord(t)<DT+s)..
meu(t) * ( inve(t) - fai*invess ) =e= 0;

TermC(t)$(ord(t)=DT+s)..
ctot(t) =e= ctotss;

Option limrow=0,limcol=0,solprint=off;
OPTION ITERLIM = 1000;
OPTION RESLIM = 500000;
OPTION DOMLIM = 1000;
OPTION DNLP = conopt4;
OPTION NLP = conopt4;
OPTION MCP = PATH;

model
IRR_lsq_om /Objective, ProdTech, EulerK_EIS, LOMK, Resource, InvConstraint_OBC/
IRR_cns_om /Objective_FSB, ProdTech, EulerK_EIS, LOMK, Resource, InvConstraint_OBC, TermC/
IRR_fsb_om /Objective_FSB, ProdTech, EulerK_EIS, LOMK, Resource, InvConstraint_OBC, TermC/
IRR_lsq_can_om /Objective, ProdTech, EulerK, LOMK, Resource, InvConstraint_OBC/
IRR_cns_can_om /Objective_FSB, ProdTech, EulerK, LOMK, Resource, InvConstraint_OBC, TermC/
IRR_fsb_can_om /Objective_FSB, ProdTech, EulerK, LOMK, Resource, InvConstraint_OBC, TermC/
IRR_lsq_df_om /Objective, ProdTech, EulerK_EIS, LOMK, Resource, InvConstraint/
IRR_cns_df_om /Objective_FSB, ProdTech, EulerK_EIS, LOMK, Resource, InvConstraint, TermC/
IRR_fsb_df_om /Objective_FSB, ProdTech, EulerK_EIS, LOMK, Resource, InvConstraint_KKT, TermC/
IRR_lsq_can_df_om /Objective, ProdTech, EulerK, LOMK, Resource, InvConstraint/
IRR_cns_can_df_om /Objective_FSB, ProdTech, EulerK, LOMK, Resource, InvConstraint, TermC/
IRR_fsb_can_df_om /Objective_FSB, ProdTech, EulerK, LOMK, Resource, InvConstraint_KKT, TermC/
;

$onecho > conopt4.opt
Num_Rounds 1
$offecho

Parameter
nnnmax / 5 /
nsuccmax / 1 /
nnn
nsucc
;

*** write model solution to csv file for matlab to read, initiation of files
Parameter
modsols(t,endonames)
ksols_ghq(t,inode)
ysols_ghq(t,inode)
invesols_ghq(t,inode)
ctotsols_ghq(t,inode)
meusols_ghq(t,inode)
om_ghq(t,inode)
moderr(inode)
moderrs(t)
regime1(t)
success(t)
success_errchk(t)
k_temp(t)
y_temp(t)
ctot_temp(t)
inve_temp(t)
meu_temp(t)
;

File model_sols /model_sols.csv/;
model_sols.pc = 5;
model_sols.ap = 0;
model_sols.pw = 30000;

put model_sols;

put 'model solutions'
put /;


*Parameters
*k_nd_1(t)
*k_nd_2(t)
*dk_nd_12(t)
*;

******* Now start simulation *******

*********** Now start formal loops *******************
******************************************************
******************************************************
loop( isim$(ord(isim) >= 1 and ord(isim) <= card(isim)),
*loop( isim$(ord(isim) >= 1 and ord(isim) <= 1),

** Give lower bound necessary to embed the occasionally binding constraints
meu.lo(t) = 0;
k.lo(t) = 1.e-2;
y.lo(t) = 1.e-2;
ctot.lo(t) = 1.e-2;
inve.lo(t) = 1.e-2;
obj.lo = -1;
k.up(t) = 10;
y.up(t) = 10;
inve.up(t) = 10;
ctot.up(t) = 10;
meu.up(t) = 10;
obj.up = 2;

* Give initial guess
*y.l(t) = inits(t,'y');
*k.l(t) = inits(t,'k');
*ctot.l(t) = inits(t,'ctot');
*inve.l(t) = inits(t,'inve');
*meu.l(t) = inits(t,'meu');
y.l(t) = yss;
k.l(t) = kss;
ctot.l(t) = ctotss;
inve.l(t) = invess;
meu.l(t) = meuss;

* Reset shocks
om(t) = omss;

* Rest success results
success(t) = 1;
success_errchk(t) = 1;

loop( tt$(ord(tt) >= 1  and ord(tt) <= Tstar),

* current period index
s = ord(tt) + 1;

* Fix s-1 variable values
k.fx(tt) = k.l(tt);

* start with the case the t+1 shock is zero
om(tt+1) = om(tt)**rou_om*exp(omshocks(isim,tt)*STD_om);
loop( t$(ord(t)>s and ord(t)<=card(t)), om(t) = om(t-1)**rou_om );
    
nsucc = 0;
nnn = 0;
repeat(
    nnn = nnn+1;
    solve IRR_lsq_can_om using dnlp minimizing obj;
    if((IRR_lsq_can_om.MODELSTAT le 2 and IRR_lsq_can_om.SOLVESTAT eq 1),
      nsucc = nsucc+1;
    else
      nsucc = 0;
);
until (nnn>=nnnmax or nsucc>=nsuccmax) );

if (nsucc<nsuccmax,
    success(tt+1) = 0;
*    continue;
);

* Now compute GHQ error
k_temp(t) = k.l(t);
y_temp(t) = y.l(t);
ctot_temp(t) = ctot.l(t);
inve_temp(t) = inve.l(t);
meu_temp(t) = meu.l(t);

ksols_ghq(t,inode)$(ord(inode) eq middle_index) = k.l(t);
ysols_ghq(t,inode)$(ord(inode) eq middle_index) = y.l(t);
ctotsols_ghq(t,inode)$(ord(inode) eq middle_index) = ctot.l(t);
invesols_ghq(t,inode)$(ord(inode) eq middle_index) = inve.l(t);
meusols_ghq(t,inode)$(ord(inode) eq middle_index) = meu.l(t);
om_ghq(t,inode)$(ord(inode) eq middle_index) = om(t);
    
k.fx(tt+1) = k.l(tt+1);
s = ord(tt) + 2;

loop(jnode$(ord(jnode) <= card(jnode)),
    om(tt+2) = om(tt+1)**rou_om*exp(ghq_nodes(jnode)*STD_om);
    loop( t$(ord(t)>s and ord(t)<=card(t)), om(t) = om(t-1)**rou_om );
    
    nsucc = 0;
    nnn = 0;
    repeat(
        nnn = nnn+1;
        solve IRR_lsq_can_om using dnlp minimizing obj;
        if((IRR_lsq_can_om.MODELSTAT le 2 and IRR_lsq_can_om.SOLVESTAT eq 1),
          nsucc = nsucc+1;
        else
          nsucc = 0;
    );
    until (nnn>=nnnmax or nsucc>=nsuccmax) );
    
    if (nsucc<nsuccmax,
        success_errchk(tt+1) = 0;
*        break;
    );
 
    ksols_ghq(t,jnode) = k.l(t);
    ysols_ghq(t,jnode) = y.l(t);
    ctotsols_ghq(t,jnode) = ctot.l(t);
    invesols_ghq(t,jnode) = inve.l(t);
    meusols_ghq(t,jnode) = meu.l(t);
    om_ghq(t,jnode) = om(t);
    
    if (ord(jnode) eq middle_index-1,
        k.l(t) = ksols_ghq(t,jnode+1);
        y.l(t) = ysols_ghq(t,jnode+1);
        ctot.l(t) = ctotsols_ghq(t,jnode+1);
        inve.l(t) = invesols_ghq(t,jnode+1);
        meu.l(t) = meusols_ghq(t,jnode+1);
    );
);

s = ord(tt) + 1;

k.lo(tt+1) = k.lo(tt+2);
k.up(tt+1) = k.up(tt+2);

k.l(t) = k_temp(t);
y.l(t) = y_temp(t);
ctot.l(t) = ctot_temp(t);
inve.l(t) = inve_temp(t);
meu.l(t) = meu_temp(t);

** Now compute error
*ctot(t)**(-gam) - (1-bet*(1-delt))*meu(t) =e= bet*ctot(t+1)**(-gam)*( 1-delt + alfa*om(t+1)*k(t)**(alfa-1) );
moderr(inode) = ( bet*ctotsols_ghq(tt+2,inode)**(-gam)*
    ( 1-delt + alfa*om_ghq(tt+2,inode)*ksols_ghq(tt+1,inode)**(alfa-1) )
    + (1-bet*(1-delt))*meusols_ghq(tt+1,inode) ) / ctotsols_ghq(tt+1,inode)**(-gam) - 1;
display moderr;
moderrs(tt+1) = abs( sum( inode, ghq_weights(inode) * moderr(inode) ) ) ;

);

* write solution to modsols
modsols(t,'k')$(ord(t)<=DT+s) = k.l(t) - kss;
modsols(t,'y')$(ord(t)<=DT+s) = y.l(t) - yss;
modsols(t,'inve')$(ord(t)<=DT+s) = inve.l(t) - invess;
modsols(t,'ctot')$(ord(t)<=DT+s) = ctot.l(t) - ctotss;
modsols(t,'meu')$(ord(t)<=DT+s) = meu.l(t) - meuss;

regime1(t)$(ord(t)<=DT+s) = inve.l(t) < 1.e-6 + fai*invess;

* write model solution to csv file for matlab to read
model_sols.ap = 1;
put model_sols;
put '';
loop(endonames, put endonames.tl;);
put 'om';
put 'rgm1';
put 'success'
put 'succ err'
put 'error'
put /;
loop(t$(ord(t) <= s),
    put t.tl::4;
    loop(endonames,
        put modsols(t,endonames)::15; 
    );
    put (om(t)-omss)::15;
    put regime1(t);
    put success(t);
    put success_errchk(t);
    put moderrs(t)::15;
    put /;
);

*succ_sim = succ_sim + 1;
*break$(succ_sim >= Nsim);

);


