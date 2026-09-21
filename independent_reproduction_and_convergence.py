import numpy as np, pandas as pd
from scipy.linalg import expm
from scipy.special import erfinv

# Independent reproduction of the final transport-aware MATLAB model.
analysis_time = 72.0
dose = 1.0

def first_primes(n):
    out=[]; cand=2
    while len(out)<n:
        ok=True
        for p in out:
            if p*p>cand: break
            if cand%p==0:
                ok=False; break
        if ok: out.append(cand)
        cand += 1
    return out

def vdc(n, base, start):
    vals=np.zeros(n)
    for row in range(n):
        integer=row+start; frac=1.0; rad=0.0
        while integer>0:
            frac/=base
            rad += frac*(integer % base)
            integer //= base
        vals[row]=rad
    return vals

def halton(n,d,start,bases=None):
    if bases is None:
        bases=first_primes(d)
    return np.column_stack([vdc(n,bases[j],start) for j in range(d)])

def ph_activation(pH,pH50):
    return 1/(1+10**(2.5*(pH-pH50)))

def ros_activation(ros,K):
    return ros**1.3/(K**1.3+ros**1.3+np.finfo(float).eps)

def relrate(pH,ros,d,mech):
    kl,kp,p50,kr,K,ki=d
    fp=ph_activation(pH,p50)
    fr=ros_activation(ros,K)
    if mech=='baseline': return kl
    if mech=='pH': return kl+kp*fp
    if mech=='ROS': return kl+kr*fr
    if mech=='dual_additive': return kl+kp*fp+kr*fr
    if mech=='dual_interaction': return kl+kp*fp+kr*fr+ki*fp*fr
    raise ValueError(mech)

transport=dict(
    k_mps=.080,k_np_loss=.006,k_h_in=.006,k_h_out=.010,k_h_loss=.005,
    k_p_in=.010,k_p_out=.004,k_p_loss=.004,k_capcore=.012,k_core_back=.002,
    k_core_loss=.003,k_d_h_in=.008,k_d_h_out=.006,k_d_p_in=.006,
    k_d_p_out=.003,k_d_capcore=.015,k_d_corecap=.004,k_elim_b=.120,
    k_elim_h=.050,k_elim_p=.035,k_elim_core=.025,k_elim_mps=.060
)

scenario=dict(
    cap_pH=6.90,cap_ROS_mM=.050,core_pH=6.70,core_ROS_mM=.100,
    k_plaque_in=.010,k_cap_core=.012
)

bounds=np.array([
    [.001,.012],[.020,.120],[6.550,6.950],[.020,.150],[.200,2.000],[0,.060]
],float)

mechs=['baseline','pH','ROS','dual_additive','dual_interaction']
labels=['Nonresponsive','pH only','ROS only','Dual additive','Dual + interaction']

def A_matrix(d,mech,sc,tr):
    tr=tr.copy()
    tr['k_p_in']=sc['k_plaque_in']
    tr['k_capcore']=sc['k_cap_core']
    kb=relrate(7.4,.005,d,mech)
    km=relrate(7.2,.020,d,mech)
    kh=relrate(7.4,.005,d,mech)
    kc=relrate(sc['cap_pH'],sc['cap_ROS_mM'],d,mech)
    kr=relrate(sc['core_pH'],sc['core_ROS_mM'],d,mech)

    A=np.zeros((15,15))
    Nb,Nm,Nh,Nc,Nco,Db,Dm,Dh,Dc,Dco,Ab,Am,Ah,Ap,Cl=range(15)

    A[Nb,Nb]=-(tr['k_mps']+tr['k_h_in']+tr['k_p_in']+tr['k_np_loss']+kb)
    A[Nb,Nh]=tr['k_h_out']; A[Nb,Nc]=tr['k_p_out']

    A[Nm,Nb]=tr['k_mps']
    A[Nm,Nm]=-(tr['k_np_loss']+km)

    A[Nh,Nb]=tr['k_h_in']
    A[Nh,Nh]=-(tr['k_h_out']+tr['k_h_loss']+kh)

    A[Nc,Nb]=tr['k_p_in']; A[Nc,Nco]=tr['k_core_back']
    A[Nc,Nc]=-(tr['k_p_out']+tr['k_capcore']+tr['k_p_loss']+kc)

    A[Nco,Nc]=tr['k_capcore']
    A[Nco,Nco]=-(tr['k_core_back']+tr['k_core_loss']+kr)

    A[Db,Nb]=kb
    A[Db,Db]=-(tr['k_elim_b']+tr['k_d_h_in']+tr['k_d_p_in'])
    A[Db,Dh]=tr['k_d_h_out']; A[Db,Dc]=tr['k_d_p_out']

    A[Dm,Nm]=km
    A[Dm,Dm]=-tr['k_elim_mps']

    A[Dh,Nh]=kh; A[Dh,Db]=tr['k_d_h_in']
    A[Dh,Dh]=-(tr['k_elim_h']+tr['k_d_h_out'])

    A[Dc,Nc]=kc; A[Dc,Db]=tr['k_d_p_in']; A[Dc,Dco]=tr['k_d_corecap']
    A[Dc,Dc]=-(tr['k_elim_p']+tr['k_d_p_out']+tr['k_d_capcore'])

    A[Dco,Nco]=kr; A[Dco,Dc]=tr['k_d_capcore']
    A[Dco,Dco]=-(tr['k_elim_core']+tr['k_d_corecap'])

    A[Ab,Db]=1; A[Am,Dm]=1; A[Ah,Dh]=1
    A[Ap,Dc]=1; A[Ap,Dco]=1

    for j,r in [
        (Nb,tr['k_np_loss']),(Nm,tr['k_np_loss']),(Nh,tr['k_h_loss']),
        (Nc,tr['k_p_loss']),(Nco,tr['k_core_loss']),(Db,tr['k_elim_b']),
        (Dm,tr['k_elim_mps']),(Dh,tr['k_elim_h']),(Dc,tr['k_elim_p']),
        (Dco,tr['k_elim_core'])
    ]:
        A[Cl,j]=r
    return A

def sim(d,mech,sc=scenario,tr=transport):
    y=np.zeros(15); y[0]=dose
    f=expm(A_matrix(d,mech,sc,tr)*analysis_time).dot(y)
    blood,mps,healthy,plaque=f[10],f[11],f[12],f[13]
    off=blood+mps+healthy
    esr=plaque/max(healthy,np.finfo(float).eps)
    eff=plaque/max(plaque+off,np.finfo(float).eps)
    return np.array([plaque,off,esr,eff])

def pareto_select(plaque,off):
    n=len(plaque); isp=np.ones(n,dtype=bool)
    for i in range(n):
        dom=(plaque>=plaque[i])&(off<=off[i])&((plaque>plaque[i])|(off<off[i]))
        if dom.any():
            isp[i]=False
    idx=np.where(isp)[0]
    p=plaque[idx]; o=off[idx]
    pn=(p-p.min())/max(p.max()-p.min(),np.finfo(float).eps)
    on=(o.max()-o)/max(o.max()-o.min(),np.finfo(float).eps)
    return idx[np.argmin(np.sqrt((1-pn)**2+(1-on)**2))]

def screen(n):
    U=halton(n,6,20,first_primes(6))
    C=bounds[:,0]+U*(bounds[:,1]-bounds[:,0])
    results={}; sels={}
    for mech,label in zip(mechs,labels):
        plaque=np.empty(n); off=np.empty(n); esr=np.empty(n); ser=np.empty(n)
        for i,d in enumerate(C):
            r=sim(d,mech)
            b=sim(d,'baseline')
            plaque[i],off[i],esr[i]=r[:3]
            ser[i]=r[2]/max(b[2],np.finfo(float).eps)
        si=pareto_select(plaque,off)
        sels[mech]=C[si].copy()
        eff=plaque[si]/(plaque[si]+off[si])
        results[mech]=dict(
            Mechanism=label,N=n,selected_index=int(si+1),
            k_leak=C[si,0],k_pH=C[si,1],pH50=C[si,2],
            k_ROS=C[si,3],K_ROS=C[si,4],k_interaction=C[si,5],
            Plaque_AUC=plaque[si],OffTarget_AUC=off[si],
            ESR=esr[si],SER=ser[si],Efficiency=eff
        )
    return results,sels

rows=[]; sel_by_n={}
for n in [1000,2000,3000]:
    res,sels=screen(n)
    sel_by_n[n]=sels
    for m in mechs:
        rows.append(res[m])
pd.DataFrame(rows).to_csv('convergence_candidate_screening.csv',index=False)

selected_dual=sel_by_n[3000]['dual_interaction']
selected_ros=sel_by_n[3000]['ROS']

def lnm(z,cv):
    s=np.sqrt(np.log(1+cv**2))
    return np.exp(-.5*s*s+s*z)

def clip(v,lo,hi):
    return min(max(v,lo),hi)

def uncertainty(N):
    U=halton(N,19,100,first_primes(19))
    z=np.sqrt(2)*erfinv(np.clip(2*U-1,-.999999999999,.999999999999))
    vals=np.zeros((N,10))

    for i in range(N):
        sc=scenario.copy(); tr=transport.copy()
        sc['cap_pH']=clip(sc['cap_pH']+.10*z[i,6],6.5,7.3)
        sc['core_pH']=clip(sc['core_pH']+.12*z[i,7],6.1,7.2)
        sc['cap_ROS_mM']*=lnm(z[i,8],.35)
        sc['core_ROS_mM']*=lnm(z[i,9],.35)
        sc['k_plaque_in']*=lnm(z[i,10],.30)
        sc['k_cap_core']*=lnm(z[i,11],.40)

        tr['k_mps']*=lnm(z[i,12],.25)
        tr['k_h_in']*=lnm(z[i,13],.25)
        tr['k_d_h_in']*=lnm(z[i,14],.25)
        tr['k_p_out']*=lnm(z[i,15],.25)
        tr['k_elim_b']*=lnm(z[i,16],.20)
        tr['k_elim_p']*=lnm(z[i,17],.20)
        tr['k_elim_mps']*=lnm(z[i,18],.20)

        du=selected_dual.copy()
        ro=selected_ros.copy()
        for j,pidx in enumerate([0,1,3,4,5]):
            mult=lnm(z[i,j],.15)
            du[pidx]*=mult
            ro[pidx]*=mult

        du[2]=clip(du[2]+.08*z[i,5],6.45,7.05)
        ro[2]=clip(ro[2]+.08*z[i,5],6.45,7.05)
        du=np.clip(du,bounds[:,0],bounds[:,1])
        ro=np.clip(ro,bounds[:,0],bounds[:,1])

        dr=sim(du,'dual_interaction',sc,tr)
        db=sim(du,'baseline',sc,tr)
        ds=dr[2]/max(db[2],np.finfo(float).eps)

        rr=sim(ro,'ROS',sc,tr)
        rb=sim(ro,'baseline',sc,tr)
        rs=rr[2]/max(rb[2],np.finfo(float).eps)

        vals[i]=[
            dr[0],dr[1],dr[2],ds,
            rr[0],rr[1],rs,
            dr[0]-rr[0],dr[1]-rr[1],ds-rs
        ]
    return vals

rows=[]; fullmc=None
for N in [500,1000,2500,5000]:
    v=uncertainty(N)
    if N==5000:
        fullmc=v
    row={'N':N}
    names=[
        'Dual_Plaque_AUC','Dual_OffTarget_AUC','Dual_ESR','Dual_SER',
        'ROS_Plaque_AUC','ROS_OffTarget_AUC','ROS_SER',
        'Delta_Plaque_AUC','Delta_OffTarget_AUC','Delta_SER'
    ]
    for j,nm in enumerate(names):
        row[nm+'_median']=np.percentile(v[:,j],50)
        row[nm+'_p2_5']=np.percentile(v[:,j],2.5)
        row[nm+'_p97_5']=np.percentile(v[:,j],97.5)
    row['P_DualPlaque_gt_ROS']=np.mean(v[:,7]>0)
    row['P_DualOff_gt_ROS']=np.mean(v[:,8]>0)
    row['P_DualSER_gt_ROS']=np.mean(v[:,9]>0)
    rows.append(row)

pd.DataFrame(rows).to_csv('convergence_uncertainty.csv',index=False)

v=fullmc
rel=np.column_stack([
    100*v[:,7]/np.maximum(v[:,4],np.finfo(float).eps),
    100*v[:,8]/np.maximum(v[:,5],np.finfo(float).eps),
    100*v[:,9]/np.maximum(v[:,6],np.finfo(float).eps)
])

eff=[]
for j,nm in enumerate([
    'Plaque_AUC_change_pct','OffTarget_AUC_change_pct','SER_change_pct'
]):
    eff.append([
        nm,*np.percentile(rel[:,j],[2.5,50,97.5]),np.mean(rel[:,j]>0)
    ])

pd.DataFrame(
    eff,
    columns=['Metric','P2_5','Median','P97_5','Probability_gt_0']
).to_csv('paired_effect_sizes_python_check.csv',index=False)

sumrows=[]
for idx,nm in [(0,'Plaque_AUC'),(1,'OffTarget_AUC'),(3,'SER')]:
    rosidx={0:4,1:5,3:6}[idx]
    dual=v[:,idx]; ros=v[:,rosidx]; diff=dual-ros
    sumrows.append([
        nm,
        np.percentile(dual,50),np.percentile(dual,2.5),np.percentile(dual,97.5),
        np.percentile(ros,50),np.percentile(ros,2.5),np.percentile(ros,97.5),
        np.percentile(diff,50),np.percentile(diff,2.5),np.percentile(diff,97.5),
        np.mean(diff>0)
    ])

pd.DataFrame(
    sumrows,
    columns=[
        'Metric','Dual_median','Dual_p2_5','Dual_p97_5',
        'ROS_median','ROS_p2_5','ROS_p97_5',
        'Difference_median','Difference_p2_5','Difference_p97_5',
        'P_dual_gt_ros'
    ]
).to_csv('principal_uncertainty_effects.csv',index=False)

sobol_names=[
    'k_leak','k_pH','pH50','k_ROS','K_ROS','k_interaction',
    'k_MPS','k_plaque_in','cap_pH','core_ROS'
]

sb=np.array([
    [max(bounds[0,0],.75*selected_dual[0]),min(bounds[0,1],1.25*selected_dual[0])],
    [max(bounds[1,0],.75*selected_dual[1]),min(bounds[1,1],1.25*selected_dual[1])],
    [max(6.45,selected_dual[2]-.15),min(7.05,selected_dual[2]+.15)],
    [max(bounds[3,0],.75*selected_dual[3]),min(bounds[3,1],1.25*selected_dual[3])],
    [max(bounds[4,0],.65*selected_dual[4]),min(bounds[4,1],1.35*selected_dual[4])],
    [max(bounds[5,0],.50*selected_dual[5]),min(bounds[5,1],1.50*selected_dual[5])],
    [.05,.12],[.006,.015],[6.75,7.05],[.05,.20]
])

pr=first_primes(20)

def eval_sob(row):
    d=selected_dual.copy()
    d[:6]=row[:6]
    tr=transport.copy(); tr['k_mps']=row[6]
    sc=scenario.copy()
    sc['k_plaque_in']=row[7]
    sc['cap_pH']=row[8]
    sc['core_ROS_mM']=row[9]
    r=sim(d,'dual_interaction',sc,tr)
    b=sim(d,'baseline',sc,tr)
    return r[2]/max(b[2],np.finfo(float).eps)

rows=[]
for N in [512,1024,2048]:
    Au=halton(N,10,50,pr[:10])
    Bu=halton(N,10,5000,pr[10:20])
    A=sb[:,0]+Au*(sb[:,1]-sb[:,0])
    B=sb[:,0]+Bu*(sb[:,1]-sb[:,0])

    fA=np.array([eval_sob(x) for x in A])
    fB=np.array([eval_sob(x) for x in B])
    var=np.var(np.r_[fA,fB])

    S1=[]; ST=[]
    for j in range(10):
        AB=A.copy()
        AB[:,j]=B[:,j]
        fAB=np.array([eval_sob(x) for x in AB])
        S1.append(np.clip(1-np.mean((fB-fAB)**2)/(2*var),0,1))
        ST.append(np.clip(np.mean((fA-fAB)**2)/(2*var),0,1))

    order=np.argsort(ST)[::-1]
    for j in range(10):
        rows.append([
            N,sobol_names[j],S1[j],ST[j],
            int(np.where(order==j)[0][0]+1)
        ])

pd.DataFrame(
    rows,
    columns=['N','Parameter','First_Order','Total_Order','Rank']
).to_csv('convergence_sobol.csv',index=False)

print('Independent convergence and uncertainty reproduction complete.')
