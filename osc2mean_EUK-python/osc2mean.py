"""EU 与 Kaula 组合的瞬时根数转平根算法。

根数顺序：[a, e, i, RAAN, omega, M]，半长轴 km，角度 deg。
EU 元素顺序：[a, h, l, i, RAAN, lambda]，角度同样为 deg。
"""
from datetime import datetime, timezone
from math import sin, cos, tan, sqrt, atan2, acos, hypot, pi, factorial, ceil
from functools import lru_cache
import warnings
import numpy as np

MU = 398600.4418
RE = 6378.1363
J2 = 1.08262668e-3

def sind(x): return sin(np.deg2rad(x))
def cosd(x): return cos(np.deg2rad(x))
def _acos(x): return acos(np.clip(x, -1.0, 1.0))

def coe2EUK_deg(coe):
    """经典根数转 EU 元素：h=e*sin(omega)，l=e*cos(omega)，lambda=omega+M。"""
    a,e,i,O,w,M = np.asarray(coe).reshape(6)
    return np.array([a,e*sind(w),e*cosd(w),i,O,w+M])

def EUK2coe_deg(x):
    """由 EU 元素恢复经典根数，保留平近点角的连续表示。"""
    a,h,l,i,O,L = np.asarray(x).reshape(6)
    e = hypot(h,l); w = np.rad2deg(atan2(h,l)) if e > 1e-15 else 0.0
    return np.array([a,e,i,O % 360,w % 360,L-w])

def coe2rv_EU(coe, mu=MU):
    """求解开普勒方程，将经典根数转为位置 km、速度 km/s。"""
    a,e,inc,O,w,M = np.asarray(coe).reshape(6)
    Mrad = np.deg2rad(M); E = Mrad if e < .8 else pi
    for _ in range(100):
        dE = -(E-e*sin(E)-Mrad)/(1-e*cos(E)); E += dE
        if abs(dE) < 1e-14: break
    nu = 2*atan2(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2))
    p = a*(1-e*e); rmag = p/(1+e*cos(nu))
    rp = np.array([rmag*cos(nu),rmag*sin(nu),0.0])
    vp = sqrt(mu/p)*np.array([-sin(nu),e+cos(nu),0.0])
    i,O,w = np.deg2rad([inc,O,w]); ci,si=cos(i),sin(i)
    cO,sO,cw,sw=cos(O),sin(O),cos(w),sin(w)
    R = np.array([[cO*cw-sO*sw*ci,-cO*sw-sO*cw*ci,sO*si],
                  [sO*cw+cO*sw*ci,-sO*sw+cO*cw*ci,-cO*si],
                  [sw*si,cw*si,ci]])
    return R@rp,R@vp

def rv2coe_EU(r,v,mu=MU):
    """利用角动量与偏心率向量，由位置速度恢复椭圆轨道根数。"""
    r=np.asarray(r).reshape(3); v=np.asarray(v).reshape(3)
    R=np.linalg.norm(r); V=np.linalg.norm(v)
    hv=np.cross(r,v); h=np.linalg.norm(hv); inc=_acos(hv[2]/h)
    nv=np.cross([0.,0.,1.],hv); n=np.linalg.norm(nv)
    ev=np.cross(v,hv)/mu-r/R; e=np.linalg.norm(ev)
    a=-mu/(2*(V*V/2-mu/R))
    O=atan2(nv[1],nv[0]) % (2*pi) if n>1e-14 else 0.
    w=_acos(np.dot(nv,ev)/(n*e)) if e>1e-14 and n>1e-14 else 0.
    if e>1e-14 and n>1e-14 and ev[2]<0: w=2*pi-w
    if e>1e-14:
        nu=_acos(np.dot(ev,r)/(e*R))
        if np.dot(r,v)<0: nu=2*pi-nu
    elif n>1e-14:
        nu=_acos(np.dot(nv,r)/(n*R))
        if r[2]<0: nu=2*pi-nu
    else: nu=atan2(r[1],r[0]) % (2*pi)
    if e>=1: raise ValueError('Elliptic orbit required (e < 1).')
    E=2*atan2(sqrt(1-e)*sin(nu/2),sqrt(1+e)*cos(nu/2)) % (2*pi)
    M=E-e*sin(E)
    return np.array([a,e,np.rad2deg(inc),np.rad2deg(O)%360,
                     np.rad2deg(w)%360,np.rad2deg(M)%360])

def eu_J2_short_period(coeEU_x, opts):
    """J2 短周期改正，仅半长轴保留纯 J2 二阶项。"""
    a0      = coeEU_x[0]
    h0      = coeEU_x[1]
    l0      = coeEU_x[2]
    inc0    = coeEU_x[3]
    lambda0 = coeEU_x[5]
    Re = opts.get('Re', RE)
    J2v = opts.get('J2', J2)
    beta0 = sind(inc0)
    # J2 为正的扁率幅值，带符号的二阶带谐系数为 -J2。
    G20 = -J2v * (Re / a0)**2
    lambda_p = 1 - 3/2 * G20 * (3 - 4 * beta0)
    L = lambda0
    da1 = -3*a0*G20/(2 * lambda_p) * ( 
    (2 - 7/2 * beta0**2) * l0 * cosd(L) 
    + (2 - 5/2 * beta0**2) * h0 * sind(L) 
    + beta0**2 * cosd(2 * L) 
    + 7/2 * beta0**2 * (l0 * cosd(3 * L) + h0 * sind(3 * L)) 
    )
    da2 = 3 * a0 / 4 * G20**2 * beta0**2 * ( 
    7 * (2 - 3 * beta0**2) * cosd(2 * L) 
    + beta0**2 * cosd(4 * L) 
    )
    da = da1 + da2
    dh = -3 / (2 * lambda_p) * G20 * ( 
    (1 - 7/4 * beta0**2) * sind(L) 
    + (1 - 3 * beta0**2) * l0 * sind(2 * L) 
    + (-3/2 + 2 * beta0**2) * h0 * cosd(2 * L) 
    + 7/12 * beta0**2 * sind(3 * L) 
    + 17/8 * beta0**2 * (l0 * sind(4 * L) - h0 * cosd(4 * L)) 
    )
    dl = -3 / (2 * lambda_p) * G20 * ( 
    (1 - 5/4 * beta0**2) * cosd(L) 
    + 1/2 * (3 - 5 * beta0**2) * l0 * cosd(2 * L) 
    + (2 - 3/2 * beta0**2) * h0 * sind(2 * L) 
    + 7/12 * beta0**2 * cosd(3 * L) 
    + 17/8 * beta0**2 * (l0 * cosd(4 * L) + h0 * sind(4 * L)) 
    )
    # cos(i) 保留逆行轨道的符号；公式角度改正最后统一转为度。
    di = -3 / (4 * lambda_p) * G20 * beta0 * cosd(inc0) * ( 
    -l0 * cosd(L) + h0 * sind(L) + cosd(2 * L) 
    + 7/3 * (l0 * cosd(3 * L) + h0 * sind(3 * L)) 
    )
    dOmega = 3 / (2 * lambda_p) * G20 * cosd(inc0) * ( 
    7/2 * l0 * sind(L) - 5/2 * h0 * cosd(L) - 1/2 * sind(2 * L) 
    - 7/6 * l0 * sind(3 * L) + 7/6 * h0 * cosd(3 * L) 
    )
    dlambda = -3 / (2 * lambda_p) * G20 * ( 
    (10 - 119/8 * beta0**2) * l0 * sind(L) 
    + (85/8 * beta0**2 - 9) * h0 * cosd(L) 
    + (2 * beta0**2 - 1/2) * sind(2 * L) 
    + (-7/6 + 119/24 * beta0**2) * (l0 * sind(3 * L) - h0 * cosd(3 * L)) 
    - (3 - 21/4 * beta0**2) * l0 * sind(L) 
    + (3 - 15/4 * beta0**2) * h0 * cosd(L) 
    - 3/4 * beta0**2 * sind(2 * L) 
    - 21/12 * beta0**2 * (l0 * sind(3 * L) - h0 * cosd(3 * L)) 
    )
    return np.array([da, dh, dl, np.rad2deg(di), np.rad2deg(dOmega), np.rad2deg(dlambda)])


def osc2mean_EU(coe_osc, epoch=None, opts=None):
    """迭代求 J2 平根，以重构瞬时位置速度的残差判断收敛。"""
    opts=dict(opts or {}); mu=opts.get('mu',MU)
    opts.setdefault('Re',RE); opts.setdefault('J2',J2)
    max_iter=opts.get('maxIter',20); tol_r=opts.get('tol_r',1e-9); tol_v=opts.get('tol_v',1e-12)
    if not isinstance(max_iter,int) or max_iter<1: raise ValueError('maxIter must be a positive integer')
    r0,v0=coe2rv_EU(coe_osc,mu); x0=coe2EUK_deg(coe_osc); x=x0.copy()
    converged=False
    for it in range(1,max_iter+1):
        r,v=coe2rv_EU(EUK2coe_deg(x+eu_J2_short_period(x,opts)),mu)
        dr,dv=r0-r,v0-v; er,ev=np.linalg.norm(dr),np.linalg.norm(dv)
        current=EUK2coe_deg(x)
        if opts.get('verbose',False): print(f'Iter {it:3d}: dr={er:.12e} km, dv={ev:.12e} km/s')
        converged=er<=tol_r and ev<=tol_v
        if converged or it==max_iter: break
        # 在当前平根对应的位置速度上施加残差，得到下一轮平根。
        rm,vm=coe2rv_EU(current,mu)
        x=coe2EUK_deg(rv2coe_EU(rm+dr,vm+dv,mu))
    return EUK2coe_deg(x), dict(iter=it, r_error=er, v_error=ev,
                               converged=bool(converged))

def read_GGM03C_gfc(filename,Lmax=60):
    """读取指定阶次的归一化重力系数；数组直接以阶 l、次 m 索引。"""
    C=np.zeros((Lmax+1,Lmax+1)); S=np.zeros_like(C); count=0
    with open(filename,encoding='utf-8') as file:
        for line in file:
            t=line.split()
            if len(t)<5 or t[0].lower()!='gfc': continue
            l,m=int(t[1]),int(t[2])
            if not 0<=m<=l<=Lmax: continue
            C[l,m]=float(t[3].replace('D','E')); S[l,m]=float(t[4].replace('D','E')); count+=1
    if not count: raise ValueError(f'No gravity coefficients in {filename}')
    return C,S

def gmst_EU(epoch):
    """格林尼治平恒星时；无时区日期按原时刻使用，不自动平移时区。"""
    if isinstance(epoch,str): epoch=datetime.fromisoformat(epoch)
    if epoch.tzinfo is not None: epoch=epoch.astimezone(timezone.utc).replace(tzinfo=None)
    JD=2451545.0+(epoch-datetime(2000,1,1,12)).total_seconds()/86400
    T=(JD-2451545.0)/36525
    sec=(67310.54841+(876600*3600+8640184.812866)*T+.093104*T*T-6.2e-6*T**3)%86400
    return 2*pi*sec/86400

def KaulaF_EU(l,m,p,i):
    """计算倾角函数及其对倾角的导数，输入倾角为弧度。"""
    if l<2 or m<0 or m>l or p<0 or p>l: return 0.,0.
    s,c=sin(i/2),cos(i/2)
    factor=factorial(l+m)/(2**l*factorial(p)*factorial(l-p))
    F=dF=0.
    for j in range(max(0,l-2*p-m),min(l-m,2*l-2*p)+1):
        if l-m-j<0 or 2*p-l+m+j<0: continue
        bin1=factorial(2*l-2*p)/(factorial(j)*factorial(2*l-2*p-j))
        bin2=factorial(2*p)/(factorial(l-m-j)*factorial(2*p-l+m+j))
        es=m-l+2*p+2*j; ec=3*l-m-2*p-2*j
        term=(-1)**j*bin1*bin2*s**es*c**ec
        F+=term; dF+=term*(es/2/tan(i/2)-ec/2*tan(i/2))
    sign=(-1)**ceil((l-m)/2)
    return sign*factor*F,sign*factor*dF

def KaulaFTable_EU(Lmax,inc_deg):
    """预计算各阶次的倾角函数表，输入倾角为度。"""
    F=np.zeros((Lmax+1,)*3); dF=np.zeros_like(F)
    for l in range(2,Lmax+1):
        for m in range(l+1):
            for p in range(l+1): F[l,m,p],dF[l,m,p]=KaulaF_EU(l,m,p,np.deg2rad(inc_deg))
    return F,dF

@lru_cache(maxsize=None)
def _binomial(a,n):
    """允许上标为负整数的广义二项式系数。"""
    if n<0 or (a>=0 and n>a): return 0.
    value=1.
    for j in range(1,n+1): value=value*(a-j+1)/j
    return value

def _G_single(l,p,e,q,N):
    """用截断级数计算单个 Hansen 偏心率系数。"""
    if e==0: return float(q==0)
    root=sqrt(1-e*e); beta=e/(1+root); beta2=beta*beta
    pp,qq=(p,q) if p<=l/2 else (l-p,-q)
    A=.5*(l-2*pp+qq)*(1+root); total=0.
    for k in range(N):
        hp,hq=(k+qq,k) if qq>0 else (k,k-qq)
        P=Q=0.
        for r in range(hp+1):
            b=_binomial(2*pp-2*l,hp-r)
            if b: P+=b*(-1)**r/factorial(r)*A**r
        for r in range(hq+1):
            b=_binomial(-2*pp,hq-r)
            if b: Q+=b/factorial(r)*A**r
        term=P*Q*beta2**k; total+=term
        if k>5 and abs(term)<1e-15*max(1,abs(total)): break
    G=(-1)**abs(q)*(1+beta2)**l*beta**abs(q)*total
    return 0. if abs(G)<1e-15 else G

def KaulaG_EU(l,p,e,qMax=5,N=50,derivative_step=None):
    """计算 q=-qMax 到 qMax 的偏心率函数及其数值导数。"""
    for value,label in [(l,'l'),(p,'p'),(qMax,'qMax'),(N,'N')]:
        if int(value)!=value or value<0: raise ValueError(f'{label} must be a nonnegative integer')
    if p>l or N<1 or not 0<=e<1: raise ValueError('Invalid Hansen coefficient arguments')
    g=np.zeros(2*qMax+1); dg=np.zeros_like(g)
    for iq,q in enumerate(range(-qMax,qMax+1)):
        g[iq]=_G_single(l,p,e,q,N)
        if e==0:
            dg[iq]=0. if q==0 else _G_single(l,p,1e-7,q,N)/1e-7
        else:
            h=derivative_step if derivative_step is not None else max(1e-8,1e-5*max(e,1e-3))
            if not np.isfinite(h) or h<=0: raise ValueError('Derivative step must be positive and finite')
            ep,em=min(e+h,1-1e-10),max(e-h,0)
            dg[iq]=(_G_single(l,p,ep,q,N)-_G_single(l,p,em,q,N))/(ep-em)
    return g,dg

def KaulaShortPeriod(coe,epoch,opts):
    """计算除 J2 外的短周期和中周期改正，返回 km 与 deg。"""
    a,e,inc,O,w,M=np.asarray(coe).reshape(6)
    mu,Re,J2v=opts.get('mu',MU),opts.get('Re',RE),opts.get('J2',J2)
    if not 0<e<1 or abs(sind(inc))<1e-12:
        raise ValueError('Classical Kaula requires 0 < e < 1 and non-equatorial inclination')
    Lmax=opts.get('Lmax',opts['C'].shape[0]-1); qMax=opts.get('qMax',5); N=opts.get('gN',8)
    C=np.array(opts['C'],copy=True,dtype=float); S=np.array(opts['S'],copy=True,dtype=float)
    if min(C.shape)<Lmax+1 or min(S.shape)<Lmax+1: raise ValueError('Gravity matrices smaller than Lmax+1')
    # 将归一化重力系数转换为本级数所使用的非归一化系数。
    if opts.get('gravityNormalized',True):
        for l in range(Lmax+1):
            for m in range(l+1):
                norm=sqrt((2*l+1)*(2-(m==0))*factorial(l-m)/factorial(l+m))
                C[l,m]*=norm; S[l,m]*=norm
    theta=gmst_EU(epoch); i=np.deg2rad(inc); n=sqrt(mu/a**3); p_orbit=a*(1-e*e)
    b=sqrt(1-e*e)
    # 用 J2 长期进动率计算各谐项的相位频率。
    Od=-1.5*J2v*n*(Re/p_orbit)**2*cos(i)
    wd=.75*J2v*n*(Re/p_orbit)**2*(4-5*sin(i)**2)
    Md=n-.75*J2v*n*(Re/p_orbit)**2*b*(3*sin(i)**2-2)
    td=7.29211514670698e-5
    F,dF=KaulaFTable_EU(Lmax,inc)
    step=opts.get('gDerivativeStep',max(1e-8,1e-4*max(e,1e-3)))
    G={}; dG={}
    for l in range(2,Lmax+1):
        for p in range(l+1): G[l,p],dG[l,p]=KaulaG_EU(l,p,e,qMax,N,step)
    total=np.zeros(6); denom=b*sin(i)
    O,w,M=np.deg2rad([O,w,M])
    for l in range(2,Lmax+1):
        ca=mu*Re**l/(n*a**(l+2)); ci=mu*Re**l/(n*a**(l+3)); ce=ci/e
        for m in range(l+1):
            # J2 已由 EU 阶段处理，避免重复改正。
            if (l==2 and m==0) or (C[l,m]==0 and S[l,m]==0): continue
            for p in range(l+1):
                f,df=F[l,m,p],dF[l,m,p]; k=l-2*p
                for iq,q in enumerate(range(-qMax,qMax+1)):
                    s=k+q
                    # 保留中周期项，跳过带谐的长期及长周期项。
                    if s==0 and m==0: continue
                    freq=k*wd+s*Md+m*(Od-td)
                    if abs(freq)<1e-16: continue
                    g,dg=G[l,p][iq],dG[l,p][iq]
                    phase=k*w+s*M+m*(O-theta)
                    cp,sp=cos(phase),sin(phase)
                    # st=-dS/dphase，保证角根数改正的积分符号正确。
                    if (l-m)%2==0: sv=C[l,m]*cp+S[l,m]*sp; st=C[l,m]*sp-S[l,m]*cp
                    else: sv=-S[l,m]*cp+C[l,m]*sp; st=-S[l,m]*sp-C[l,m]*cp
                    total[0]+=ca*2*f*g*s*sv/freq
                    total[1]+=ce*b*f*g*(b*s-k)*sv/freq
                    total[2]+=ci/denom*f*g*(k*cos(i)-m)*sv/freq
                    total[3]+=ci/denom*df*g*st/freq
                    total[4]+=ci*(b/e*f*dg-1/tan(i)/b*df*g)*st/freq
                    total[5]+=ci*(-(1-e*e)/e*dg+2*(l+1)*g)*f*st/freq
    total[2:]=np.rad2deg(total[2:])
    return total

def osc2mean_EUK(coe_osc,epoch,opts):
    """返回 EUK 平根及诊断信息，其中 coe_J2mean 为 EU 阶段结果。"""
    mean,info=osc2mean_EU(coe_osc,epoch,opts)
    if not info['converged']: warnings.warn('J2 inversion did not converge',RuntimeWarning)
    dk=KaulaShortPeriod(mean,epoch,opts); result=mean-dk; result[3:]%=360
    info.update(coe_J2mean=mean,dKaula=dk)
    return result,info

