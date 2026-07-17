"""
Eye-aging extension of the epigenetic-clock instance (paper Section 4.2), on a new tissue and modality.

Data: GSE314970, rat neural-retina bulk RNA-seq, 85 samples across 7 ages (6-27 months), one sex, one
tissue (Shavlakadze et al. 2025, Calico/Regeneron aging atlas). Human retinal aging data is single-cell
with too few donors (6-18) to compare clocks, so a well-powered rat bulk time-course is the tractable
substrate. This is transcriptomic, not methylation, so it EXTENDS the clock instance to a new tissue and
modality rather than adding a methylation clock.

Concept = transcriptomic aging state of the retina. Competing MEASURES (all oriented concept-positive,
higher = older / more aged):
  SenMayo         senescence-associated secretory signature (SAUL_SEN_MAYO)
  Fridman_up      FRIDMAN_SENESCENCE_UP
  Fridman_dn_neg  negated FRIDMAN_SENESCENCE_DN (down-with-senescence genes, flipped)
  Inflammatory    HALLMARK_INFLAMMATORY_RESPONSE (an inflammaging proxy; reused from the inflammation instance)
  Clock           elastic-net transcriptomic age predictor (out-of-fold cross_val_predict on age)
External anchor = chronological age in months (continuous). Deterministic (seeded CV and sampler).

Honest scope: a single dataset, so no cross-cohort reproducibility here (the study's RPE/choroid subset
is the obvious second cohort, not downloaded). Rat, one sex. Treated as a qualitative extension.
"""
import json
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr
from numpy.linalg import lstsq
from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_predict, KFold

import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/epigenetic-clock-desiderata/eye_aging/analysis/outputs/eye_clocks.txt"
class _Tee:
    def __init__(self, p): self._f = open(p, "w"); self._o = _sys.stdout
    def write(self, s): self._o.write(s); self._f.write(s)
    def flush(self): self._o.flush(); self._f.flush()
_sys.stdout = _Tee(_OUT)

D = "/Users/polina/Documents/BioInfStuff/epigenetic-clock-desiderata/eye_aging/data"
INFL = D  # hallmark_inflammatory.json copied locally so this instance is self-contained
counts = pd.read_parquet(f"{D}/GSE314970_rat_retina_counts.parquet")     # genes x samples
meta = pd.read_csv(f"{D}/GSE314970_rat_retina_metadata.csv")
meta = meta.set_index("gsm").reindex(counts.columns)
age = meta["age_months"].to_numpy(float)
print(f"=== Eye-aging (GSE314970 rat retina): {counts.shape[0]} genes x {counts.shape[1]} samples ===")
print(f"ages (months): {dict(pd.Series(age).value_counts().sort_index().astype(int))}\n")

# normalize: log2 CPM, then gene z-scores across samples
cpm = counts / counts.sum(0) * 1e6
logcpm = np.log2(cpm + 1.0)
Z = logcpm.sub(logcpm.mean(1), axis=0).div(logcpm.std(1).replace(0, np.nan), axis=0)
sym_upper = pd.Index([str(s).upper() for s in Z.index])                  # rat symbol -> upper for human-set matching

def load_msig(name):
    d = json.load(open(f"{D}/sig_{name}.json")); k = list(d)[0]
    return [g.upper() for g in d[k]["geneSymbols"]]
def load_hallmark():
    d = json.load(open(f"{INFL}/hallmark_inflammatory.json")); k = list(d)[0]
    return [g.upper() for g in d[k]["geneSymbols"]]
def score(genes):
    m = sym_upper.isin(set(genes))
    return Z.loc[m].mean(0).to_numpy(), int(m.sum()), len(genes)

sen, n1, t1 = score(load_msig("SAUL_SEN_MAYO"))
fup, n2, t2 = score(load_msig("FRIDMAN_SENESCENCE_UP"))
fdn, n3, t3 = score(load_msig("FRIDMAN_SENESCENCE_DN"))
infl, n4, t4 = score(load_hallmark())
print(f"signature gene coverage (matched/total): SenMayo {n1}/{t1}, Fridman_up {n2}/{t2}, "
      f"Fridman_dn {n3}/{t3}, Inflammatory {n4}/{t4}")

# fitted elastic-net transcriptomic clock, out-of-fold predictions on the top-variance genes
X = logcpm.loc[logcpm.var(1).sort_values().index[-2000:]].T.to_numpy()   # samples x 2000 genes
Xs = StandardScaler().fit_transform(X)
en = ElasticNetCV(l1_ratio=[.1,.5,.9], cv=5, max_iter=5000, random_state=0)
clock = cross_val_predict(en, Xs, age, cv=KFold(5, shuffle=True, random_state=0))
print(f"fitted clock: out-of-fold predicted age vs true age Spearman = {spearmanr(clock, age).correlation:.3f}\n")

# concept-positive measures (higher = older / more aged)
meas = {"SenMayo": sen, "Fridman_up": fup, "Fridman_dn_neg": -fdn, "Inflammatory": infl, "Clock": clock}
names = list(meas); M = np.column_stack([meas[n] for n in names]); N = len(age)

print("(families) Spearman among the 5 aging measures:")
S = pd.DataFrame(M, columns=names).corr("spearman")
print("               " + " ".join(f"{n[:9]:>9s}" for n in names))
for n in names: print(f"  {n:14s} " + " ".join(f"{S.loc[n,c]:+.2f}" for c in names))

print("\n(anchor) each measure vs chronological age (Spearman):")
for n in names: print(f"  {n:14s} {spearmanr(meas[n], age).correlation:+.3f}")
cons = np.column_stack([rankdata(M[:,k]) for k in range(len(names))]).mean(1)
print(f"  consensus (mean rank) vs age: {spearmanr(cons, age).correlation:+.3f}")
cp = (rankdata(cons)-1)/(N-1); q1,q2 = np.quantile(cp,[1/3,2/3])
print(f"  consensus vs age by tertile (mean age months): bottom {age[cp<=q1].mean():.1f} | "
      f"middle {age[(cp>q1)&(cp<q2)].mean():.1f} | top {age[cp>=q2].mean():.1f}")

print("\n(near-tie) pairwise discordance among the 5 measures vs consensus separation:")
Rk = np.column_stack([rankdata(M[:,k]) for k in range(len(names))])
i,j = np.triu_indices(N,1); s0 = np.sign(Rk[j,0]-Rk[i,0]); disc = np.zeros(len(i),bool)
for k in range(1,len(names)): disc |= (np.sign(Rk[j,k]-Rk[i,k]) != s0)
ax = cp; sep = np.abs(ax[i]-ax[j]); midpos = 1-2*np.abs((ax[i]+ax[j])/2-0.5)
print(f"  overall discordance {disc.mean():.3f}")
qs = np.quantile(sep,[0,.2,.4,.6,.8,1.0])
for a,b in zip(qs[:-1],qs[1:]):
    m=(sep>=a)&(sep<b if b!=qs[-1] else sep<=b)
    print(f"    |Δaxis| in [{a:.3f},{b:.3f}]: discordance={disc[m].mean():.3f}  (n={m.sum()})")
y=disc.astype(float); X0=np.column_stack([sep,np.ones_like(sep)]); X1=np.column_stack([sep,midpos,np.ones_like(sep)])
b0=lstsq(X0,y,rcond=None)[0]; b1=lstsq(X1,y,rcond=None)[0]
print(f"    linear-prob R^2: sep-only={1-((y-X0@b0).var()/y.var()):.4f}, "
      f"sep+midpos={1-((y-X1@b1).var()/y.var()):.4f}; midpos coef={b1[1]:+.4f}")

# consensus poset + average-rank canonical
geq = np.ones((N,N),bool)
for k in range(M.shape[1]): xk=M[:,k]; geq &= (xk[:,None]>=xk[None,:])
comp = geq|geq.T; np.fill_diagonal(comp,True)
print(f"\n(consensus poset) incomparable fraction = {1-(comp.sum()-N)/(N*(N-1)):.3f}")
def sar(seed,burn=400_000,nsamp=5000,gap=300):
    rng=np.random.default_rng(seed); order=list(np.argsort(-M.mean(1))); rs=np.zeros(N); ns=0
    for step in range(burn+nsamp*gap):
        a=rng.integers(N-1); x=order[a]; z=order[a+1]
        if not comp[x,z] and rng.random()<0.5: order[a],order[a+1]=z,x
        if step>=burn and (step-burn)%gap==0:
            for p,v in enumerate(order): rs[v]+=(p+1)
            ns+=1
    return rs/ns
ar1=sar(1); ar2=sar(2); avg=(ar1+ar2)/2
print(f"  Bubley-Dyer chains Spearman = {spearmanr(ar1,ar2).correlation:.4f}")
def dr(s): return rankdata(-s,method="average")
print("  Spearman(measure, canonical average-rank), and each measure vs age:")
for k,n in enumerate(names):
    print(f"    {n:14s} vs canonical {spearmanr(dr(M[:,k]),avg).correlation:+.3f}")
print(f"  canonical average-rank vs age: {spearmanr(-avg, age).correlation:+.3f}")

# save per-sample scores
out = meta[["age_months"]].copy()
for n in names: out[n]=meas[n]
out["consensus_rank"]=cons
out.to_csv(f"{D}/eye_scores.csv")
print("\nwrote eye_scores.csv")
print("\nHonest reading: several transcriptomic aging measures on one tissue disagree in the expected")
print("near-tie pattern, and the consensus tracks true age. Single cohort (no cross-cohort reproducibility),")
print("rat, one sex, and the clock is fitted on these same samples (out-of-fold) — a qualitative extension")
print("of the clock instance to a new tissue and modality, not a standalone domain.")
