#!/usr/bin/env python3
"""MCD + unweighted betaMNTD/betaNTI + incidence RCBray for two marker sets.
Definitions follow Danczak Generate_MCD / picante comdistnt / RCBC_PresAbs.
"""
from pathlib import Path
import numpy as np,pandas as pd,scipy,scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist,squareform
import json,hashlib,time,sys,os,itertools
B=Path(__file__).resolve().parent
D=next((p for p in [B.parent,B.parent.parent,Path('/workspace/DON')] if (p/'分子转化分析/全部样本_分子转化次数.csv').exists()),B.parent)
NREP=9999;BASE_SEED=20261107;TRAITS=['C','H','O','N','S','P','DBE','AI_Mod','kdefect']
PROC=['variable_selection','homogeneous_selection','homogenizing_dispersal','dispersal_limitation','drift']
CN={'variable_selection':'异质选择','homogeneous_selection':'同质选择','homogenizing_dispersal':'均质化扩散','dispersal_limitation':'扩散限制','drift':'漂变及其他未主导过程','undefined':'无法判定'}
read=lambda p:pd.read_csv(p,float_precision='round_trip')

def classify(z,rc):
 if not np.isfinite(z):return 'undefined'
 if z>2:return 'variable_selection'
 if z< -2:return 'homogeneous_selection'
 if not np.isfinite(rc):return 'undefined'
 if rc< -.95:return 'homogenizing_dispersal'
 if rc> .95:return 'dispersal_limitation'
 return 'drift'

def mntd_batch(dm,com,pairs):
 """picante abundance.weighted=FALSE: pooled mean of both directional minima.
 This is richness-weighted across directions, NOT 0.5*(two separate means).
 """
 single=dm.ndim==2
 if single:dm=dm[None,:,:]
 sizes=np.array([len(x) for x in com]);directional=np.empty((len(dm),len(com),len(com)))
 for j,inds in enumerate(com):
  nearest=dm[:,:,inds].min(axis=2)
  for i,indi in enumerate(com):directional[:,i,j]=nearest[:,indi].sum(axis=1)
 out=np.stack([(directional[:,i,j]+directional[:,j,i])/(sizes[i]+sizes[j]) for i,j in pairs],axis=1)
 return out[0] if single else out

def mntd_reference(dm,com,pairs):
 vals=[]
 for i,j in pairs:
  sub=dm[np.ix_(com[i],com[j])]
  vals.append(np.mean(np.r_[sub.min(axis=1),sub.min(axis=0)]))
 return np.array(vals)

def newick_from_Z(Z,labels,keep=None):
 # ape as.phylo.hclust divides hclust heights by two; pair path = cophenetic height.
 root=sch.to_tree(Z)
 def visit(node):
  if node.is_leaf():return None if keep is not None and node.id not in keep else (labels[node.id],0.)
  left=visit(node.left);right=visit(node.right);height=node.dist/2
  if left is None:return right
  if right is None:return left
  ls,lh=left;rs,rh=right
  return (f'({ls}:{height-lh:.15g},{rs}:{height-rh:.15g})',height)
 return visit(root)[0]+';'

def beta_null(Dfull,observed_indices,communities,pairs,nrep,seed,pool_mode,outpath):
 rng=np.random.default_rng(seed);nobs=len(observed_indices);pool=len(Dfull) if pool_mode=='full_marker_library' else nobs
 result=np.empty((nrep,len(pairs)));batch=96;t=time.time()
 for start in range(0,nrep,batch):
  bs=min(batch,nrep-start)
  perm=np.stack([rng.permutation(pool) for _ in range(bs)])
  mapped=perm[:,observed_indices] if pool_mode=='full_marker_library' else observed_indices[perm]
  dm=Dfull[mapped[:,:,None],mapped[:,None,:]]
  result[start:start+bs]=mntd_batch(dm,communities,pairs)
  if start==0:
   for k in range(min(bs,3)):np.testing.assert_allclose(result[k],mntd_reference(dm[k],communities,pairs),rtol=1e-13,atol=1e-13)
  if (start//batch)%25==0:print(outpath.parent.name,pool_mode,'beta null',min(start+bs,nrep),'/',nrep,'elapsed',round(time.time()-t,1),flush=True)
 np.savez_compressed(outpath,betaMNTD=result)
 return result

def rc_null(presence,pairs,nrep,seed,outpath):
 """Independent weighted sampling without replacement preserving each richness.
 Exponential races yield the same sequential PPS/Plackett-Luce sampling law
 as R sample(replace=FALSE, prob=occurrence). Compare integer shared counts,
 avoiding floating-point ties for binary Bray-Curtis.
 """
 rng=np.random.default_rng(seed);occ=presence.sum(axis=0).astype(float);assert (occ>0).all()
 richness=presence.sum(axis=1).astype(int);gamma=len(occ);null_shared=np.empty((nrep,len(pairs)),dtype=np.int16);observed=[];rc=[];time0=time.time()
 def draw(k):
  keys=rng.exponential(size=(nrep,gamma))/occ[None,:]
  inds=np.argpartition(keys,k-1,axis=1)[:,:k]
  out=np.zeros((nrep,gamma),dtype=bool);np.put_along_axis(out,inds,True,axis=1)
  assert (out.sum(axis=1)==k).all()
  return out
 for q,(i,j) in enumerate(pairs):
  a=draw(richness[i]);b=draw(richness[j]);shared=np.logical_and(a,b).sum(axis=1);obs=int(np.logical_and(presence[i],presence[j]).sum());observed.append(obs);null_shared[:,q]=shared
  value=2*(np.count_nonzero(shared>obs)+.5*np.count_nonzero(shared==obs))/nrep-1;rc.append(value)
  if q%10==0:print(outpath.parent.name,'RC pair',q+1,'/',len(pairs),'elapsed',round(time.time()-time0,1),flush=True)
 np.savez_compressed(outpath,shared=null_shared,observed_shared=np.array(observed),richness=richness)
 return np.array(rc),null_shared,np.array(observed)

def write_matrix(path,values,ids):pd.DataFrame(values,index=ids,columns=ids).to_csv(path,index_label='sample_id',encoding='utf-8-sig')
def pairs_matrix(v,ids,pairs,diagonal=np.nan):
 out=np.full((len(ids),len(ids)),np.nan);np.fill_diagonal(out,diagonal)
 for value,(i,j) in zip(v,pairs):out[i,j]=out[j,i]=value
 return out

fullpath=D/'分子转化分析/全部样本_分子转化次数.csv';raw=read(fullpath);samplemeta=read(D/'分子转化分析/样本汇总.csv');lake=raw[raw.endmember=='Lake'].copy();meta=samplemeta[samplemeta.endmember=='Lake'].sort_values('sample_id');ids=meta.sample_id.tolist();names=meta.set_index('sample_id')['sample'].to_dict();pairs=list(itertools.combinations(range(len(ids)),2));assert len(ids)==10 and len(pairs)==45
configs=[('RI_le1',D/'DON_marker_N大于0_HC小于1.5/候选marker_来源独有RI.csv'),('RI_eq0',D/'DON_marker_RI等于0_N大于0_HC小于1.5/RI等于0_候选marker.csv')]
all_summ=[];all_process=[];all_pairs=[];audit={'seed':BASE_SEED,'randomizations':NREP,'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'pandas':pd.__version__,'definitions':{'MCD':'z-score C,H,O,N,S,P,DBE,AI_Mod,kdefect then Euclidean + average-linkage UPGMA; zero-variance traits omitted','betaMNTD':'picante pooled concatenation of two directional nearest distances, abundance.weighted=FALSE, exclude.conspecifics=FALSE','betaNTI':'(observed-mean null)/sample SD null, ddof=1','beta_primary_null_pool':'marker formulas detected in at least one of the 10 Lake samples','beta_sensitivity_pool':'all source-exclusive candidates in the same RI marker library','RCBray':'binary Bray-Curtis; per-pair independent occurrence-weighted without-replacement sampling preserving richness; half-tie correction; no rounding before classification'},'input_sha256':{str(p.relative_to(D)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [fullpath]+[c[1] for c in configs]},'sets':{}}
for q,(setname,mpath) in enumerate(configs):
 out=B/setname;out.mkdir(exist_ok=True);markers=read(mpath).sort_values('formula').reset_index(drop=True);forms=markers.formula.tolist();markerids=markers.marker_id.tolist();n=len(markers)
 traits=markers[['marker_id','source','formula','C','H','O','N','S','P']].copy();masses=raw[raw.formula.isin(forms)].groupby('formula').mz.median();traits['representative_mz']=traits.formula.map(masses)
 C,H,O,N,S,P=[traits[k].to_numpy(float) for k in ['C','H','O','N','S','P']]
 traits['DBE']=1+C-.5*H+.5*N+.5*P;num=1+C-.5*O-S-.5*(N+P+H);den=C-.5*O-S-N-P
 traits['AI_Mod']=np.divide(num,den,out=np.zeros(n),where=(num>0)&(den>0));km=traits.representative_mz.to_numpy()*14/14.01565;traits['kdefect']=np.ceil(km)-km
 X=traits[TRAITS].to_numpy(float);assert np.isfinite(X).all();sd=X.std(axis=0,ddof=1);used=np.array(TRAITS)[sd>0].tolist();dropped=np.array(TRAITS)[sd==0].tolist();scaled=(X[:,sd>0]-X[:,sd>0].mean(axis=0))/sd[sd>0]
 traits.to_csv(out/'全部候选marker_MCD分子特征.csv',index=False,encoding='utf-8-sig');pd.DataFrame(scaled,index=markerids,columns=used).to_csv(out/'MCD特征Z分数.csv',index_label='marker_id',encoding='utf-8-sig');pd.DataFrame({'feature':TRAITS,'mean':X.mean(axis=0),'sample_SD':sd,'included':sd>0}).to_csv(out/'MCD特征标准化记录.csv',index=False,encoding='utf-8-sig')
 Z=sch.linkage(pdist(scaled,'euclidean'),method='average');Dfull=squareform(sch.cophenet(Z));assert np.isfinite(Dfull).all() and np.allclose(Dfull,Dfull.T);np.savez_compressed(out/'完整marker树距离与linkage.npz',distance=Dfull,linkage=Z,marker_ids=np.array(markerids),formulas=np.array(forms))
 incidence=pd.crosstab(lake.formula,lake.sample_id).reindex(index=forms,columns=ids,fill_value=0).gt(0).astype(int);observed_indices=np.flatnonzero(incidence.sum(axis=1).to_numpy()>0);obsforms=[forms[i] for i in observed_indices];obsids=[markerids[i] for i in observed_indices];presence=incidence.iloc[observed_indices].T.to_numpy(bool);assert (presence.sum(axis=1)>0).all();com=[np.flatnonzero(row) for row in presence]
 incidence.T.set_axis(markerids,axis=1).to_csv(out/'Lake_完整marker检出矩阵.csv',index_label='sample_id',encoding='utf-8-sig');pd.DataFrame(presence.astype(int),index=ids,columns=obsids).to_csv(out/'Lake_实际检出marker矩阵.csv',index_label='sample_id',encoding='utf-8-sig');traits.iloc[observed_indices].to_csv(out/'Lake检出marker_特征与身份.csv',index=False,encoding='utf-8-sig')
 sampleinfo=meta[['sample_id','sample']].copy();sampleinfo['marker_richness']=presence.sum(axis=1);sampleinfo.to_csv(out/'Lake样本marker检出数.csv',index=False,encoding='utf-8-sig')
 (out/'全部候选marker_MCD.nwk').write_text(newick_from_Z(Z,markerids));(out/'Lake检出marker_修剪MCD.nwk').write_text(newick_from_Z(Z,markerids,set(observed_indices)))
 Dobserved=Dfull[np.ix_(observed_indices,observed_indices)];pd.DataFrame(Dobserved,index=obsids,columns=obsids).to_csv(out/'Lake检出marker_MCD树距离.csv',index_label='marker_id',encoding='utf-8-sig')
 # independent cophenetic check: LCA heights equal patristic distances from exported tree construction.
 root=sch.to_tree(Z);parents={};heights={}
 def treewalk(node):
  heights[node.id]=node.dist/2
  if not node.is_leaf():
   for ch in [node.left,node.right]:parents[ch.id]=node.id;treewalk(ch)
 treewalk(root)
 for a,b in [(0,n-1),(n//3,n//2),(observed_indices[0],observed_indices[-1])]:
  ancestors=set();x=a
  while True:
   ancestors.add(x)
   if x not in parents:break
   x=parents[x]
  x=b
  while x not in ancestors:x=parents[x]
  assert np.isclose(Dfull[a,b],2*heights[x],atol=1e-12)
 obs=mntd_batch(Dobserved,com,pairs);np.testing.assert_allclose(obs,mntd_reference(Dobserved,com,pairs),rtol=1e-13,atol=1e-13)
 rc,rcshared,obsshared=rc_null(presence,pairs,NREP,BASE_SEED+100*q+3,out/'RCBray_9999次随机交集.npz')
 nprimary=beta_null(Dfull,observed_indices,com,pairs,NREP,BASE_SEED+100*q+1,'observed_Lake_markers',out/'betaMNTD_主分析_9999次.npz')
 nsens=beta_null(Dfull,observed_indices,com,pairs,NREP,BASE_SEED+100*q+2,'full_marker_library',out/'betaMNTD_完整来源库敏感性_9999次.npz')
 prefix_rows=[];pairframes={}
 for mode,nulls,poolsize in [('primary_observed_Lake_pool',nprimary,len(observed_indices)),('sensitivity_full_source_pool',nsens,n)]:
  nullmean=nulls.mean(axis=0);nullsd=nulls.std(axis=0,ddof=1);z=np.divide(obs-nullmean,nullsd,out=np.full(len(pairs),np.nan),where=nullsd>0);rows=[]
  for k,(i,j) in enumerate(pairs):
   bray=1-2*obsshared[k]/(presence[i].sum()+presence[j].sum());proc=classify(z[k],rc[k]);rcround=round(float(rc[k]),2)
   rows.append({'marker_set':setname,'null_pool':mode,'sample_i':ids[i],'name_i':names[ids[i]],'sample_j':ids[j],'name_j':names[ids[j]],'richness_i':int(presence[i].sum()),'richness_j':int(presence[j].sum()),'shared_markers':int(obsshared[k]),'observed_betaMNTD':obs[k],'null_mean_betaMNTD':nullmean[k],'null_SD_betaMNTD':nullsd[k],'betaNTI':z[k],'observed_Bray_Curtis':bray,'RCBray':rc[k],'RCBray_rounded2_reference_only':rcround,'process':proc,'process_Chinese':CN[proc],'stochastic':int(proc in PROC[2:]),'rounding_would_change_process':int(classify(z[k],rcround)!=proc),'betaNTI_distance_to_abs2':abs(abs(z[k])-2),'RCBray_distance_to_abs095':abs(abs(rc[k])-.95)})
  frame=pd.DataFrame(rows);pairframes[mode]=frame;prefix='主分析' if mode.startswith('primary') else '完整来源库敏感性';frame.to_csv(out/(prefix+'_45样本对结果.csv'),index=False,encoding='utf-8-sig');all_pairs.append(frame)
  write_matrix(out/(prefix+'_betaNTI矩阵.csv'),pairs_matrix(z,ids,pairs),ids);write_matrix(out/(prefix+'_betaMNTD矩阵.csv'),pairs_matrix(obs,ids,pairs,diagonal=0),ids)
  counts=frame.process.value_counts();valid=int((frame.process!='undefined').sum());stoch=int(frame.stochastic.sum());det=int(frame.process.isin(PROC[:2]).sum())
  sr={'marker_set':setname,'analysis':mode,'source_candidate_markers':n,'Lake_detected_markers':len(observed_indices),'null_pool_size':poolsize,'Lake_samples':len(ids),'sample_pairs':len(pairs),'valid_pairs':valid,'undefined_pairs':len(pairs)-valid,'randomizations':NREP,'stochastic_pairs':stoch,'stochastic_percent_of_valid':100*stoch/valid if valid else None,'deterministic_pairs':det,'deterministic_percent_of_valid':100*det/valid if valid else None,'mean_betaNTI':float(np.nanmean(z)),'min_betaNTI':float(np.nanmin(z)),'max_betaNTI':float(np.nanmax(z)),'rounding_sensitive_pairs':int(frame.rounding_would_change_process.sum()),'dropped_constant_traits':';'.join(dropped)};all_summ.append(sr)
  for proc in PROC+['undefined']:
   k=int(counts.get(proc,0));all_process.append({'marker_set':setname,'analysis':mode,'process':proc,'process_Chinese':CN[proc],'pairs':k,'percent_of_all_pairs':100*k/len(pairs)})
  for cutoff in [999,1999,4999,NREP]:
   nn=nulls[:cutoff];zz=(obs-nn.mean(axis=0))/nn.std(axis=0,ddof=1);rr=np.array([2*(np.count_nonzero(rcshared[:cutoff,k]>obsshared[k])+.5*np.count_nonzero(rcshared[:cutoff,k]==obsshared[k]))/cutoff-1 for k in range(len(pairs))]);pp=[classify(a,b) for a,b in zip(zz,rr)]
   for k,(i,j) in enumerate(pairs):prefix_rows.append({'analysis':mode,'iterations':cutoff,'sample_i':ids[i],'sample_j':ids[j],'betaNTI':zz[k],'RCBray':rr[k],'process':pp[k],'matches_9999_classification':int(pp[k]==frame.iloc[k].process)})
  print('RESULT',json.dumps(sr,ensure_ascii=False),flush=True)
 pd.DataFrame(prefix_rows).to_csv(out/'随机化次数收敛对照.csv',index=False,encoding='utf-8-sig');write_matrix(out/'RCBray矩阵.csv',pairs_matrix(rc,ids,pairs),ids)
 pa=pairframes['primary_observed_Lake_pool'];sa=pairframes['sensitivity_full_source_pool'];cmp=pa[['sample_i','sample_j','betaNTI','process','stochastic']].merge(sa[['sample_i','sample_j','betaNTI','process','stochastic']],on=['sample_i','sample_j'],suffixes=('_primary','_full_library'));cmp['same_process']=cmp.process_primary==cmp.process_full_library;cmp['both_stochastic']=(cmp.stochastic_primary==1)&(cmp.stochastic_full_library==1);cmp.to_csv(out/'零模型分子库敏感性逐对比较.csv',index=False,encoding='utf-8-sig')
 audit['sets'][setname]={'source_candidate_markers':n,'Lake_observed_markers':len(observed_indices),'traits_used':used,'constant_traits_omitted':dropped,'per_sample_richness':dict(zip(ids,presence.sum(axis=1).astype(int).tolist())),'primary_beta_seed':BASE_SEED+100*q+1,'full_library_beta_seed':BASE_SEED+100*q+2,'RC_seed':BASE_SEED+100*q+3,'null_pool_process_disagreements':int((~cmp.same_process).sum()),'both_pool_stochastic_pairs':int(cmp.both_stochastic.sum()),'verifications':['finite traits; z-score sample SD; zero-variance trait handling','UPGMA cophenetic=LCA patristic distance','observed betaMNTD vs independent per-pair implementation','first 3 randomized betaMNTD vs independent per-pair implementation for each pool','all RC randomized communities preserve observed richness','45 valid pair identities and 9999 null repetitions saved']}
summary=pd.DataFrame(all_summ);process=pd.DataFrame(all_process);pairresults=pd.concat(all_pairs,ignore_index=True)
summary.to_csv(B/'两套marker_保守性验证汇总.csv',index=False,encoding='utf-8-sig');process.to_csv(B/'两套marker_组装过程占比.csv',index=False,encoding='utf-8-sig');pairresults.to_csv(B/'两套marker_全部样本对结果.csv',index=False,encoding='utf-8-sig');meta[['sample_id','sample','source_file']].to_csv(B/'Lake样本身份表.csv',index=False,encoding='utf-8-sig');(B/'参数与数值核对.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2))
with pd.ExcelWriter(B/'DON_两套marker_湖泊保守性验证.xlsx',engine='openpyxl') as w:
 summary.to_excel(w,sheet_name='两套汇总',index=False);process.to_excel(w,sheet_name='过程占比',index=False);meta[['sample_id','sample','source_file']].to_excel(w,sheet_name='Lake样本身份',index=False)
 for setname,_ in configs:
  for title,file in [('主分析','主分析_45样本对结果.csv'),('全库敏感性','完整来源库敏感性_45样本对结果.csv'),('检出数','Lake样本marker检出数.csv'),('分子库对照','零模型分子库敏感性逐对比较.csv'),('迭代收敛','随机化次数收敛对照.csv'),('marker特征','全部候选marker_MCD分子特征.csv')]:read(B/setname/file).to_excel(w,sheet_name=setname+'_'+title,index=False)
print('COMPLETE',summary.to_string(index=False),flush=True)
