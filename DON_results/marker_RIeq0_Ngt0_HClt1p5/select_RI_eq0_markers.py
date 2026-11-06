from pathlib import Path
from collections import defaultdict
import pandas as pd,numpy as np,json,hashlib,zipfile
B=Path(__file__).resolve().parent
D=next((p for p in [B.parent,B.parent.parent,Path('/workspace/DON')] if (p/'分子转化分析/全部样本_分子转化次数.csv').exists()),B.parent)
read=lambda p:pd.read_csv(p,float_precision='round_trip')
fullpath=D/'分子转化分析/全部样本_分子转化次数.csv';ripath=D/'RI_N大于0_HC小于1.5/全部筛后分子_RI标记.csv';summarypath=D/'分子转化分析/样本汇总.csv';previouspath=D/'DON_marker_N大于0_HC小于1.5/候选marker_来源独有RI.csv'
f=read(fullpath);r=read(ripath);samplemeta=read(summarypath);previous=read(previouspath)
sources=sorted(set(r.endmember)-{'Lake'});sourcefull=f[f.endmember.isin(sources)];lake=f[f.endmember=='Lake'];sourcefiltered=r[r.endmember.isin(sources)];zero=r[r.transformation_count==0].copy()
assert (r.N>0).all() and (r.H/r.C<1.5).all() and np.array_equal(r.RI_eq0,(r.transformation_count==0).astype(int))
source_sets={s:set(sourcefull.loc[sourcefull.endmember==s,'formula']) for s in sources};nsource=samplemeta.groupby('endmember').sample_id.nunique().to_dict()
prevmap=previous.set_index('formula').to_dict('index');rows=[];audit=[]
for (owner,form),obs in sourcefiltered.groupby(['endmember','formula']):
 z=obs[obs.transformation_count==0]
 if z.empty:continue
 others=[s for s in sources if s!=owner and form in source_sets[s]]
 audit.append({'source':owner,'formula':form,'source_RI_eq0_samples':z.sample_id.nunique(),'source_RI_eq0_peak_observations':len(z),'other_source_detected_count':len(others),'other_sources_detected':';'.join(others),'source_exclusive':int(not others),'exclusion_reason':'detected_in_other_source' if others else ''})
 if others:continue
 assert form in prevmap and prevmap[form]['source']==owner
 lc=lake[lake.formula==form];first=obs.iloc[0]
 rows.append({'marker_id':prevmap[form]['marker_id'],'source':owner,'formula':form,**{e:int(first[e]) for e in ['C','H','O','N','P','S']},'H_C':float(first.H/first.C),'RI_threshold':0,'source_total_samples':nsource[owner],'source_detected_samples':obs.sample_id.nunique(),'source_detected_fraction':obs.sample_id.nunique()/nsource[owner],'source_RI_eq0_samples':z.sample_id.nunique(),'source_RI_eq0_sample_fraction':z.sample_id.nunique()/nsource[owner],'source_peak_observations':len(obs),'source_RI_eq0_peak_observations':len(z),'source_min_count':int(obs.transformation_count.min()),'source_max_count':int(obs.transformation_count.max()),'source_mean_count':float(obs.transformation_count.mean()),'RI_eq0_all_source_observations':int((obs.transformation_count==0).all()),'source_detected_at_least_2_samples':int(obs.sample_id.nunique()>=2),'source_RI_eq0_at_least_2_samples':int(z.sample_id.nunique()>=2),'source_mz_min':float(obs.mz.min()),'source_mz_median':float(obs.mz.median()),'source_mz_max':float(obs.mz.max()),'other_source_detected_count':0,'lake_total_samples':nsource['Lake'],'lake_detected_samples':lc.sample_id.nunique(),'lake_detected_fraction':lc.sample_id.nunique()/nsource['Lake'],'lake_detected':int(not lc.empty),'source_sample_ids':';'.join(sorted(obs.sample_id.unique())),'source_RI_eq0_sample_ids':';'.join(sorted(z.sample_id.unique())),'lake_sample_ids':';'.join(sorted(lc.sample_id.unique()))})
m=pd.DataFrame(rows).sort_values(['source','formula']).reset_index(drop=True);audit=pd.DataFrame(audit).sort_values(['source','formula']);assert m.formula.is_unique
sel=set(m.formula);lakem=m[m.lake_detected==1].copy();allzero=m[m.RI_eq0_all_source_observations==1].copy();repeat=m[m.source_RI_eq0_at_least_2_samples==1].copy()
summ=[]
for owner in sources:
 a=audit[audit.source==owner];x=m[m.source==owner]
 summ.append({'source':owner,'source_total_samples':nsource[owner],'RI_eq0_candidate_formulas':len(a),'excluded_shared_with_other_sources':int((a.source_exclusive==0).sum()),'RI_eq0_source_exclusive_markers':len(x),'markers_detected_in_Lake':int(x.lake_detected.sum()),'all_source_observations_eq0_markers':int(x.RI_eq0_all_source_observations.sum()),'RI_eq0_in_at_least_2_source_samples_markers':int(x.source_RI_eq0_at_least_2_samples.sum()),'previous_RI_le1_markers':int((previous.source==owner).sum()),'removed_from_previous_markers':int((previous.source==owner).sum())-len(x)})
summ=pd.DataFrame(summ)
idmap=m.set_index('formula').marker_id.to_dict();ownmap=m.set_index('formula').source.to_dict()
evidence=r[r.formula.isin(sel)].copy();evidence['marker_id']=evidence.formula.map(idmap);evidence['marker_source']=evidence.formula.map(ownmap);evidence['role']=np.where(evidence.endmember=='Lake','receptor','source');evidence['is_source_RI_eq0_evidence']=((evidence.endmember==evidence.marker_source)&(evidence.transformation_count==0)).astype(int)
fullsel=f[f.formula.isin(sel)].copy();fullsel['marker_id']=fullsel.formula.map(idmap);ids=samplemeta.sample_id.tolist();lakeids=samplemeta.loc[samplemeta.endmember=='Lake','sample_id'].tolist()
presence=pd.crosstab(fullsel.marker_id,fullsel.sample_id).reindex(index=m.marker_id,columns=ids,fill_value=0).gt(0).astype(int).reset_index()
intens=fullsel[fullsel.endmember=='Lake'].pivot_table(index='marker_id',columns='sample_id',values='intensity',aggfunc='sum',fill_value=0).reindex(index=m.marker_id,columns=lakeids,fill_value=0).reset_index()
for t in [presence,intens]:t.insert(1,'source',t.marker_id.map(m.set_index('marker_id').source));t.insert(2,'formula',t.marker_id.map(m.set_index('marker_id').formula))
comparison=previous[['marker_id','source','formula','lake_detected']].copy();comparison['retained_under_RI_eq0']=comparison.formula.isin(sel).astype(int)
exports={'RI等于0_候选marker':m,'Lake检出_RI等于0候选marker':lakem,'全部来源记录均等于0_候选marker':allzero,'来源至少2样RI等于0_候选marker':repeat,'9来源汇总':summ,'RI等于0候选_排重全记录':audit,'候选marker_逐峰证据':evidence,'候选marker_60样本检出矩阵':presence,'候选marker_Lake原始强度矩阵':intens,'与RI小于等于1版对照':comparison,'全部样本_RI等于0分子峰':zero}
for name,t in exports.items():t.to_csv(B/(name+'.csv'),index=False,encoding='utf-8-sig')
with pd.ExcelWriter(B/'DON_RI等于0_候选marker结果.xlsx',engine='openpyxl') as x:
 for sheet,t in [('9来源汇总',summ),('RI等于0候选marker',m),('Lake检出marker',lakem),('来源全部记录均为0',allzero),('来源至少2样为0',repeat),('RI等于0排重审计',audit),('60样本检出矩阵',presence),('Lake原始强度',intens),('逐峰证据',evidence),('与RI小于等于1对照',comparison),('全部样本RI等于0分子峰',zero)]:t.to_excel(x,sheet_name=sheet,index=False)
 for owner in sources:m[m.source==owner].to_excel(x,sheet_name=owner,index=False)
# Independent set subtraction from the full source background, plus exact comparison with previous strict flags.
for owner in sources:
 zforms=set(sourcefiltered.loc[(sourcefiltered.endmember==owner)&(sourcefiltered.transformation_count==0),'formula']);background=set(sourcefull.loc[sourcefull.endmember!=owner,'formula'])
 assert set(m.loc[m.source==owner,'formula'])==zforms-background
 assert set(lakem.loc[lakem.source==owner,'formula'])==(zforms-background)&set(lake.formula)
assert sel==set(previous.loc[previous.RI_eq0_any_source_observation==1,'formula'])
assert summ.RI_eq0_source_exclusive_markers.sum()==len(m) and summ.markers_detected_in_Lake.sum()==len(lakem)
assert (m.source_min_count==0).all() and (m.N>0).all() and (m.H_C<1.5).all()
report={'RI_primary_rule':'transformation_count == 0 in at least one detected source peak','N_HC_filters':'N>0 and H/C<1.5 before constructing the existing per-sample network','sources':sources,'receptor':'Lake','RI_eq0_peak_observations_all_60_samples':len(zero),'RI_eq0_unique_formulas_all_60_samples':zero.formula.nunique(),'source_RI_eq0_candidate_combinations':len(audit),'source_shared_candidates_excluded':int((audit.source_exclusive==0).sum()),'source_exclusive_RI_eq0_markers':len(m),'markers_detected_in_Lake':len(lakem),'all_source_observations_eq0_markers':len(allzero),'source_RI_eq0_at_least_2_samples_markers':len(repeat),'source_detected_at_least_2_samples_markers':int(m.source_detected_at_least_2_samples.sum()),'previous_RI_le1_markers':len(previous),'removed_from_previous_markers':len(previous)-len(m),'checks':'PASS: all 9 independent source set differences; Lake intersections; exact subset of previous strict RI flags; summary counts; N and HC conditions','input_sha256':{str(p.relative_to(D)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [fullpath,ripath,summarypath,previouspath]}}
(B/'参数与核对记录.json').write_text(json.dumps(report,ensure_ascii=False,indent=2));print(json.dumps(report,ensure_ascii=False,indent=2));print(summ.to_string(index=False))
