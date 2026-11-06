#!/usr/bin/env python3
"""Source-exclusive RI candidate markers; Lake is receptor, not a source."""
from pathlib import Path
from collections import defaultdict
import pandas as pd,numpy as np,json,hashlib,zipfile
B=Path(__file__).resolve().parent
D=next((p for p in [B.parent,B.parent.parent,Path('/workspace/DON')] if (p/'分子转化分析/全部样本_分子转化次数.csv').exists()),B.parent)
read=lambda p:pd.read_csv(p,float_precision='round_trip')
fullpath=D/'分子转化分析/全部样本_分子转化次数.csv';ripath=D/'RI_N大于0_HC小于1.5/全部筛后分子_RI标记.csv';summarypath=D/'分子转化分析/样本汇总.csv'
f=read(fullpath);r=read(ripath);samplemeta=read(summarypath)
source_names=sorted(set(r.endmember)-{'Lake'});assert len(source_names)==9
source_full=f[f.endmember.isin(source_names)].copy();source_ri=r[r.endmember.isin(source_names)].copy();lake=f[f.endmember=='Lake'].copy();lake_ri=r[r.endmember=='Lake'].copy()
# Original valid detected peaks, including all non-RI peaks, define source exclusivity.
support=source_full.groupby('formula').endmember.agg(lambda x:tuple(sorted(set(x)))).to_dict()
samples_per_source=samplemeta.groupby('endmember').sample_id.nunique().to_dict();lake_n=samples_per_source['Lake']
for k in ['C','H','O','N','P','S']:
 assert (f.groupby('formula')[k].nunique()==1).all(),('inconsistent formula',k)
assert (source_ri.N>0).all() and (source_ri.H/source_ri.C<1.5).all()
assert np.array_equal(r.RI_le1,(r.transformation_count<=1).astype(int))
source_groups={(s,form):x for (s,form),x in source_ri.groupby(['endmember','formula'])};lake_groups={form:x for form,x in lake.groupby('formula')}
rows=[];audits=[]
for (owner,formula),obs in source_groups.items():
 riobs=obs[obs.RI_le1==1]
 if riobs.empty:continue
 others=[s for s in support[formula] if s!=owner]
 audits.append({'source':owner,'formula':formula,'candidate_RI_samples':riobs.sample_id.nunique(),'detected_source_count':len(support[formula]),'detected_sources':';'.join(support[formula]),'other_sources_detected':';'.join(others),'source_exclusive':int(not others),'exclusion_reason':'detected_in_other_source' if others else ''})
 if others:continue
 assert support[formula]==(owner,)
 lo=lake_groups.get(formula);ls=0 if lo is None else lo.sample_id.nunique();first=obs.iloc[0];i=len(rows)+1
 record={'marker_id':f'M{i:04d}','source':owner,'formula':formula,**{e:int(first[e]) for e in ['C','H','O','N','P','S']},'H_C':float(first.H/first.C),'source_total_samples':samples_per_source[owner],'source_detected_samples':obs.sample_id.nunique(),'source_detected_fraction':obs.sample_id.nunique()/samples_per_source[owner],'source_RI_le1_samples':riobs.sample_id.nunique(),'source_RI_le1_sample_fraction':riobs.sample_id.nunique()/samples_per_source[owner],'source_peak_observations':len(obs),'source_RI_le1_peak_observations':len(riobs),'source_RI_eq0_samples':obs.loc[obs.RI_eq0==1,'sample_id'].nunique(),'source_min_count':int(obs.transformation_count.min()),'source_max_count':int(obs.transformation_count.max()),'source_mean_count':float(obs.transformation_count.mean()),'RI_le1_all_source_observations':int((obs.transformation_count<=1).all()),'RI_eq0_any_source_observation':int((obs.transformation_count==0).any()),'RI_eq0_all_source_observations':int((obs.transformation_count==0).all()),'source_at_least_2_samples':int(obs.sample_id.nunique()>=2),'source_mz_min':float(obs.mz.min()),'source_mz_median':float(obs.mz.median()),'source_mz_max':float(obs.mz.max()),'other_source_detected_count':0,'lake_total_samples':lake_n,'lake_detected_samples':ls,'lake_detected_fraction':ls/lake_n,'lake_detected':int(ls>0),'source_sample_ids':';'.join(sorted(obs.sample_id.unique())),'source_RI_le1_sample_ids':';'.join(sorted(riobs.sample_id.unique())),'lake_sample_ids':'' if lo is None else ';'.join(sorted(lo.sample_id.unique()))}
 rows.append(record)
markers=pd.DataFrame(rows).sort_values(['source','formula']).reset_index(drop=True);markers['marker_id']=[f'M{i:04d}' for i in range(1,len(markers)+1)]
audit=pd.DataFrame(audits).sort_values(['source','formula']).reset_index(drop=True)
assert markers.formula.is_unique
lake_detected=markers[markers.lake_detected==1].copy();strict=markers[markers.RI_eq0_any_source_observation==1].copy();allRI=markers[markers.RI_le1_all_source_observations==1].copy();replicate=markers[markers.source_at_least_2_samples==1].copy()
# Summary includes every source, including sources with no candidates.
summary=[]
for owner in source_names:
 a=audit[audit.source==owner];ms=markers[markers.source==owner]
 summary.append({'source':owner,'source_total_samples':samples_per_source[owner],'RI_candidate_formulas':len(a),'excluded_shared_with_other_sources':int((a.source_exclusive==0).sum()),'source_exclusive_RI_markers':len(ms),'markers_detected_in_Lake':int(ms.lake_detected.sum()),'markers_absent_in_Lake':int((ms.lake_detected==0).sum()),'strict_RI_eq0_any_markers':int(ms.RI_eq0_any_source_observation.sum()),'RI_le1_all_source_observations_markers':int(ms.RI_le1_all_source_observations.sum()),'detected_in_at_least_2_source_samples_markers':int(ms.source_at_least_2_samples.sum())})
summ=pd.DataFrame(summary)
# Presence of ALL valid detected formulas in 9 sources and Lake, not just RI.
presence=pd.crosstab(f.formula,f.endmember).reindex(columns=source_names+['Lake'],fill_value=0).gt(0).astype(int).reset_index()
source_map=markers.set_index('formula').source.to_dict();id_map=markers.set_index('formula').marker_id.to_dict();selected=set(markers.formula)
evidence=r[r.formula.isin(selected)].copy();evidence['marker_id']=evidence.formula.map(id_map);evidence['marker_source']=evidence.formula.map(source_map);evidence['role']=np.where(evidence.endmember=='Lake','receptor','source')
# A formula's composition fixes N and H/C, so qualifying marker presence in Lake agrees across full/screened inputs.
assert selected&set(lake.formula)==selected&set(lake_ri.formula)
source_ids=samplemeta.loc[samplemeta.endmember.isin(source_names),'sample_id'].tolist();lake_ids=samplemeta.loc[samplemeta.endmember=='Lake','sample_id'].tolist()
full_selected=f[f.formula.isin(selected)].copy();full_selected['marker_id']=full_selected.formula.map(id_map)
sample_presence=pd.crosstab(full_selected.marker_id,full_selected.sample_id).reindex(index=markers.marker_id,columns=source_ids+lake_ids,fill_value=0).gt(0).astype(int).reset_index()
sample_presence.insert(1,'source',sample_presence.marker_id.map(markers.set_index('marker_id').source));sample_presence.insert(2,'formula',sample_presence.marker_id.map(markers.set_index('marker_id').formula))
# Marker × receptor sample intensity; sum peak intensities when a formula has multiple peaks; raw, unnormalized.
intens=full_selected[full_selected.endmember=='Lake'].pivot_table(index='marker_id',columns='sample_id',values='intensity',aggfunc='sum',fill_value=0).reindex(index=markers.marker_id,columns=lake_ids,fill_value=0).reset_index()
intens.insert(1,'source',intens.marker_id.map(markers.set_index('marker_id').source));intens.insert(2,'formula',intens.marker_id.map(markers.set_index('marker_id').formula))
exports={'9来源端元汇总':summ,'候选marker_来源独有RI':markers,'Lake中检出的候选marker':lake_detected,'严格RI等于0_候选marker':strict,'全部来源检出记录均为RI_候选marker':allRI,'来源至少2样检出_候选marker':replicate,'RI候选排重全记录':audit,'9来源和Lake_全部分子式检出矩阵':presence,'候选marker_逐峰证据':evidence,'候选marker_60样本检出矩阵':sample_presence,'候选marker_Lake原始强度矩阵':intens,'样本身份表':samplemeta[['sample_id','endmember','sample','source_file']]}
for name,x in exports.items():x.to_csv(B/(name+'.csv'),index=False,encoding='utf-8-sig')
with pd.ExcelWriter(B/'DON_候选marker筛选结果.xlsx',engine='openpyxl') as x:
 for sheet,t in [('9来源汇总',summ),('来源独有RI候选marker',markers),('Lake检出marker',lake_detected),('严格RI等于0候选',strict),('来源全部记录均RI',allRI),('来源至少2样检出',replicate),('RI候选排重审计',audit),('60样本检出矩阵',sample_presence),('Lake原始强度',intens),('逐峰证据',evidence),('样本身份',exports['样本身份表'])]:t.to_excel(x,sheet_name=sheet,index=False)
 for owner in source_names:markers[markers.source==owner].to_excel(x,sheet_name=owner,index=False)
# Independent source-set exclusion check on raw source formulas, not using the support dictionary.
source_sets={s:set(source_full.loc[source_full.endmember==s,'formula']) for s in source_names}
for owner in source_names:
 ri_set=set(source_ri.loc[(source_ri.endmember==owner)&(source_ri.RI_le1==1),'formula']);other_union=set().union(*(source_sets[s] for s in source_names if s!=owner))
 expected=ri_set-other_union
 actual=set(markers.loc[markers.source==owner,'formula']);assert actual==expected
 assert set(lake_detected.loc[lake_detected.source==owner,'formula'])==expected&set(lake.formula)
assert summ.source_exclusive_RI_markers.sum()==len(markers) and summ.markers_detected_in_Lake.sum()==len(lake_detected)
assert (markers.N>0).all() and (markers.H_C<1.5).all() and markers.other_source_detected_count.eq(0).all()
assert len(summ)==9 and len(sample_presence.columns)-3==60
report={'sources':source_names,'receptor':'Lake','source_sample_count':len(source_ids),'receptor_sample_count':len(lake_ids),'RI_candidate_source_formula_combinations':len(audit),'shared_candidates_excluded':int((audit.source_exclusive==0).sum()),'source_exclusive_RI_markers':len(markers),'markers_detected_in_Lake':len(lake_detected),'markers_absent_in_Lake':len(markers)-len(lake_detected),'strict_RI_eq0_any_markers':len(strict),'RI_le1_all_source_observations_markers':len(allRI),'detected_in_at_least_2_source_samples_markers':len(replicate),'criteria':['source RI degree<=1 in at least one observed source peak','N>0 and H/C<1.5 before network construction','formula absent from all detected valid peaks in other 8 source types','Lake excluded from source-exclusivity comparison and used for receptor detection'],'independent_set_difference_check':'PASS for every source','limitations':'Candidate markers for the current sampled sources; no fixed prevalence threshold, statistical validation, random forest, XGBoost, SHAP or Bayesian source fractions inferred','input_sha256':{str(p.relative_to(D)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [fullpath,ripath,summarypath]}}
(B/'参数与核对记录.json').write_text(json.dumps(report,ensure_ascii=False,indent=2));print(json.dumps(report,ensure_ascii=False,indent=2));print(summ.to_string(index=False))
