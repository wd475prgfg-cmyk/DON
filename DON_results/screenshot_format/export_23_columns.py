from pathlib import Path
import csv,collections,json,zipfile,shutil,math
from openpyxl import Workbook,load_workbook
from openpyxl.styles import Font,PatternFill
from openpyxl.cell import WriteOnlyCell
from openpyxl.utils import get_column_letter
B=Path(__file__).resolve().parent;D=B.parent
FIELDS=['sample','peak.x','C_x','H_x','O_x','N_x','P_x','S_x','formula_x','class_x','peak.y','C_y','H_y','O_y','N_y','P_y','S_y','formula_y','class_y','Dist','Dist.plus','Dist.minus','Trans.name']
configs=[('最新_N大于0_HC小于1.5',D/'RI_N大于0_HC小于1.5','全部筛后分子_RI标记.csv',205000),('旧版_未加N和HC条件',D/'分子转化分析','全部样本_分子转化次数.csv',1215666)]
def reader(p):
 with p.open(encoding='utf-8-sig',newline='') as f:yield from csv.DictReader(f)
def num(v):
 x=float(v);return int(x) if x.is_integer() else x
summaries=[];checks=[]
for version,A,nodename,expected_total in configs:
 out=B/version;(out/'逐样本CSV').mkdir(parents=True,exist_ok=True);(out/'10端元CSV').mkdir(exist_ok=True)
 nodes=collections.defaultdict(dict)
 for r in reader(A/nodename):nodes[r['sample_id']][r['source_row']]=r
 writers={};handles=[];worksheets={};counts=collections.Counter();total=0
 allfile=(out/'全部样本_23列转化明细.csv').open('w',encoding='utf-8-sig',newline='');allwriter=csv.writer(allfile);allwriter.writerow(FIELDS);handles.append(allfile)
 wb=Workbook(write_only=True) if version.startswith('最新') else None
 index=[]
 for f in sorted((A/'逐样本转化峰对').glob('*.csv')):
  sid=f.name[:3];sample_nodes=nodes[sid];n=0
  with (out/'逐样本CSV'/f.name).open('w',encoding='utf-8-sig',newline='') as dest:
   sw=csv.writer(dest);sw.writerow(FIELDS)
   for e in reader(f):
    x,y=sample_nodes[e['source_row_x']],sample_nodes[e['source_row_y']];group=e['endmember'];sample='Sample_'+sid+'_'+group+'_'+e['sample']
    assert x['formula']==e['formula_x'] and y['formula']==e['formula_y']
    assert float(x['mz'])==float(e['mz_x']) and float(y['mz'])==float(e['mz_y'])
    delta=float(e['mass_difference_Da']);fmt=lambda v:format(v,'.10f').rstrip('0').rstrip('.')
    row=[sample,e['mz_x'],*[int(float(x[k])) for k in ['C','H','O','N','P','S']],x['formula'],x['class'],e['mz_y'],*[int(float(y[k])) for k in ['C','H','O','N','P','S']],y['formula'],y['class'],fmt(delta),fmt(delta+1e-5),fmt(delta-1e-5),e['transformation']]
    assert len(row)==23
    if version.startswith('最新'):
     assert row[5]>0 and row[14]>0 and row[3]/row[2]<1.5 and row[12]/row[11]<1.5
    if group not in writers:
     h=(out/'10端元CSV'/(group+'.csv')).open('w',encoding='utf-8-sig',newline='');handles.append(h);gw=csv.writer(h);gw.writerow(FIELDS);writers[group]=gw
     if wb:
      ws=wb.create_sheet(group);ws.freeze_panes='A2'
      for i,name in enumerate(FIELDS,1):ws.column_dimensions[get_column_letter(i)].width=46 if i==1 else 22 if name in ['formula_x','formula_y','Trans.name'] else 15
      head=[]
      for name in FIELDS:
       c=WriteOnlyCell(ws,value=name);c.font=Font(bold=True,color='FFFFFF');c.fill=PatternFill('solid',fgColor='305496');head.append(c)
      ws.append(head);worksheets[group]=ws
    sw.writerow(row);allwriter.writerow(row);writers[group].writerow(row)
    if wb:
     values=[float(v) if i in [1,10,19,20,21] else v for i,v in enumerate(row)]
     worksheets[group].append(values)
    n+=1;total+=1;counts[group]+=1
  index.append([sid,f.stem,n]);print(version,sid,n,flush=True)
 for h in handles:h.close()
 assert total==expected_total
 with (out/'样本索引.csv').open('w',encoding='utf-8-sig',newline='') as h:
  w=csv.writer(h);w.writerow(['sample_id','sample_file','matched_pairs']);w.writerows(index)
 if wb:
  for group,ws in worksheets.items():ws.auto_filter.ref=f'A1:W{counts[group]+1}'
  ws=wb.create_sheet('样本索引');ws.append(['sample_id','sample_file','matched_pairs'])
  for r in index:ws.append(r)
  wb.save(out/'DON_最新_截图格式23列转化表.xlsx')
 # Re-read all per-sample CSVs and verify exported tolerances, exact schema and source annotations.
 seen=0
 for f in (out/'逐样本CSV').glob('*.csv'):
  with f.open(encoding='utf-8-sig',newline='') as h:
   r=csv.DictReader(h);assert r.fieldnames==FIELDS
   for row in r:
    assert len(row)==23
    delta=float(row['Dist']);assert math.isclose(delta,float(row['peak.x'])-float(row['peak.y']),rel_tol=0,abs_tol=1e-9)
    assert math.isclose(float(row['Dist.plus'])-delta,1e-5,rel_tol=0,abs_tol=1e-9)
    assert math.isclose(delta-float(row['Dist.minus']),1e-5,rel_tol=0,abs_tol=1e-9);seen+=1
 assert seen==total
 if wb:
  x=load_workbook(out/'DON_最新_截图格式23列转化表.xlsx',read_only=True,data_only=True)
  excelcount=0
  for group in counts:
   it=iter(x[group].values);assert list(next(it))==FIELDS;excelcount+=sum(1 for _ in it)
  assert excelcount==total;x.close()
 report={'version':version,'samples':len(index),'endmembers':len(counts),'pairs':total,'columns':FIELDS,'checks':'PASS: source row IDs/masses/formulas mapped; all CSV rows, difference bounds, counts; latest Excel schema and all rows checked'}
 (out/'核对记录.json').write_text(json.dumps(report,ensure_ascii=False,indent=2));summaries.append(report)
 print('VERIFIED',version,total,flush=True)
(B/'总核对记录.json').write_text(json.dumps(summaries,ensure_ascii=False,indent=2))
