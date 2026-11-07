from pathlib import Path
import sys,json,re,math,itertools
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Patch
from scipy.stats import skew,kurtosis
import openpyxl
from openpyxl.styles import Font,PatternFill
from openpyxl.utils import get_column_letter
B=Path(__file__).resolve().parent;read=lambda p:pd.read_csv(p,float_precision='round_trip')
PROCESS=['variable_selection','homogeneous_selection','homogenizing_dispersal','dispersal_limitation','drift'];PCN=['异质选择','同质选择','均质化扩散','扩散限制','漂变与未主导过程'];colors=['#D55E00','#E69F00','#56B4E9','#009E73','#8DA0CB']
summ=read(B/'两套marker_保守性验证汇总.csv');proc=read(B/'两套marker_组装过程占比.csv');checks={};boundary=[]
def classify(z,rc):
 if z>2:return PROCESS[0]
 if z< -2:return PROCESS[1]
 if rc< -.95:return PROCESS[2]
 if rc>.95:return PROCESS[3]
 return PROCESS[4]
# Check weighted sampling law against an enumerated exact small example.
rng=np.random.default_rng(873261);w=np.array([1.,2.,3.,4.]);draws=200000;keys=rng.exponential(size=(draws,4))/w;chosen=np.sort(np.argpartition(keys,1,axis=1)[:,:2],axis=1)
maxerror=0
for i,j in itertools.combinations(range(4),2):
 p=w[i]/w.sum()*w[j]/(w.sum()-w[i])+w[j]/w.sum()*w[i]/(w.sum()-w[j]);phat=np.mean((chosen[:,0]==i)&(chosen[:,1]==j));maxerror=max(maxerror,abs(p-phat));assert abs(p-phat)<.004
checks['weighted_without_replacement_vs_exact_R_successive_sampling']={'draws':draws,'max_absolute_frequency_error':maxerror,'status':'PASS'}
for setname in ['RI_le1','RI_eq0']:
 out=B/setname;presence=read(out/'Lake_实际检出marker矩阵.csv').set_index('sample_id');ids=presence.index.tolist();C=presence.to_numpy(dtype=bool);pairs=list(itertools.combinations(range(10),2));rca=np.load(out/'RCBray_9999次随机交集.npz');null_shared=rca['shared'];obshared=rca['observed_shared'];assert null_shared.shape==(9999,45)
 rc=[]
 for k,(i,j) in enumerate(pairs):
  observed=int(np.logical_and(C[i],C[j]).sum());assert observed==obshared[k]
  assert ((null_shared[:,k]>=0)&(null_shared[:,k]<=min(C[i].sum(),C[j].sum()))).all()
  y=(null_shared[:,k]>observed).astype(float)+.5*(null_shared[:,k]==observed);rc.append(2*y.mean()-1)
 rc=np.array(rc)
 for mode,file,filename in [('primary','主分析_45样本对结果.csv','betaMNTD_主分析_9999次.npz'),('sensitivity','完整来源库敏感性_45样本对结果.csv','betaMNTD_完整来源库敏感性_9999次.npz')]:
  p=read(out/file);null=np.load(out/filename)['betaMNTD'];assert null.shape==(9999,45) and np.isfinite(null).all()
  np.testing.assert_allclose(p.null_mean_betaMNTD,null.mean(axis=0),rtol=1e-12,atol=1e-12);np.testing.assert_allclose(p.null_SD_betaMNTD,null.std(axis=0,ddof=1),rtol=1e-12,atol=1e-12)
  z=(p.observed_betaMNTD.to_numpy()-null.mean(axis=0))/null.std(axis=0,ddof=1);np.testing.assert_allclose(p.betaNTI,z,rtol=1e-12,atol=1e-12);np.testing.assert_allclose(p.RCBray,rc,rtol=1e-12,atol=1e-12)
  assert list(p.process)==[classify(a,b) for a,b in zip(z,rc)]
  zse=np.sqrt(np.maximum(0,1+z*skew(null,axis=0,bias=False)+z*z/4*(kurtosis(null,axis=0,fisher=False,bias=False)-1))/9999)
  for k,(i,j) in enumerate(pairs):
   y=(null_shared[:,k]>obshared[k]).astype(float)+.5*(null_shared[:,k]==obshared[k]);rcse=2*y.std(ddof=1)/np.sqrt(9999)
   zlow,zhi=z[k]-1.96*zse[k],z[k]+1.96*zse[k];rclow,rchi=max(-1,rc[k]-1.96*rcse),min(1,rc[k]+1.96*rcse)
   boundary.append({'marker_set':setname,'analysis':mode,'sample_i':ids[i],'sample_j':ids[j],'betaNTI':z[k],'betaNTI_MC_SE_delta_method':zse[k],'betaNTI_MC_lower95':zlow,'betaNTI_MC_upper95':zhi,'betaNTI_MC_interval_crosses_threshold':int((zlow<=-2<=zhi)or(zlow<=2<=zhi)),'RCBray':rc[k],'RCBray_MC_SE':rcse,'RCBray_MC_lower95':rclow,'RCBray_MC_upper95':rchi,'RCBray_MC_interval_crosses_threshold':int((rclow<=-.95<=rchi)or(rclow<=.95<=rchi)),'process_at_9999':p.iloc[k].process})
  checks[setname+'_'+mode]='PASS: all 45 saved null distributions, betaNTI, RCBray, process classifications'
 # Newick independently parsed to verify patristic distances, including unary-node pruning.
 def parse_newick(text):
  tok=re.findall(r'[(),:;]|[^(),:;\s]+',text);pos=0;serial=0;parent={};edge={};leaves={}
  def node():
   nonlocal pos,serial
   serial+=1;me=serial
   if tok[pos]=='(':
    pos+=1
    while True:
     child,length=node();parent[child]=me;edge[child]=length
     if tok[pos]==',':pos+=1;continue
     assert tok[pos]==')';pos+=1;break
   else:leaves[tok[pos]]=me;pos+=1
   length=0.
   if tok[pos]==':':pos+=1;length=float(tok[pos]);pos+=1
   return me,length
  root,_=node();assert tok[pos]==';'
  def distance(a,b):
   anc={};x=leaves[a];d=0
   while True:
    anc[x]=d
    if x not in parent:break
    d+=edge[x];x=parent[x]
   x=leaves[b];d=0
   while x not in anc:d+=edge[x];x=parent[x]
   return d+anc[x]
  return leaves,distance
 for treefile in ['全部候选marker_MCD.nwk','Lake检出marker_修剪MCD.nwk']:
  leaves,distance=parse_newick((out/treefile).read_text());saved=np.load(out/'完整marker树距离与linkage.npz');lookup={x:i for i,x in enumerate(saved['marker_ids'])};dm=saved['distance'];labs=list(leaves)
  rng=np.random.default_rng(183);tests=[tuple(rng.choice(labs,2,replace=False)) for _ in range(150)]
  for a,b in tests:assert math.isclose(distance(a,b),dm[lookup[a],lookup[b]],rel_tol=1e-11,abs_tol=1e-11)
 checks[setname+'_Newick_export']='PASS: independently parsed full and pruned trees, 150 pairs each'
boundary=pd.DataFrame(boundary);boundary.to_csv(B/'蒙特卡洛阈值敏感性.csv',index=False,encoding='utf-8-sig')
# Plot with vector text, and declare the user's required fonts in the delivered SVG.
# Fonts are not embedded. Preview PNGs use available serif fonts and remain outside result exports.
plt.rcParams.update({'font.family':'DejaVu Serif','svg.fonttype':'none','axes.unicode_minus':False,'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
zh=FontProperties(family='Noto Serif CJK SC',fname='/usr/share/fonts/opentype/noto/NotoSerifCJK-Regular.ttc');eng=FontProperties(family='DejaVu Serif')
def save_svg(fig,path,preview):
 fig.savefig(path,format='svg',bbox_inches='tight')
 s=path.read_text();s=s.replace('DejaVu Serif','Times New Roman')
 for family in ['Noto Serif CJK JP','Noto Serif CJK SC','Noto Serif CJK TC','Noto Serif CJK HK','Noto Serif CJK KR']:s=s.replace(family,'SimSun')
 assert 'Times New Roman' in s and 'SimSun' in s
 path.write_text(s);fig.savefig(preview,dpi=150,bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(2,2,figsize=(12,9),gridspec_kw={'width_ratios':[1,1.1]});rng=np.random.default_rng(724)
for row,(mode,filename,label) in enumerate([('primary_observed_Lake_pool','主分析_45样本对结果.csv','主分析'),('sensitivity_full_source_pool','完整来源库敏感性_45样本对结果.csv','完整来源分子库敏感性分析')]):
 vals=[read(B/sn/filename).betaNTI.to_numpy() for sn in ['RI_le1','RI_eq0']];a=ax[row,0]
 a.axhspan(-2,2,color='#009E73',alpha=.08);a.axhline(2,c='black',lw=.8,ls='--');a.axhline(-2,c='black',lw=.8,ls='--')
 bp=a.boxplot(vals,positions=[1,2],widths=.35,patch_artist=True,showfliers=False,medianprops={'color':'black'})
 for box,c in zip(bp['boxes'],['#56B4E9','#D55E00']):box.set_facecolor(c);box.set_alpha(.35)
 for pos,v,c in zip([1,2],vals,['#0072B2','#D55E00']):a.scatter(pos+rng.uniform(-.1,.1,len(v)),v,s=13,color=c,alpha=.65,edgecolors='none')
 a.set_xticks([1,2],['RI ≤ 1','RI = 0']);a.set_ylabel('βNTI');a.set_ylim(-4.6,3.9);a.set_xlim(.5,2.5);a.set_title(label+'的最近邻偏离指数',fontproperties=zh,fontsize=12);a.grid(axis='y',alpha=.15)
 b=ax[row,1];bottom=np.zeros(2)
 for process,color in zip(PROCESS,colors):
  data=np.array([proc.loc[(proc.marker_set==sn)&(proc.analysis==mode)&(proc.process==process),'percent_of_all_pairs'].iloc[0] for sn in ['RI_le1','RI_eq0']]);b.bar([0,1],data,bottom=bottom,width=.52,color=color,edgecolor='white',linewidth=.4)
  for x,value,bot in zip([0,1],data,bottom):
   if value>=5:b.text(x,bot+value/2,f'{value:.1f}%',ha='center',va='center',fontsize=10)
  bottom+=data
 for x,sn in enumerate(['RI_le1','RI_eq0']):
  value=summ.loc[(summ.marker_set==sn)&(summ.analysis==mode),'stochastic_percent_of_valid'].iloc[0];b.text(x,103,f'Stochastic: {value:.2f}%',ha='center',fontsize=10)
 b.set_xticks([0,1],['RI ≤ 1','RI = 0']);b.set_ylabel('Sample-pair proportion (%)');b.set_ylim(0,112);b.set_yticks([0,25,50,75,100]);b.set_title(label+'的组装过程比例',fontproperties=zh,fontsize=12)
fig.suptitle('两套分子标志物的湖泊组装过程验证',fontproperties=zh,fontsize=16,y=.99)
fig.legend(handles=[Patch(facecolor=c,label=t) for c,t in zip(colors,PCN)],loc='lower center',bbox_to_anchor=(.5,.035),ncol=5,prop=zh,frameon=False)
fig.text(.5,.014,'9,999 randomizations; 10 Lake samples; 45 non-independent sample pairs.',ha='center',fontproperties=eng,fontsize=10)
fig.tight_layout(rect=[0,.09,1,.96]);save_svg(fig,B/'两套marker_保守性验证对比.svg',Path('/workspace/don-setup/conservatism_preview.png'))
fig,axes=plt.subplots(2,2,figsize=(12,10))
for row,sn in enumerate(['RI_le1','RI_eq0']):
 p=read(B/sn/'主分析_45样本对结果.csv');meta=read(B/'Lake样本身份表.csv');ids=meta.sample_id.tolist();labels=meta['sample'].tolist();lookup={x:i for i,x in enumerate(ids)}
 for col,(field,vmax) in enumerate([('betaNTI',4.),('RCBray',1.)]):
  matrix=np.full((10,10),np.nan)
  for x in p.itertuples():i,j=lookup[x.sample_i],lookup[x.sample_j];matrix[i,j]=matrix[j,i]=getattr(x,field)
  a=axes[row,col];im=a.imshow(matrix,cmap='RdBu_r',vmin=-vmax,vmax=vmax);a.set_xticks(range(10),labels,rotation=45,ha='right');a.set_yticks(range(10),labels);title=('RI ≤ 1' if sn=='RI_le1' else 'RI = 0')+': '+('βNTI' if field=='betaNTI' else 'RCBray');a.set_title(title)
  for i in range(10):
   for j in range(i+1,10):a.text(j,i,f'{matrix[i,j]:.3f}',ha='center',va='center',fontsize=6.6,color='white' if abs(matrix[i,j])>.6*vmax else 'black')
  fig.colorbar(im,ax=a,fraction=.046,pad=.04)
fig.suptitle('主分析的样本对指标矩阵',fontproperties=zh,fontsize=16,y=.99);fig.tight_layout(rect=[0,0,1,.96]);save_svg(fig,B/'主分析_betaNTI_RCBray热图.svg',Path('/workspace/don-setup/conservatism_heatmap_preview.png'))
# Format Excel and add methods/uncertainty, without changing numeric results.
xpath=B/'DON_两套marker_湖泊保守性验证.xlsx';w=openpyxl.load_workbook(xpath)
if 'MC阈值敏感性' not in w.sheetnames:
 ws=w.create_sheet('MC阈值敏感性');ws.append(boundary.columns.tolist())
 for row in boundary.itertuples(index=False,name=None):ws.append(row)
if '方法与解释' in w.sheetnames:del w['方法与解释']
notes=w.create_sheet('方法与解释');notes.append(['项目','说明'])
notes_rows=[('主分析分子库','仅在Lake至少一次检出的marker：RI≤1为173种，RI=0为101种；保留原完整marker树的分子特征尺度'),('敏感性分子库','完整来源marker候选库：RI≤1为726种，RI=0为484种；与主分析使用同一完整来源特征树'),('MCD构建','C/H/O/N/S/P、DBE、AI_Mod、kdefect按样本SD做Z标准化，欧氏距离、UPGMA；P为常数0，剔除该零方差列'),('强度处理','不做峰强度归一化；零模型使用检出/未检出。MCD的分子特征标准化不等于强度归一化'),('βMNTD','按picante非丰度加权实现：将两个方向的最近邻距离拼接后取均值，保留共享分子的零距离'),('随机化','各方案βNTI和RCBray均为9999次，固定种子，保存全部零模型分布'),('过程分类','βNTI>2异质选择，<-2同质选择；其余按RC<-0.95、>0.95、其余分别为均质化扩散、扩散限制、漂变/未主导'),('阈值精度','使用未经四舍五入的数值判定；MC阈值敏感性表中的区间只反映有限随机化误差，不是生态学置信区间'),('解释边界','结果为45个样本对的组装过程比例，不是每个marker的保守性分数，也不是通过验证的分子百分比；样本对不独立'),('与WR关系','依据正文2.3.2、3.3及作者公开代码进行Python数值实现；未获得文章Text S2/S3及原始R会话，不能声称逐项复现原文'),('图表字体','SVG以文本保存，指定中文SimSun宋体、英文数字Times New Roman；不嵌入字体，查看电脑需安装相应字体')]
for row in notes_rows:notes.append(row)
engfont=Font(name='Times New Roman',size=11);cnfont=Font(name='SimSun',size=11);head=Font(name='Times New Roman',bold=True,color='FFFFFF');fill=PatternFill('solid',fgColor='305496')
for ws in w:
 ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions
 for cell in ws[1]:
  cell.font=Font(name='SimSun',bold=True,color='FFFFFF') if isinstance(cell.value,str) and re.search('[\u4e00-\u9fff]',cell.value) else head;cell.fill=fill
 for row in ws.iter_rows(min_row=2):
  for cell in row:cell.font=cnfont if isinstance(cell.value,str) and re.search('[\u4e00-\u9fff]',cell.value) else engfont
 for i,cell in enumerate(ws[1],1):ws.column_dimensions[get_column_letter(i)].width=min(42,max(15,len(str(cell.value))+2))
notes.column_dimensions['A'].width=25;notes.column_dimensions['B'].width=110;w.save(xpath)
wx=openpyxl.load_workbook(xpath,read_only=True,data_only=True);assert wx['两套汇总'].max_row==5 and wx['RI_le1_主分析'].max_row==46 and wx['RI_eq0_主分析'].max_row==46;wx.close()
checks['Excel_saved_rows']='PASS summary=4 analyses, 45 pairs per marker/pool';checks['SVG_font_declarations']='Times New Roman and SimSun; text, not embedded fonts';(B/'独立数值与导出核对.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2));print(json.dumps(checks,ensure_ascii=False,indent=2),flush=True)
