# DON 来源端元特异候选 marker

**[下载候选marker Excel](https://github.com/wd475prgfg-cmyk/DON/raw/refs/heads/main/DON_results/marker_Ngt0_HClt1p5/DON_candidate_markers.xlsx)** · [下载完整结果、CSV及代码](https://github.com/wd475prgfg-cmyk/DON/raw/refs/heads/main/DON_results/marker_Ngt0_HClt1p5/DON_candidate_marker_results.zip) · [查看详细方法](方法与结果说明.md) · [查看来源汇总](source_summary.csv)

## 筛选口径

- 沿用N>0、H/C<1.5筛后网络的RI≤1名单。
- Lake作为受体，其他9类作为来源（50个来源样本、10个Lake样本）。
- 某分子式在对应来源至少一次为RI，且在其他8类来源的**全部有效检出分子**中均未出现，保留为来源独有候选。
- Lake用于检出验证，不参与来源间排重。Lake检出不要求该分子在Lake中也属于RI。
- 不擅自增加检出率阈值；提供每个候选的来源检出率、重复检出、全部来源记录均RI和严格=0标记。

## 结果

| 来源 | 候选marker | Lake检出 |
| --- | ---: | ---: |
| 藻类 Algae | 22 | 12 |
| 磷矿 Apatit | 27 | 8 |
| 水产养殖尾水 Aquaculture | 83 | 28 |
| 畜禽粪便 Excrement | 75 | 21 |
| 肥料 Fertilizer | 61 | 2 |
| 植物 Plant | 60 | 10 |
| 污水处理厂尾水 STP | 39 | 20 |
| 底泥 Sediment | 176 | 37 |
| 土壤 Soil | 183 | 35 |
| 合计 | **726** | **173** |

4,204个来源×RI候选组合中，3,478个因其他来源检出而排除。主名单726种分子式互不重复；484种至少一次满足严格RI=0，672种全部来源检出记录均为RI≤1，75种在对应来源至少两个样本检出。

这些是当前采样条件下的候选marker；651种仅在单个来源样本检出，特异性尚需验证。Excel保留了可追溯的源行号、检出率、60样本0/1矩阵和Lake原始强度矩阵。173种Lake检出不能直接解释为来源贡献比例。

完整集合差独立核对、全部726种候选的来源/Lake样本数及60样本检出矩阵、Excel21个工作表均已核对通过。
