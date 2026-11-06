# DON 全部分析结果下载

完整结果按内容拆为三个独立 ZIP，每个可单独解压。下载并解压到同一文件夹即可得到本次全部结果，无需拼接文件。

| 下载 | 内容 |
| --- | --- |
| [分包1：原始数据与10端元合并表](https://github.com/wd475prgfg-cmyk/DON/raw/refs/heads/main/DON_results/DON_results_part01_data_and_merged.zip) | DON1、DON2、合并表、转化代码、总说明及校验清单 |
| [分包2：全部分子转化分析](https://github.com/wd475prgfg-cmyk/DON/raw/refs/heads/main/DON_results/DON_results_part02_transformations.zip) | 60个样本的基础转化峰对、逐分子次数、汇总、方法及核对记录 |
| [分包3：两版RI结果](https://github.com/wd475prgfg-cmyk/DON/raw/refs/heads/main/DON_results/DON_results_part03_RI.zip) | 原版RI及新增N>0、H/C<1.5后的重算网络、RI名单与核对记录 |

## DON 候选 marker（RI＝0，新增）

[下载RI＝0候选marker Excel](https://github.com/wd475prgfg-cmyk/DON/raw/refs/heads/main/DON_results/marker_RIeq0_Ngt0_HClt1p5/DON_RI_eq0_candidate_markers.xlsx) · [完整结果及口径说明](marker_RIeq0_Ngt0_HClt1p5/README.md)

保持N>0、H/C<1.5及9来源排重口径，RI门槛改为0次：484种来源独有候选，101种Lake检出；另列441种全部来源检出记录均为0次的更严格名单。

## DON 候选 marker（RI≤1，对照版）

[下载候选marker Excel](https://github.com/wd475prgfg-cmyk/DON/raw/refs/heads/main/DON_results/marker_Ngt0_HClt1p5/DON_candidate_markers.xlsx) · [完整结果、筛选规则及9来源汇总](marker_Ngt0_HClt1p5/README.md)

最新N>0、H/C<1.5、RI≤1结果中，按9来源间全部检出分子式排重，筛得726种来源独有候选，其中173种在Lake受体检出。Lake不参与来源间排重。

## PPT截图格式的分子转化明细（新增）

[下载最新23列转化Excel](https://github.com/wd475prgfg-cmyk/DON/raw/refs/heads/main/DON_results/screenshot_format/DON_latest_23_columns.xlsx) · [查看60样本明细、旧版表及PPT对应说明](screenshot_format/README.md)

## 最新结果

[直接下载最新RI Excel](https://github.com/wd475prgfg-cmyk/DON/raw/refs/heads/main/DON_results/DON_RI_Ngt0_HClt1p5.xlsx) · [查看10端元汇总](latest_RI_endmember_summary.csv)

最新口径：实测m/z 200–900 Da、强度>0、元素组成有效、去除明确标注同位素峰；新增N>0、H/C<1.5后，逐样本重新计算转化网络，RI为转化次数≤1，另保留=0的严格标记。

- 60个样本、10个端元，筛后47,959条逐样本分子峰记录、205,000对候选转化峰对。
- RI≤1为6,999条逐样本记录，跨样本按分子式去重3,306种。
- RI=0为3,471条逐样本记录。

解压后优先查看 `RI_N大于0_HC小于1.5` 文件夹；旧版结果保留供比较。所有结果已完成核对。现已完成候选marker筛选，详见上方新增marker目录；尚未进行定量贡献建模。
