# 截图/PPT第5页格式：23列分子转化明细

## 下载

- **[最新23列Excel（优先查看）](https://github.com/wd475prgfg-cmyk/DON/raw/refs/heads/main/DON_results/screenshot_format/DON_latest_23_columns.xlsx)**：已筛N>0、H/C<1.5，205,000个候选峰对，按10端元分表，sample列区分60个样本，另附样本索引。
- [最新完整包：Excel、60样本CSV、10端元CSV、全部样本CSV及说明](https://github.com/wd475prgfg-cmyk/DON/raw/refs/heads/main/DON_results/screenshot_format/DON_latest_23_columns_complete.zip)
- [旧版23列表：未加N/H/C条件的60样本CSV及说明](https://github.com/wd475prgfg-cmyk/DON/raw/refs/heads/main/DON_results/screenshot_format/DON_previous_23_columns_60_samples.zip)：1,215,666个峰对；按样本分别提供，以免超过Excel单表行数上限。
- [详细方法与PPT对应说明](详细方法与PPT对应说明.md)

## 与截图一致的列顺序

`sample, peak.x, C_x, H_x, O_x, N_x, P_x, S_x, formula_x, class_x, peak.y, C_y, H_y, O_y, N_y, P_y, S_y, formula_y, class_y, Dist, Dist.plus, Dist.minus, Trans.name`

两端元素数、分子式及类别已按样本和原始行号关联并核对；质量差上下限为±0.000010 Da。只是补齐明细格式，没有改变转化网络或RI结果。

**类别说明：** 原始DON表class是N1O3等元素类别，输出原样保留。截图中的LP/LG/AP是另一分类体系，原始数据和PPT没有给出其分类边界，本次没有擅自转换标签。x/y质量排序不证明真实反应方向。

PPT第5页对应本次转化明细；第6页为逐分子次数与RI（≤1，另保留=0），然后才筛端元特异marker。随机森林、XGBoost/SHAP、贝叶斯定量以及组装分析为后续步骤，本次未自动开始。
