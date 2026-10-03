# chromatogram-internal-standard

色谱内标定量工具（仅依赖 Python 标准库）。根据标样的校正峰面积比与浓度拟合
直线，反解未知样浓度。

## 用法

```
python -m chromatogram 输入.json --out 输出目录
```

成功时在输出目录生成：

- `results.csv`：逐未知样一行，列为 `id, analyte_area, internal_standard_area,
  area_ratio, concentration, extrapolated`。
- `calibration.json`：校准报告，含 `slope`、`intercept`、`concentration_range`
  及各标样的面积、面积比、拟合值与残差（残差 = 实测 - 拟合）。

输入非法或整批无法定量（如内标面积非正、拟合斜率非正）时，以非零状态码退出，
且不生成任何报告文件。

## 输入格式

```json
{
  "analyte_window": [t0, t1],
  "internal_standard_window": [t0, t1],
  "standards": [
    {"id": "std1", "concentration": 1.0, "curve": [[t, signal], ...]}
  ],
  "unknowns": [
    {"id": "sample1", "curve": [[t, signal], ...]}
  ]
}
```

约束：

- 样品编号在标样与未知样中全局唯一。
- 每条曲线至少两个有限采样点，时间严格递增，且完整覆盖两个时间窗口。
- 两个窗口均为正长度且互不重叠。
- 标样浓度有限且非负，至少三个不同浓度的标样；未知样不得提供浓度。
- 假定所有样品加入相同浓度的内标，因此直接用面积比作纵轴。

## 计算口径

**基线与积分**：曲线在节点间线性变化，窗口边界不在节点时线性插值。每个窗口
以两端信号（插值后）的连线为基线，扣除基线后只积分正值部分。跨零线段在零点
解析切分后再积分，而不是把负节点截零后直接梯形积分。

**校准**：以全部标样浓度为横轴、待测峰面积/内标峰面积比为纵轴，拟合带截距的
普通最小二乘直线；斜率非正则整批失败。任一样本内标面积非正也整批失败。

**定量与外推**：未知样由直线反解浓度 `concentration = (area_ratio -
intercept) / slope`，保留负值不归零。`extrapolated` 表示反解浓度是否超出标样
浓度范围，恰好等于边界不算外推。

## 测试

```
python -m unittest discover -s tests
```
