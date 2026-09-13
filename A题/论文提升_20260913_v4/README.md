# 基于守恒有限体积法的药材烘干热湿传递与收缩效应分析

本目录保留论文、AI工具使用说明及电子支撑材料。计算模型、程序、数据与数值结果保持不变。

## 文件入口

- `完整论文.pdf`、`完整论文-LaTeX/`：论文及完整源码。
- `AI工具使用详情.pdf`、`AI工具说明-LaTeX/`：AI工具使用说明及源码；当前说明仍标注为修改稿，提交前须核对。
- `solve.py`、`run_all.py`、`reproduce.py`及其他脚本、`utils/`、`requirements.txt`：完整程序和依赖。
- `results/`：全部结果、过程轨迹及既有数值验证记录。
- `figures/`：原有数据图；`新增示意图/`：三张示意图的YAML、draw.io、SVG和PDF源文件。
- `A题.pdf`、`附件/`、`assets/`：原题、原始附件及字体和许可。
- `题目分析报告.md`、`术语表格.md`：保留的前期技术材料，当前结构以论文为准。

正文一级标题为：一、问题重述；二、问题分析；三、符号说明；四、数据处理和模型准备；五、模型的建立与求解；六、模型评价与推广；七、结论。模型假设位于4.1，四问位于5.1—5.4，综合检验位于5.5。

## 主要结果

问题三完成时间57.4728 h，问题四51.0908 h，问题四终点半径1.2000 cm；总体时间缩短11.1043%，固定问题四物性时收缩降幅60.6524%。整秒输出不等同于亚秒绝对精度。

`result2.xlsx`两张工作表均保留1—206902 s的全过程逐秒结果；问题三、四每60 s输出并补入各自终点。所有完整代码保留在电子支撑材料和论文附录B中。

## 运行与编译

在本目录准备Python依赖（版本见`requirements.txt`）后运行：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B reproduce.py
```

脚本建立独立复现输出，不覆盖正式结果。进入`完整论文-LaTeX/`或`AI工具说明-LaTeX/`，使用XeLaTeX编译两遍；需要ctex、fvextra、常用数学表格宏包及Noto CJK字体：

```bash
xelatex -no-shell-escape -interaction=nonstopmode -halt-on-error main.tex
xelatex -no-shell-escape -interaction=nonstopmode -halt-on-error main.tex
```

本地依赖环境、编辑过程、编译输出、旧版Word稿和重复预览已移出仓库。运行与编译前需自行安装依赖。Git历史和推送配置保留，未执行提交或推送。

AI使用声明完整保留，使用者须独立审阅论文并确认披露与提交要求。
