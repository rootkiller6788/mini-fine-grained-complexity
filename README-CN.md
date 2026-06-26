# Mini Fine-Grained Complexity（迷你细粒度复杂度）

**从零开始、零依赖的 C 语言实现**，涵盖细粒度复杂度理论——研究基本计算问题的精确多项式时间下界。每个模块对应 MIT 与 Stanford 的高级计算复杂度课程，将条件性下界、等价类与难度猜想（SETH、ETH、OVC、3SUM、APSP）转化为可运行的 C 代码。

## 子模块

| 子模块 | 主题 | 参考课程 |
|-----------|--------|-------------|
| [mini-3sum-apsp-conjectures](mini-3sum-apsp-conjectures/) | 3SUM 猜想、APSP 猜想、亚立方等价性、极小加乘积、归约链 | MIT 6.8410, Stanford CS367 |
| [mini-conditional-lower-bounds](mini-conditional-lower-bounds/) | 条件性下界框架、SETH/ETH/OVC/3SUM/APSP 假设、归约网络、字符串难度 | MIT 6.8410, Stanford CS254 |
| [mini-edit-distance-lcs-hardness](mini-edit-distance-lcs-hardness/) | 编辑距离、最长公共子序列、Fréchet 距离、序列比对、SETH 下的平方下界 | MIT 6.8410, Stanford CS367 |
| [mini-equivalence-classes](mini-equivalence-classes/) | 亚立方等价类、亚平方等价类、3SUM 等价类、细粒度归约 | MIT 6.8410, Stanford CS367 |
| [mini-k-clique-hardness](mini-k-clique-hardness/) | k-团参数化难度、ETH/SETH 下界、FPT 算法、颜色编码 | MIT 6.8410, Stanford CS367 |
| [mini-orthogonal-vectors](mini-orthogonal-vectors/) | 正交向量问题、SETH→OV 归约、OV→编辑距离/LCS/直径归约、模式匹配 | MIT 6.8410, Stanford CS367 |
| [mini-polynomial-method-consequences](mini-polynomial-method-consequences/) | 多项式方法、代数复杂度、通过代数技术超越暴力搜索 | MIT 6.8420, Stanford CS254 |
| [mini-seth-strong-eth](mini-seth-strong-eth/) | SETH/ETH 形式化、指数级复杂度、CNF 生成、稀疏化 | MIT 6.8410, Stanford CS254 |

## 设计理念

- **零外部依赖** — 纯 C（C99/C11），仅使用 `libc` 和 `libm`
- **模块自包含** — 每个目录自带 `include/`、`src/`、`tests/` 和 `Makefile`
- **理论到代码的映射** — 每个模块包含 `docs/` 目录，内有猜想形式化说明和课程对齐注释
- **实用基准测试** — 经验时间复杂度测量、缩放分析与归约验证

## 构建方式

每个模块相互独立。进入模块目录后运行：

```bash
cd mini-3sum-apsp-conjectures
make all    # 构建全部
make test   # 运行测试
```

需要 **GCC** 和 **GNU Make**。

## 项目结构

```
mini-fine-grained-complexity/
├── mini-3sum-apsp-conjectures/         # 3SUM 与 APSP 猜想、亚立方等价性
├── mini-conditional-lower-bounds/      # 条件性下界框架、归约网络
├── mini-edit-distance-lcs-hardness/    # 编辑距离、LCS、Fréchet 距离的难度
├── mini-equivalence-classes/           # 亚立方、亚平方、3SUM 等价类
├── mini-k-clique-hardness/             # k-团参数化复杂度、ETH/SETH 下界
├── mini-orthogonal-vectors/            # 正交向量问题、SETH 联系、平方归约
├── mini-polynomial-method-consequences/# 细粒度复杂度中的多项式方法
└── mini-seth-strong-eth/               # SETH 与 ETH 形式化、指数级复杂度
```

## 许可证

MIT
