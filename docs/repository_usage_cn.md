# 仓库使用说明（中文）

本仓库是《ICI》相关章节的在线计算资源包，用于替代书中冗长的程序附录。建议在书中仅保留简短的“计算资源与在线代码库”说明，并将完整程序、SOP、数据模板和示例放在本仓库中维护。

## 推荐书内说明

```latex
\appendix

\chapter{计算资源与在线代码库}

本书所涉及的完整计算附录作为开放在线资源维护。该资源包括 ICI 框架的 Python 实现、参数提取标准操作程序、不确定性传播程序、跨域验证工具、可视化脚本、完整示例以及配套参数数据库。

为保持纸质书和电子书正文的可读性，本书仅保留核心公式、方法说明和资源索引。完整代码、数据模板、示例程序和更新文档可在以下地址获取：

\[
\texttt{https://github.com/EQTResearcher/ICI-computational-resources}
\]

读者可通过在线代码库获取可执行程序、更新数据、复现实验说明和版本记录。
```

## 上传到 GitHub 的步骤

1. 在 GitHub 账号 `EQTResearcher` 下新建仓库：`ICI-computational-resources`。
2. 将本文件包解压后的所有文件上传到仓库根目录。
3. 在仓库设置中启用 GitHub Pages（可选）：Settings → Pages → Deploy from branch → main / root。
4. 发布第一个 release：`v0.1.0`。
5. 如需 DOI，可将 GitHub 仓库连接 Zenodo，并对 release 归档。

## 版本建议

- `v0.1.0`：从书稿附录整理出的初始公开资源包。
- `v1.0.0`：随正式出版版本冻结的可引用版本。
- 后续更新以 `CHANGELOG.md` 记录。
