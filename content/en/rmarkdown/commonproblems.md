---
title: "Common Problems"
linkTitle: "Common Problems"
weight: 6
type: docs
---

### 1. During site building: 

`Error: Error building site: POSTCSS: failed to transform "scss/main.css" (text/css):resource "scss/scss/main.scss_4853eb546e7a6c0898ed71feae7357c0" not found in file cache`.

```bash
cd YOUR_REPO
npm audit fix
```
More [information](https://github.com/google/docsy/issues/235)

### 2. During site building: 

`Error: Error building site: "~/content/en/_index.html:6:1": failed to extract shortcode: template for shortcode "blocks/cover" not found`.

```bash
cd YOUR_REPO
git submodule update --init --recursive
```

### 3. Rmarkdown rendering

```r
The extension tex_math_dollars is not supported for gfm
Error: pandoc document conversion failed with error 23
Execution halted
Error : Failed to render content/en/docs/test.Rmarkdown
```
Update pandoc >= 2.12 solves the problem. Use `rmarkdown::pandoc_version()` to check version.

If you use Rstudio, simply update Rstudio to the latest version which comes with pandoc >= 2.12 and 
rmarkdown will use it automatically. Otherwise, you need to [update pandoc](https://pandoc.org/installing.html) by yourself.
