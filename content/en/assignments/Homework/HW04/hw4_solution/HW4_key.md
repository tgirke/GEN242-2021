---
title: "HW4: Pairwise Alignments" 
author: "Author: Your Name"
date: "Last update: 22 April, 2021" 
output:
  html_document:
    toc: true
    toc_float:
        collapsed: true
        smooth_scroll: true
    toc_depth: 3
    fig_caption: yes
    code_folding: show
    number_sections: false
fontsize: 14pt
bibliography: bibtex.bib
type: docs
weight: 5
---

<script src="/rmarkdown-libs/kePrint/kePrint.js"></script>
<link href="/rmarkdown-libs/lightable/lightable.css" rel="stylesheet" />
<script src="/rmarkdown-libs/kePrint/kePrint.js"></script>

<link href="/rmarkdown-libs/lightable/lightable.css" rel="stylesheet" />

<!---
- Compile from command-line
Rscript -e "rmarkdown::render('HW4_key.Rmd', c('html_document'), clean=FALSE)"
-->

<div style="text-align: right">

Source code downloads:    
\[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/assignments/Homework/HW04/hw4_solution/HW4_key.Rmd) \]    
\[ [pairwiseAlign\_Fct.R](https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/assignments/Homework/HW04/hw4_solution/pairwiseAlign_Fct.R) \]

</div>

# Rendering Instructions

To render this R Markdown document, one needs to download the following files to the same directory.

-   [HW4\_key.Rmd](https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/assignments/Homework/HW04/hw4_solution/HW4_key.Rmd)
-   [pairwiseAlign\_Fct.R](https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/assignments/Homework/HW04/hw4_solution/pairwiseAlign_Fct.R)
-   [bibtex.bib](https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/assignments/Homework/HW04/hw4_solution/bibtex.bib)

After this the HTML report can be rendered from the corresponding directory with `rmarkdown::render()` as follows.

``` r
rmarkdown::render('HW4_key.Rmd') # From R
Rscript -e "rmarkdown::render('HW4_key.Rmd')" # From command-line
```

# A. Choice of Sequence Type

Task 1: Which sequence type - amino acid or nucleotide - is more appropriate to search databases
for remotely related sequences? Provide at least three reasons for your decision.

Answer: When coding sequences are expected to have weak similarities then one
should use protein sequences rather than DNA sequences for database searching,
because of (1) their higher information content (20 versus 4 letter alphabet),
as well as (2) the better scoring and (3) functional classification systems available
for amino acids.

# B. Dynamic Programming for Pairwise Alignments

Task 2: Create manually (or write an R script for it) one global and one local
alignment for the following two protein sequences using the Needleman-Wusch and
Smith-Waterman algorithms, respectively (Smith and Waterman 1981; Needleman and Wunsch 1970).

``` sh
O15528: PFGFGKRSCMGRRLA
P98187: FIPFSAGPRNCIGQK
```

## Source functions

All alignment functions used in the following sections are defined in the downloaded R script file
that is named `pairwiseAlign_Fct.R`. These functions are loaded with the `source()` command below.

``` r
source("pairwiseAlign_Fct.R")
```

## Input sequences

Define within R or import them (here former).

``` r
S1 <- "PFGFGKRSCMGRRLA"
S2 <- "FIPFSAGPRNCIGQK"
```

Additional test sequences

``` r
# S1 <- "HEAGAWGHEE"
# S2 <- "PAWHEAE"
```

## Global alignment

The alignment type choice is passed on to all following functions.

``` r
align_type <- "global"
# align_type <- "local"
```

### Dynamic programming matrices

``` r
dynMA <- dynProgMatrix(S1, S2, align_method=align_type, gap_penalty=8, substitutionMA="BLOSUM50")
```

The matrices are stored in a list and returned below.

``` r
dynMA
```

    ## $glob_ma
    ##      gp    P   F   G   F   G   K   R   S   C   M   G   R    R    L    A
    ## gp    0   -8 -16 -24 -32 -40 -48 -56 -64 -72 -80 -88 -96 -104 -112 -120
    ## F    -8   -4   0  -8 -16 -24 -32 -40 -48 -56 -64 -72 -80  -88  -96 -104
    ## I   -16  -11  -4  -4  -8 -16 -24 -32 -40 -48 -54 -62 -70  -78  -86  -94
    ## P   -24   -6 -12  -6  -8 -10 -17 -25 -33 -41 -49 -56 -64  -72  -80  -87
    ## F   -32  -14   2  -6   2  -6 -14 -20 -28 -35 -41 -49 -57  -65  -71  -79
    ## S   -40  -22  -6   2  -6   2  -6 -14 -15 -23 -31 -39 -47  -55  -63  -70
    ## A   -48  -30 -14  -6  -1  -6   1  -7 -13 -16 -24 -31 -39  -47  -55  -58
    ## G   -56  -38 -22  -6  -9   7  -1  -2  -7 -15 -19 -16 -24  -32  -40  -48
    ## P   -64  -46 -30 -14 -10  -1   6  -2  -3 -11 -18 -21 -19  -27  -35  -41
    ## R   -72  -54 -38 -22 -17  -9   2  13   5  -3 -11 -19 -14  -12  -20  -28
    ## N   -80  -62 -46 -30 -25 -17  -6   5  14   6  -2 -10 -18  -15  -16  -21
    ## C   -88  -70 -54 -38 -32 -25 -14  -3   6  27  19  11   3   -5  -13  -17
    ## I   -96  -78 -62 -46 -38 -33 -22 -11  -2  19  29  21  13    5   -3  -11
    ## G  -104  -86 -70 -54 -46 -30 -30 -19 -10  11  21  37  29   21   13    5
    ## Q  -112  -94 -78 -62 -54 -38 -28 -27 -18   3  13  29  38   30   22   14
    ## K  -120 -102 -86 -70 -62 -46 -32 -25 -26  -5   5  21  32   41   33   25
    ## 
    ## $glob_ma_path
    ##    gp P F G F G K R S C M G R R L A
    ## gp  0 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
    ## F   2 1 1 3 1 3 3 3 3 3 3 3 3 3 3 3
    ## I   2 1 1 1 1 3 3 3 3 3 1 3 3 3 1 3
    ## P   2 1 2 1 1 1 1 3 1 3 3 1 3 3 3 1
    ## F   2 2 1 3 1 3 1 1 1 1 1 3 3 3 1 3
    ## S   2 2 2 1 2 1 1 3 1 3 3 3 3 3 3 1
    ## A   2 2 2 1 1 1 1 3 1 1 1 1 3 3 3 1
    ## G   2 2 2 1 2 1 3 1 1 3 1 1 3 3 3 3
    ## P   2 1 2 2 1 2 1 3 1 1 1 1 1 1 3 1
    ## R   2 2 2 2 1 2 1 1 3 3 3 3 1 1 3 3
    ## N   2 2 2 2 2 1 2 2 1 3 3 3 3 1 1 1
    ## C   2 2 2 2 1 2 2 2 2 1 3 3 3 3 3 1
    ## I   2 2 2 2 1 2 2 2 2 2 1 3 3 3 1 3
    ## G   2 2 2 1 2 1 2 2 2 2 2 1 3 3 3 3
    ## Q   2 2 2 2 2 2 1 2 2 2 2 2 1 1 3 3
    ## K   2 2 2 2 2 2 1 1 2 2 2 2 1 1 3 3

### Compute alignment

The following `alignList` stores all relevant results in a list, including dynamic
programming matrices, as well as the coordinates (named `path_coor`) to highlight path in
dynamic progamming matrix (see below).

``` r
alignList <- alignmentTraceback(ma=dynMA[[1]], ma_path=dynMA[[2]], align_method=align_type) 
names(alignList)
```

    ## [1] "ma"        "ma_path"   "path_coor" "as1"       "consensus" "as2"      
    ## [7] "score"

``` r
# alignList$ma # dyn ma with scores
# alignList$ma_path # dyn ma with path
# alignList$path_coor # coordinates for path to auto highlight path in HTML/PDF table
```

### Return results

#### Traceback in matrix

The following prints the fully populated dynamic programming matrix where the traceback
path is highlighted in color.

``` r
printColMa(alignList)
```

<table class=" lightable-paper lightable-striped" style="font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:left;">
gp
</th>
<th style="text-align:left;">
P
</th>
<th style="text-align:left;">
F
</th>
<th style="text-align:left;">
G
</th>
<th style="text-align:left;">
F
</th>
<th style="text-align:left;">
G
</th>
<th style="text-align:left;">
K
</th>
<th style="text-align:left;">
R
</th>
<th style="text-align:left;">
S
</th>
<th style="text-align:left;">
C
</th>
<th style="text-align:left;">
M
</th>
<th style="text-align:left;">
G
</th>
<th style="text-align:left;">
R
</th>
<th style="text-align:left;">
R
</th>
<th style="text-align:left;">
L
</th>
<th style="text-align:left;">
A
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
gp
</td>
<td style="text-align:left;">
<span style="     color: red !important;">0</span>
</td>
<td style="text-align:left;">
-8
</td>
<td style="text-align:left;">
-16
</td>
<td style="text-align:left;">
-24
</td>
<td style="text-align:left;">
-32
</td>
<td style="text-align:left;">
-40
</td>
<td style="text-align:left;">
-48
</td>
<td style="text-align:left;">
-56
</td>
<td style="text-align:left;">
-64
</td>
<td style="text-align:left;">
-72
</td>
<td style="text-align:left;">
-80
</td>
<td style="text-align:left;">
-88
</td>
<td style="text-align:left;">
-96
</td>
<td style="text-align:left;">
-104
</td>
<td style="text-align:left;">
-112
</td>
<td style="text-align:left;">
-120
</td>
</tr>
<tr>
<td style="text-align:left;">
F
</td>
<td style="text-align:left;">
<span style="     color: red !important;">-8</span>
</td>
<td style="text-align:left;">
-4
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
-8
</td>
<td style="text-align:left;">
-16
</td>
<td style="text-align:left;">
-24
</td>
<td style="text-align:left;">
-32
</td>
<td style="text-align:left;">
-40
</td>
<td style="text-align:left;">
-48
</td>
<td style="text-align:left;">
-56
</td>
<td style="text-align:left;">
-64
</td>
<td style="text-align:left;">
-72
</td>
<td style="text-align:left;">
-80
</td>
<td style="text-align:left;">
-88
</td>
<td style="text-align:left;">
-96
</td>
<td style="text-align:left;">
-104
</td>
</tr>
<tr>
<td style="text-align:left;">
I
</td>
<td style="text-align:left;">
<span style="     color: red !important;">-16</span>
</td>
<td style="text-align:left;">
-11
</td>
<td style="text-align:left;">
-4
</td>
<td style="text-align:left;">
-4
</td>
<td style="text-align:left;">
-8
</td>
<td style="text-align:left;">
-16
</td>
<td style="text-align:left;">
-24
</td>
<td style="text-align:left;">
-32
</td>
<td style="text-align:left;">
-40
</td>
<td style="text-align:left;">
-48
</td>
<td style="text-align:left;">
-54
</td>
<td style="text-align:left;">
-62
</td>
<td style="text-align:left;">
-70
</td>
<td style="text-align:left;">
-78
</td>
<td style="text-align:left;">
-86
</td>
<td style="text-align:left;">
-94
</td>
</tr>
<tr>
<td style="text-align:left;">
P
</td>
<td style="text-align:left;">
-24
</td>
<td style="text-align:left;">
<span style="     color: red !important;">-6</span>
</td>
<td style="text-align:left;">
-12
</td>
<td style="text-align:left;">
-6
</td>
<td style="text-align:left;">
-8
</td>
<td style="text-align:left;">
-10
</td>
<td style="text-align:left;">
-17
</td>
<td style="text-align:left;">
-25
</td>
<td style="text-align:left;">
-33
</td>
<td style="text-align:left;">
-41
</td>
<td style="text-align:left;">
-49
</td>
<td style="text-align:left;">
-56
</td>
<td style="text-align:left;">
-64
</td>
<td style="text-align:left;">
-72
</td>
<td style="text-align:left;">
-80
</td>
<td style="text-align:left;">
-87
</td>
</tr>
<tr>
<td style="text-align:left;">
F
</td>
<td style="text-align:left;">
-32
</td>
<td style="text-align:left;">
-14
</td>
<td style="text-align:left;">
<span style="     color: red !important;">2</span>
</td>
<td style="text-align:left;">
-6
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
-6
</td>
<td style="text-align:left;">
-14
</td>
<td style="text-align:left;">
-20
</td>
<td style="text-align:left;">
-28
</td>
<td style="text-align:left;">
-35
</td>
<td style="text-align:left;">
-41
</td>
<td style="text-align:left;">
-49
</td>
<td style="text-align:left;">
-57
</td>
<td style="text-align:left;">
-65
</td>
<td style="text-align:left;">
-71
</td>
<td style="text-align:left;">
-79
</td>
</tr>
<tr>
<td style="text-align:left;">
S
</td>
<td style="text-align:left;">
-40
</td>
<td style="text-align:left;">
-22
</td>
<td style="text-align:left;">
-6
</td>
<td style="text-align:left;">
<span style="     color: red !important;">2</span>
</td>
<td style="text-align:left;">
-6
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
-6
</td>
<td style="text-align:left;">
-14
</td>
<td style="text-align:left;">
-15
</td>
<td style="text-align:left;">
-23
</td>
<td style="text-align:left;">
-31
</td>
<td style="text-align:left;">
-39
</td>
<td style="text-align:left;">
-47
</td>
<td style="text-align:left;">
-55
</td>
<td style="text-align:left;">
-63
</td>
<td style="text-align:left;">
-70
</td>
</tr>
<tr>
<td style="text-align:left;">
A
</td>
<td style="text-align:left;">
-48
</td>
<td style="text-align:left;">
-30
</td>
<td style="text-align:left;">
-14
</td>
<td style="text-align:left;">
-6
</td>
<td style="text-align:left;">
<span style="     color: red !important;">-1</span>
</td>
<td style="text-align:left;">
-6
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
-7
</td>
<td style="text-align:left;">
-13
</td>
<td style="text-align:left;">
-16
</td>
<td style="text-align:left;">
-24
</td>
<td style="text-align:left;">
-31
</td>
<td style="text-align:left;">
-39
</td>
<td style="text-align:left;">
-47
</td>
<td style="text-align:left;">
-55
</td>
<td style="text-align:left;">
-58
</td>
</tr>
<tr>
<td style="text-align:left;">
G
</td>
<td style="text-align:left;">
-56
</td>
<td style="text-align:left;">
-38
</td>
<td style="text-align:left;">
-22
</td>
<td style="text-align:left;">
-6
</td>
<td style="text-align:left;">
-9
</td>
<td style="text-align:left;">
<span style="     color: red !important;">7</span>
</td>
<td style="text-align:left;">
-1
</td>
<td style="text-align:left;">
-2
</td>
<td style="text-align:left;">
-7
</td>
<td style="text-align:left;">
-15
</td>
<td style="text-align:left;">
-19
</td>
<td style="text-align:left;">
-16
</td>
<td style="text-align:left;">
-24
</td>
<td style="text-align:left;">
-32
</td>
<td style="text-align:left;">
-40
</td>
<td style="text-align:left;">
-48
</td>
</tr>
<tr>
<td style="text-align:left;">
P
</td>
<td style="text-align:left;">
-64
</td>
<td style="text-align:left;">
-46
</td>
<td style="text-align:left;">
-30
</td>
<td style="text-align:left;">
-14
</td>
<td style="text-align:left;">
-10
</td>
<td style="text-align:left;">
-1
</td>
<td style="text-align:left;">
<span style="     color: red !important;">6</span>
</td>
<td style="text-align:left;">
-2
</td>
<td style="text-align:left;">
-3
</td>
<td style="text-align:left;">
-11
</td>
<td style="text-align:left;">
-18
</td>
<td style="text-align:left;">
-21
</td>
<td style="text-align:left;">
-19
</td>
<td style="text-align:left;">
-27
</td>
<td style="text-align:left;">
-35
</td>
<td style="text-align:left;">
-41
</td>
</tr>
<tr>
<td style="text-align:left;">
R
</td>
<td style="text-align:left;">
-72
</td>
<td style="text-align:left;">
-54
</td>
<td style="text-align:left;">
-38
</td>
<td style="text-align:left;">
-22
</td>
<td style="text-align:left;">
-17
</td>
<td style="text-align:left;">
-9
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
<span style="     color: red !important;">13</span>
</td>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
-3
</td>
<td style="text-align:left;">
-11
</td>
<td style="text-align:left;">
-19
</td>
<td style="text-align:left;">
-14
</td>
<td style="text-align:left;">
-12
</td>
<td style="text-align:left;">
-20
</td>
<td style="text-align:left;">
-28
</td>
</tr>
<tr>
<td style="text-align:left;">
N
</td>
<td style="text-align:left;">
-80
</td>
<td style="text-align:left;">
-62
</td>
<td style="text-align:left;">
-46
</td>
<td style="text-align:left;">
-30
</td>
<td style="text-align:left;">
-25
</td>
<td style="text-align:left;">
-17
</td>
<td style="text-align:left;">
-6
</td>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
<span style="     color: red !important;">14</span>
</td>
<td style="text-align:left;">
6
</td>
<td style="text-align:left;">
-2
</td>
<td style="text-align:left;">
-10
</td>
<td style="text-align:left;">
-18
</td>
<td style="text-align:left;">
-15
</td>
<td style="text-align:left;">
-16
</td>
<td style="text-align:left;">
-21
</td>
</tr>
<tr>
<td style="text-align:left;">
C
</td>
<td style="text-align:left;">
-88
</td>
<td style="text-align:left;">
-70
</td>
<td style="text-align:left;">
-54
</td>
<td style="text-align:left;">
-38
</td>
<td style="text-align:left;">
-32
</td>
<td style="text-align:left;">
-25
</td>
<td style="text-align:left;">
-14
</td>
<td style="text-align:left;">
-3
</td>
<td style="text-align:left;">
6
</td>
<td style="text-align:left;">
<span style="     color: red !important;">27</span>
</td>
<td style="text-align:left;">
19
</td>
<td style="text-align:left;">
11
</td>
<td style="text-align:left;">
3
</td>
<td style="text-align:left;">
-5
</td>
<td style="text-align:left;">
-13
</td>
<td style="text-align:left;">
-17
</td>
</tr>
<tr>
<td style="text-align:left;">
I
</td>
<td style="text-align:left;">
-96
</td>
<td style="text-align:left;">
-78
</td>
<td style="text-align:left;">
-62
</td>
<td style="text-align:left;">
-46
</td>
<td style="text-align:left;">
-38
</td>
<td style="text-align:left;">
-33
</td>
<td style="text-align:left;">
-22
</td>
<td style="text-align:left;">
-11
</td>
<td style="text-align:left;">
-2
</td>
<td style="text-align:left;">
19
</td>
<td style="text-align:left;">
<span style="     color: red !important;">29</span>
</td>
<td style="text-align:left;">
21
</td>
<td style="text-align:left;">
13
</td>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
-3
</td>
<td style="text-align:left;">
-11
</td>
</tr>
<tr>
<td style="text-align:left;">
G
</td>
<td style="text-align:left;">
-104
</td>
<td style="text-align:left;">
-86
</td>
<td style="text-align:left;">
-70
</td>
<td style="text-align:left;">
-54
</td>
<td style="text-align:left;">
-46
</td>
<td style="text-align:left;">
-30
</td>
<td style="text-align:left;">
-30
</td>
<td style="text-align:left;">
-19
</td>
<td style="text-align:left;">
-10
</td>
<td style="text-align:left;">
11
</td>
<td style="text-align:left;">
21
</td>
<td style="text-align:left;">
<span style="     color: red !important;">37</span>
</td>
<td style="text-align:left;">
29
</td>
<td style="text-align:left;">
21
</td>
<td style="text-align:left;">
13
</td>
<td style="text-align:left;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
Q
</td>
<td style="text-align:left;">
-112
</td>
<td style="text-align:left;">
-94
</td>
<td style="text-align:left;">
-78
</td>
<td style="text-align:left;">
-62
</td>
<td style="text-align:left;">
-54
</td>
<td style="text-align:left;">
-38
</td>
<td style="text-align:left;">
-28
</td>
<td style="text-align:left;">
-27
</td>
<td style="text-align:left;">
-18
</td>
<td style="text-align:left;">
3
</td>
<td style="text-align:left;">
13
</td>
<td style="text-align:left;">
29
</td>
<td style="text-align:left;">
<span style="     color: red !important;">38</span>
</td>
<td style="text-align:left;">
30
</td>
<td style="text-align:left;">
22
</td>
<td style="text-align:left;">
14
</td>
</tr>
<tr>
<td style="text-align:left;">
K
</td>
<td style="text-align:left;">
-120
</td>
<td style="text-align:left;">
-102
</td>
<td style="text-align:left;">
-86
</td>
<td style="text-align:left;">
-70
</td>
<td style="text-align:left;">
-62
</td>
<td style="text-align:left;">
-46
</td>
<td style="text-align:left;">
-32
</td>
<td style="text-align:left;">
-25
</td>
<td style="text-align:left;">
-26
</td>
<td style="text-align:left;">
-5
</td>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
21
</td>
<td style="text-align:left;">
32
</td>
<td style="text-align:left;">
<span style="     color: red !important;">41</span>
</td>
<td style="text-align:left;">
<span style="     color: red !important;">33</span>
</td>
<td style="text-align:left;">
<span style="     color: red !important;">25</span>
</td>
</tr>
</tbody>
</table>

#### Alignment and score

``` r
printAlign(x=alignList)
```

    ## 
    ##  S1:  --PFGFGKRSCMGRRLA 
    ##         ||  | | | |     
    ##  S2:  FIPFSAGPRNCIGQK-- 
    ## 
    ##  Score of alignment: 25

## Local alignment

The alignment type choice is passed on to all following functions.

``` r
# align_type <- "global"
align_type <- "local"
```

### Dynamic programming matrices

``` r
dynMA <- dynProgMatrix(S1, S2, align_method=align_type, gap_penalty=8, substitutionMA="BLOSUM50")
```

The matrices are stored in a list and returned below.

``` r
dynMA
```

    ## $loc_ma
    ##    gp  P  F  G  F  G  K  R  S  C  M  G  R  R  L  A
    ## gp  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    ## F   0  0  8  0  8  0  0  0  0  0  0  0  0  0  1  0
    ## I   0  0  0  4  0  4  0  0  0  0  2  0  0  0  2  0
    ## P   0 10  2  0  0  0  3  0  0  0  0  0  0  0  0  1
    ## F   0  2 18 10  8  0  0  0  0  0  0  0  0  0  1  0
    ## S   0  0 10 18 10  8  0  0  5  0  0  0  0  0  0  2
    ## A   0  0  2 10 15 10  7  0  1  4  0  0  0  0  0  5
    ## G   0  0  0 10  7 23 15  7  0  0  1  8  0  0  0  0
    ## P   0 10  2  2  6 15 22 14  6  0  0  0  5  0  0  0
    ## R   0  2  7  0  0  7 18 29 21 13  5  0  7 12  4  0
    ## N   0  0  0  7  0  0 10 21 30 22 14  6  0  6  8  3
    ## C   0  0  0  0  5  0  2 13 22 43 35 27 19 11  4  7
    ## I   0  0  0  0  0  1  0  5 14 35 45 37 29 21 13  5
    ## G   0  0  0  8  0  8  0  0  6 27 37 53 45 37 29 21
    ## Q   0  0  0  0  4  0 10  2  0 19 29 45 54 46 38 30
    ## K   0  0  0  0  0  2  6 13  5 11 21 37 48 57 49 41
    ## 
    ## $loc_ma_path
    ##    gp P F G F G K R S C M G R R L A
    ## gp  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    ## F   0 1 1 3 1 3 1 1 1 1 1 1 1 1 1 1
    ## I   0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ## P   0 1 3 1 1 1 1 1 1 1 1 1 1 1 1 1
    ## F   0 2 1 3 1 3 1 1 1 1 1 1 1 1 1 1
    ## S   0 1 2 1 3 1 1 1 1 1 1 1 1 1 1 1
    ## A   0 1 2 1 1 1 1 3 1 1 1 1 1 1 1 1
    ## G   0 1 1 1 2 1 3 3 1 1 1 1 3 1 1 1
    ## P   0 1 3 2 1 2 1 3 1 3 1 2 1 1 1 1
    ## R   0 2 1 1 1 2 1 1 3 3 3 1 1 1 3 1
    ## N   0 1 2 1 3 1 2 2 1 3 3 3 1 1 1 1
    ## C   0 1 1 2 1 1 2 2 2 1 3 3 3 3 1 1
    ## I   0 1 1 1 1 1 1 2 2 2 1 3 3 3 1 3
    ## G   0 1 1 1 3 1 3 1 2 2 2 1 3 3 3 3
    ## Q   0 1 1 2 1 2 1 3 1 2 2 2 1 1 3 3
    ## K   0 1 1 1 1 1 1 1 3 2 2 2 1 1 3 3

### Compute alignment

Note: `alignList` stores all relevant results in a list, including dynamic
programming matrices, as well as the coordinates (named `path_coor`) to highlight path in
dynamic progamming matrix. This way on can easily show single dynamic programming matrix
with path highlighted in color or arrows in an HTML or PDF document.

``` r
alignList <- alignmentTraceback(ma=dynMA[[1]], ma_path=dynMA[[2]], align_method=align_type) 
names(alignList)
```

    ## [1] "ma"        "ma_path"   "path_coor" "as1"       "consensus" "as2"      
    ## [7] "score"

``` r
# alignList$ma # dyn ma with scores
# alignList$ma_path # dyn ma with path
# alignList$path_coor # coordinates for path to auto highlight path in HTML/PDF table
```

### Return results

#### Traceback in matrix

The following prints the fully populated dynamic programming matrix where the traceback
path is highlighted in color.

``` r
printColMa(alignList)
```

<table class=" lightable-paper lightable-striped" style="font-family: &quot;Arial Narrow&quot;, arial, helvetica, sans-serif; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:left;">
gp
</th>
<th style="text-align:left;">
P
</th>
<th style="text-align:left;">
F
</th>
<th style="text-align:left;">
G
</th>
<th style="text-align:left;">
F
</th>
<th style="text-align:left;">
G
</th>
<th style="text-align:left;">
K
</th>
<th style="text-align:left;">
R
</th>
<th style="text-align:left;">
S
</th>
<th style="text-align:left;">
C
</th>
<th style="text-align:left;">
M
</th>
<th style="text-align:left;">
G
</th>
<th style="text-align:left;">
R
</th>
<th style="text-align:left;">
R
</th>
<th style="text-align:left;">
L
</th>
<th style="text-align:left;">
A
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
gp
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
F
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
8
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
8
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
I
</td>
<td style="text-align:left;">
<span style="     color: red !important;">0</span>
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
4
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
4
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
P
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
<span style="     color: red !important;">10</span>
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
3
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
F
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
<span style="     color: red !important;">18</span>
</td>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
8
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
S
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
<span style="     color: red !important;">18</span>
</td>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
8
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
A
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
<span style="     color: red !important;">15</span>
</td>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
7
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
4
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
G
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
7
</td>
<td style="text-align:left;">
<span style="     color: red !important;">23</span>
</td>
<td style="text-align:left;">
15
</td>
<td style="text-align:left;">
7
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
8
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
P
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
6
</td>
<td style="text-align:left;">
15
</td>
<td style="text-align:left;">
<span style="     color: red !important;">22</span>
</td>
<td style="text-align:left;">
14
</td>
<td style="text-align:left;">
6
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
R
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
7
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
7
</td>
<td style="text-align:left;">
18
</td>
<td style="text-align:left;">
<span style="     color: red !important;">29</span>
</td>
<td style="text-align:left;">
21
</td>
<td style="text-align:left;">
13
</td>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
7
</td>
<td style="text-align:left;">
12
</td>
<td style="text-align:left;">
4
</td>
<td style="text-align:left;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
N
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
7
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
21
</td>
<td style="text-align:left;">
<span style="     color: red !important;">30</span>
</td>
<td style="text-align:left;">
22
</td>
<td style="text-align:left;">
14
</td>
<td style="text-align:left;">
6
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
6
</td>
<td style="text-align:left;">
8
</td>
<td style="text-align:left;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
C
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
13
</td>
<td style="text-align:left;">
22
</td>
<td style="text-align:left;">
<span style="     color: red !important;">43</span>
</td>
<td style="text-align:left;">
35
</td>
<td style="text-align:left;">
27
</td>
<td style="text-align:left;">
19
</td>
<td style="text-align:left;">
11
</td>
<td style="text-align:left;">
4
</td>
<td style="text-align:left;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
I
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
14
</td>
<td style="text-align:left;">
35
</td>
<td style="text-align:left;">
<span style="     color: red !important;">45</span>
</td>
<td style="text-align:left;">
37
</td>
<td style="text-align:left;">
29
</td>
<td style="text-align:left;">
21
</td>
<td style="text-align:left;">
13
</td>
<td style="text-align:left;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
G
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
8
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
8
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
6
</td>
<td style="text-align:left;">
27
</td>
<td style="text-align:left;">
37
</td>
<td style="text-align:left;">
<span style="     color: red !important;">53</span>
</td>
<td style="text-align:left;">
45
</td>
<td style="text-align:left;">
37
</td>
<td style="text-align:left;">
29
</td>
<td style="text-align:left;">
21
</td>
</tr>
<tr>
<td style="text-align:left;">
Q
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
4
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
19
</td>
<td style="text-align:left;">
29
</td>
<td style="text-align:left;">
45
</td>
<td style="text-align:left;">
<span style="     color: red !important;">54</span>
</td>
<td style="text-align:left;">
46
</td>
<td style="text-align:left;">
38
</td>
<td style="text-align:left;">
30
</td>
</tr>
<tr>
<td style="text-align:left;">
K
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
0
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
6
</td>
<td style="text-align:left;">
13
</td>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
11
</td>
<td style="text-align:left;">
21
</td>
<td style="text-align:left;">
37
</td>
<td style="text-align:left;">
48
</td>
<td style="text-align:left;">
<span style="     color: red !important;">57</span>
</td>
<td style="text-align:left;">
49
</td>
<td style="text-align:left;">
41
</td>
</tr>
</tbody>
</table>

#### Alignment and score

``` r
printAlign(x=alignList)
```

    ## 
    ##  S1:  PFGFGKRSCMGRR 
    ##       ||  | | | |   
    ##  S2:  PFSAGPRNCIGQK 
    ## 
    ##  Score of alignment: 57

# C. Different Substitution Matrices

Task 1: Load the Biostrings package in R, import the following two cytochrome
P450 sequences O15528 and P98187 from NCBI (save as myseq.fasta), and create a
global alignment with the pairwiseAlignment function from Biostrings as
follows.

``` r
library(Biostrings)
myseq <- readAAStringSet("myseq.fasta", "fasta")
(p <- pairwiseAlignment(myseq[[1]], myseq[[2]], type="global", substitutionMatrix="BLOSUM50"))
writePairwiseAlignments(p)
```

Your answers should address the following items:

Record the scores for the scoring matrices BLOSUM50, BLOSUM62 and BLOSUM80.
How and why do the scores differ for the three scoring matrices?

Answer 1: The scores for the three BLOSUM substitutions matrices are:

-   BLOSUM50: 227
-   BLOSUM62: 54
-   BLOSUM80: -52

Answer 2: Since the two sequences are relatively dissimilar (as determined by alignment view from
`writePairwiseAlignments(p)`) it is expected that the BLOSUM matrices trained on more dissimilar
sequences (*e.g.* BLSOSOM50) will result in higher scores than those trained on less similar
sequences (*e.g.* BLOSOM 80).

# Session Info

``` r
sessionInfo()
```

    ## R version 4.0.5 (2021-03-31)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Debian GNU/Linux 10 (buster)
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.8.0
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.8.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ## [1] kableExtra_1.3.4    Biostrings_2.58.0   XVector_0.30.0     
    ## [4] IRanges_2.24.0      S4Vectors_0.28.0    BiocGenerics_0.36.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bslib_0.2.4       compiler_4.0.5    jquerylib_0.1.3   tools_4.0.5      
    ##  [5] zlibbioc_1.36.0   digest_0.6.27     viridisLite_0.3.0 jsonlite_1.7.1   
    ##  [9] evaluate_0.14     lifecycle_0.2.0   rlang_0.4.8       rstudioapi_0.13  
    ## [13] yaml_2.2.1        blogdown_1.2      xfun_0.22         stringr_1.4.0    
    ## [17] httr_1.4.2        knitr_1.30        xml2_1.3.2        sass_0.3.1       
    ## [21] systemfonts_1.0.1 webshot_0.5.2     svglite_2.0.0     glue_1.4.2       
    ## [25] R6_2.5.0          rmarkdown_2.7     bookdown_0.21     magrittr_2.0.1   
    ## [29] scales_1.1.1      htmltools_0.5.1.1 rvest_0.3.6       colorspace_2.0-0 
    ## [33] stringi_1.5.3     munsell_0.5.0     crayon_1.3.4

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Needleman1970-md" class="csl-entry">

Needleman, S B, and C D Wunsch. 1970. “<span class="nocase">A general method applicable to the search for similarities in the amino acid sequence of two proteins</span>.” *J. Mol. Biol.* 48 (3): 443–53. <https://doi.org/10.1016/0022-2836(70)90057-4>.

</div>

<div id="ref-Smith1981-ax" class="csl-entry">

Smith, T F, and M S Waterman. 1981. “<span class="nocase">Identification of common molecular subsequences</span>.” *J. Mol. Biol.* 147 (1): 195–97. <http://www.ncbi.nlm.nih.gov/pubmed/7265238>.

</div>

</div>
