# This folder includes the R scripts for each sequential experiment  of this project.

There are **6 sections** which can be opened in [RStudio](https://rstudio.com/) 

Each script is formed by main components: fecundity analysis, cleaning up the raw datasets, recombination rate results, general statistics for each analysis. Cleanedup datasets can also be found in the ["raw_data_files"](https://github.com/StevisonLab/Peak-Plasticity-Project/tree/master/raw_data_files)

The analysis used the following packages. Make sure to read the corresponding documentation for each to understand what the code is doing.
- library([ggplot2](https://ggplot2.tidyverse.org/))
- library([ggthemes](https://cran.r-project.org/web/packages/ggthemes/ggthemes.pdf))
- library([emmeans](https://cran.r-project.org/web/packages/emmeans/index.html))
- library([lme4](https://cran.r-project.org/web/packages/lme4/lme4.pdf))
- library([lmerTest](https://cran.r-project.org/web/packages/lmerTest/index.html))
- library([doBy](https://cran.r-project.org/web/packages/doBy/doBy.pdf))
- library([reshape2](https://cran.r-project.org/web/packages/reshape2/reshape2.pdf))
- library([car](https://cran.r-project.org/web/packages/car/car.pdf))

Instructions to run code on your computer:
1. Once you download the github repository, navigate to the ["Scripts"](https://github.com/StevisonLab/Peak-Plasticity-Project/tree/master/Scripts) folder
2. Open R studio and set your working directory to the ["code"](https://github.com/StevisonLab/Peak-Plasticity-Project/tree/master/Scripts) folder
3. Raw data files are in the ["raw_data_files"](https://github.com/StevisonLab/Peak-Plasticity-Project/tree/master/raw_data_files) folder.
4. Open the R script in the same folder as the raw data files and run the code.
5. The script allows to analyze each experiment separately as well as sections in R can be run together.  
6. Reproducibility compares between Experiment 2 and 3 for sd-y, y-se, intervals in the same timepoint, 7-9 days post-mating. 
