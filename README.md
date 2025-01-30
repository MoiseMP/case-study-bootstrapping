# case-study-bootstrapping
Econometrics Case Study Period 3 2024-2025

## Tested on
R version 4.4.1 (2024-06-14)

## TODO: How to run

**Step 1: Setup**

*Option 1*: From command line
```{bash}
git clone https://github.com/username/case-study-bootstrapping.git
```
Open project in Rstudio

*Option 2*:
Inside of Rstudio -> File -> New Project -> Version Control -> Git -> Repository Url: https://github.com/MoiseMP/case-study-bootstrapping

**Step 2**: Install dependencies

```{R}
install.packages("renv")
renv::restore()
```
**Step 3**: Run project
Settings can be adjusted inside of script/analysis.R, file can be excecuted using
```{r}
source("script/analysis.R")
```

## Project Structure
Monte Carlo results are stored inside of output folder.